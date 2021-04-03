
// Build pSSAlib with MATLAB support:
// CC=$HOME/bin/gcc CXX=$HOME/bin/g++ CFLAGS="-fPIC" CXXFLAGS="-fPIC" ./configure --prefix=$HOME/.local --disable-sbml --disable-optimisation
// make clean; make -j && make install


// To compile:
// mex ./mex/mexpssa_lyapunov.cpp -I~/.local/include -L~/.local/lib64 -lboost_thread -lboost_chrono -lboost_system -lut -l:libpssa.a -l:libgsl.a -l:libgslcblas.a -v -DHAVE_CONFIG_H -g

/*
To test this code:

tc = 'ca';
k = [2 2.1 0.1 1.0 0.01 0.1 15 0 0];

% steady state solution
S_ss = zeros(2, 1);
S_ss(1) = (sqrt(k(5)^2 + 8 * k(2) * k(3) * k(6) / (k(4) + k(6))) - k(5))*(k(4) + k(6))/(4 * k(3) * k(6));
S_ss(2) = k(3) * S_ss(1)^2 / (k(4) + k(6));

% Jacobian
[-4*k(3)*S_ss(1)-k(5), 2*k(4); 2*k(3)*S_ss(1), -k(6)-k(4)]

% Rate functions
F = [k(2); k(3)*S_ss(1)^2; k(4)*S_ss(2); k(5)*S_ss(1); k(6)*S_ss(2)]

S_ss = S_ss * k(7);
[J, Q, S] = mexpssa_lyapunov(tc, S_ss, k)

C_anl = lyap(J, k(7) * S * diag(F) * S');
C_num = lyap(J, Q);

% Variance as estimated by LNA
S_ss .* diag(C)
*/

#include "mex.h"

#include "../include/models.h"
#include "../include/tools.h"

#include <boost/scoped_array.hpp>

void
computeQ(
  pssalib::datamodel::detail::Model & model,
  const double * pardPopulation,
  double * pardQ
)
{
  boost::scoped_array<double> ardF(new double[model.getReactionsCount()]);
  computeF(model, pardPopulation, ardF.get());

  for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
  {
    pssalib::datamodel::detail::Reaction * r = model.getReaction(ri);

    // Scale with compartment volume
    ardF[ri] *= model.getCompartmentVolume();

    for(UINTEGER sri = 0; sri < r->getSpeciesReferencesCount(); ++sri)
    {
      const pssalib::datamodel::detail::SpeciesReference * srI = r->getSpeciesReferenceAt(sri);

      for(UINTEGER srj = 0; srj < r->getSpeciesReferencesCount(); ++srj)
      {
        const pssalib::datamodel::detail::SpeciesReference * srJ = r->getSpeciesReferenceAt(srj);

        // fill the stoichiometry matrix S
        if(!srI->isReservoir()&&!srJ->isReservoir())
        {
          double dq = ardF[ri] * srI->getStoichiometryAbs() * srJ->getStoichiometryAbs();
          if((sri < r->getReactantsCount())^(srj < r->getReactantsCount()))
            pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] -= dq;
          else
            pardQ[srI->getIndex() + srJ->getIndex() * model.getSpeciesCount()] += dq;
        }
      }
    }
  }
}

void
mexFunction(
    int nlhs,
    mxArray *plhs[],
    int nrhs,
    const mxArray *prhs[])
{

  /* Check for proper number of input and output arguments */    
  if (nrhs < 3) {
      mexErrMsgIdAndTxt( "MATLAB:mexpssa:maxrhs",
              "Three input arguments required.");
  }
  if (nlhs < 2) {
      mexErrMsgIdAndTxt( "MATLAB:mexpssa:maxlhs",
              "Two output arguments are required.");
  }

  /*
    *Check if the input is of proper type
    */
  if(!mxIsChar(prhs[0])) { // not sa character string
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "First argument has to be a character string.");
  }
  size_t testCase = std::numeric_limits<size_t>::max();
  {
    char *pData = (char*)mxArrayToString(prhs[0]);
    if(NULL != pData)
    {
      if(0 == std::strcmp("clc",pData))
        testCase = 0;
      else if(0 == std::strcmp("ca",pData))
        testCase = 1;
      else if(0 == std::strcmp("homo",pData))
        testCase = 2;
      else if(0 == std::strcmp("sbd",pData))
        testCase = 3;
      else if(0 == std::strcmp("tcs",pData))
        testCase = 4;
      else if(0 == std::strcmp("ed",pData))
        testCase = 5;
      mxFree(pData);
    }
    if(5 < testCase) { // not a recognized character string
        mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                          "First argument has to be either 'ca' for Colloidal Aggregation, 'clc' for Cyclic Linear Chain, 'sbd' for SingleBirthDeath, 'homo' for Homoreaction, 'tcs' for TwoComponentSystem or 'ed' for Enzymatic Degradation.");
    }
  }
  /* species populations vector */
  if(!mxIsNumeric(prhs[1]) || // not numeric
      mxIsComplex(prhs[1])) { // or complex
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "Second argument has to be an integer vector with number of elements equal to number of species in the system.");
  }
  double * arPopulation = mxGetPr(prhs[1]);
  /* reaction rates vector */
  if(!mxIsNumeric(prhs[2]) || // not numeric
      mxIsComplex(prhs[2])) { // or complex
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "Third argument has to be a real vector.");
  }
  double * arParams = mxGetPr(prhs[2]);
  size_t szParams = mxGetNumberOfElements(prhs[2]);

  /* local variables */
  mxArray * mxarF, * mxarJ, * mxarStoichiometry, * mxarQ;
  pssalib::datamodel::detail::Model model;

  switch(testCase)
  {
  case 0:
    generateCyclicLinearChain(model, arParams, szParams);
    break;
  case 1:
    generateColloidalAggregation(model, arParams, szParams);
    break;
  case 2:
    generateHomoreaction(model, arParams, szParams);
    break;
  case 3:
    generateSingleBirthDeath(model, arParams, szParams);
    break;
  case 4:
    generateTwoComponentSystem(model, arParams, szParams);
    break;
  case 5:
    generateEnzymaticDegradation(model, arParams, szParams);
    break;
  default:
    mexErrMsgIdAndTxt("MATLAB:mexpssa:internal",
                      "Invalid test case identifier.");
  }

  // normalize the model
  model.normalize();

  // allocate arrays
  mxarQ = mxCreateDoubleMatrix(model.getSpeciesCount(), model.getSpeciesCount(), mxREAL);
  mxarJ = mxCreateDoubleMatrix(model.getSpeciesCount(), model.getSpeciesCount(), mxREAL);

  computeJ(model, arPopulation, mxGetPr(mxarJ));
  computeQ(model, arPopulation, mxGetPr(mxarQ));

  // output
  plhs[0] = mxarJ;
  plhs[1] = mxarQ;

  if(nlhs > 2)
  {
    mwSize nnz = 0, nnzi = 0;
    for(UINTEGER ri = 0; ri < model.getReactionsCount(); ++ri)
      nnz += model.getReaction(ri)->getSpeciesReferencesCount();

    mxarStoichiometry = mxCreateSparse(model.getSpeciesCount(), model.getReactionsCount(), nnz, mxREAL);
    generateS(model, arPopulation, &mxarStoichiometry);
    plhs[2] = mxarStoichiometry;

    if(nlhs > 3)
    {
      mxarF = mxCreateDoubleMatrix(model.getReactionsCount(), 1, mxREAL);
      computeF(model, arPopulation, mxGetPr(mxarF));
      plhs[3] = mxarF;
    }
  }
}

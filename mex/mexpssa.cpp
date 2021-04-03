
// Build pSSAlib with MATLAB support:
// CC=$HOME/bin/gcc CXX=$HOME/bin/g++ CFLAGS="-fPIC" CXXFLAGS="-fPIC" ./configure --prefix=$HOME/.local --disable-sbml --disable-optimisation
// make clean; make -j && make install


// To compile:
// mex ./mex/mexpssa.cpp -I~/.local/include -L~/.local/lib64 -lboost_thread -lboost_chrono -lboost_system -lut -l:libpssa.a -l:libgsl.a -l:libgslcblas.a -v -DHAVE_CONFIG_H -g

#include "mex.h"

#ifndef MX_API_VER
#  define mxIsScalar(x) (1 >= mxGetN(x)*mxGetM(x))
#endif

#include <time.h>       /* time_t, struct tm, time, localtime */

#include <pssalib/PSSA.h>
#include <pssalib/util/Timing.h>

#include "../include/models.h"

/* Global definitions for pSSAlib */
pssalib::datamodel::SimulationInfo simInfo;

#ifdef MULTITHREADED

#ifdef BOOST_NO_CXX11_HDR_THREAD
#  include <boost/thread/thread_only.hpp>
#  include <boost/chrono/time_point.hpp>
#else
#  include <thread>
#endif

/*
 * utIsInterruptPending(): "undocumented MATLAB API implemented in
 * libut.so, libut.dll, and included in the import library
 * libut.lib. To use utIsInterruptPending in a mex-file, one must
 * manually declare bool utIsInterruptPending() because this function
 * is not included in any header files shipped with MATLAB. Since
 * libut.lib, by default, is not linked by mex, one must explicitly
 * tell mex to use libut.lib." -- Wotao Yin, 
 * http://www.caam.rice.edu/~wy1/links/mex_ctrl_c_trick/
 *
 */
#ifdef __cplusplus
    extern "C" bool utIsInterruptPending();
    extern "C" bool utSetInterruptPending(bool);
#else
    extern     bool utIsInterruptPending();
    extern     bool utSetInterruptPending(bool);
#endif

// output variables
bool bProcessing, bResult;

void runSimulation(pssalib::datamodel::SimulationInfo *ptrSimInfo)
{
  pssalib::PSSA simEngine;

  simEngine.setMethod(pssalib::PSSA::M_PDM);

  bResult = simEngine.run(ptrSimInfo);

  bProcessing = false;
}
#endif /* HAVE_BOOST_THREAD */

void
mexFunction(
    int nlhs,
    mxArray *plhs[],
    int nrhs,
    const mxArray *prhs[])
{
  /* Check for proper number of input and output arguments */    
  if (nrhs < 2) {
      mexErrMsgIdAndTxt( "MATLAB:mexpssa:maxrhs",
              "Two input arguments required.");
  }
  if (nlhs > 2) {
      mexErrMsgIdAndTxt( "MATLAB:mexpssa:maxlhs",
              "Too many output arguments.");
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
  }
  if(5 < testCase) { // not a recognized character string
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "First argument has to be either 'ca' for Colloidal Aggregation, 'clc' for Cyclic Linear Chain, 'sbd' for SingleBirthDeath, 'homo' for Homoreaction, 'tcs' for TwoComponentSystem or 'ed' for Enzymatic Degradation.");
  }
  /* parameters vector */
  if(!mxIsNumeric(prhs[1]) || // not numeric
      mxIsComplex(prhs[1])) { // or complex
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "Second argument has to be a real vector.");
  }
  double * arParams = mxGetPr(prhs[1]);
  size_t szParams = mxGetNumberOfElements(prhs[1]);
  /* optional arguments */
  double timeStart = 0.0;
  double timeStep = 1e-1;
  double timeEnd = 1.0;
  if(nrhs > 2)
  {
    /* final time of the simulation */
    if(!mxIsNumeric(prhs[2]) || // not numeric
        mxIsComplex(prhs[2]) || // or complex
        !mxIsScalar(prhs[2]))  { // or not scalar
        mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                          "Third argument has to be a real scalar.");
    }
    timeEnd = mxGetScalar(prhs[2]);
    if(nrhs > 3)
    {
      /* time step for the simulation */
      if(!mxIsNumeric(prhs[3]) || // not numeric
          mxIsComplex(prhs[3]) || // or complex
          !mxIsScalar(prhs[3]))  { // or not scalar
          mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                            "Fourth argument has to be a real scalar.");
      }
      timeStep = mxGetScalar(prhs[3]);
      if(nrhs > 4)
      {
        /* time step for the simulation */
        if(!mxIsNumeric(prhs[4]) || // not numeric
            mxIsComplex(prhs[4]) || // or complex
            !mxIsScalar(prhs[4]))  { // or not scalar
            mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                              "Fifth argument has to be a real scalar.");
        }
        timeStart = mxGetScalar(prhs[4]);
      }
    }
  }

  pssalib::datamodel::detail::Model & model = simInfo.getModel();
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

  // setup the simulation params
  simInfo.unSamplesTotal = 1;
  simInfo.dTimeStart = timeStart;
  simInfo.dTimeEnd = timeEnd;
  simInfo.dTimeStep = timeStep;

  size_t szTimePoints = pssalib::timing::getNumTimePoints(simInfo.dTimeStart, simInfo.dTimeEnd, simInfo.dTimeStep);

  // storage for trajectories
  boost::scoped_ptr<UINTEGER> arunTrajectory(new UINTEGER[model.getSpeciesCount()*(szTimePoints)]);
#ifndef PSSALIB_ENGINE_CHECK
  simInfo.ptrarRawTrajectory = arunTrajectory.get();

  // setup output
  simInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
    | pssalib::datamodel::SimulationInfo::ofRawTrajectory
    | pssalib::datamodel::SimulationInfo::ofTrajectory
//       | pssalib::datamodel::SimulationInfo::ofStatus
    | pssalib::datamodel::SimulationInfo::ofLog
//       | pssalib::datamodel::SimulationInfo::ofTrace
    | pssalib::datamodel::SimulationInfo::ofInfo
//       | pssalib::datamodel::SimulationInfo::ofWarning
    | pssalib::datamodel::SimulationInfo::ofError
//       | pssalib::datamodel::SimulationInfo::eofModuleGrouping
//       | pssalib::datamodel::SimulationInfo::eofModuleSampling
//       | pssalib::datamodel::SimulationInfo::eofModuleUpdate
    ;
  simInfo.resetOutput();

  simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
//  simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cout.rdbuf());
#else

#warning "debug build"
  // setup output
  simInfo.unOutputFlags = pssalib::datamodel::SimulationInfo::ofNone
//       | pssalib::datamodel::SimulationInfo::ofRawTrajectory
    | pssalib::datamodel::SimulationInfo::ofTrajectory
//       | pssalib::datamodel::SimulationInfo::ofStatus
    | pssalib::datamodel::SimulationInfo::ofLog
    | pssalib::datamodel::SimulationInfo::ofTrace
    | pssalib::datamodel::SimulationInfo::ofInfo
    | pssalib::datamodel::SimulationInfo::ofWarning
    | pssalib::datamodel::SimulationInfo::ofError
    | pssalib::datamodel::SimulationInfo::eofModuleGrouping
    | pssalib::datamodel::SimulationInfo::eofModuleSampling
    | pssalib::datamodel::SimulationInfo::eofModuleUpdate
    ;
  simInfo.resetOutput();

  simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofLog, std::cerr.rdbuf());
  simInfo.setOutputStreamBuf(pssalib::datamodel::SimulationInfo::ofTrajectory, std::cout.rdbuf());
#endif

  {
    time_t rawtime;
    struct tm * timeinfo;

    time (&rawtime);
    timeinfo = localtime (&rawtime);
    std::cerr << "\n\t[MEXPSSA]: Simulation started at: " << asctime(timeinfo) << "\n";
  }

  double * dest = NULL;
#ifdef MULTITHREADED
  bProcessing = true;
#ifdef BOOST_NO_CXX11_HDR_THREAD
  boost::thread
#else
  std::thread
#endif
    simThread(&runSimulation, &simInfo);
      
  while(bProcessing)
  {
#ifdef BOOST_NO_CXX11_HDR_THREAD
    boost::this_thread::sleep_for(boost::chrono::seconds(1));
#else      
    std::this_thread::sleep_for(std::chrono::seconds(1));
#endif
    if(utIsInterruptPending()&&(!simInfo.bInterruptRequested))
    {
      std::cerr << "\n\t[MEXPSSA]: Simulation cancelled by keyboard interrupt!\n";
      utSetInterruptPending(false); // Consume the event
      simInfo.bInterruptRequested = true;
    }
  }

  if(bResult)
#else
  pssalib::PSSA simEngine;

  simEngine.setMethod(pssalib::PSSA::M_PDM);

  if(simEngine.run(&simInfo))
#endif
  {
    plhs[0] = mxCreateDoubleMatrix(szTimePoints, model.getSpeciesCount(), mxREAL);
    dest = mxGetPr(plhs[0]);
    UINTEGER * src = arunTrajectory.get();

    for(size_t i = 0; i < szTimePoints; ++i)
      for(size_t j = 0, in = i*model.getSpeciesCount(), jn = model.getSpeciesCount(); j < jn; ++j)
        *(dest+i+j*szTimePoints) = *(src+in+j);

    // produce time points
    if(nlhs > 1)
    {
      plhs[1] = mxCreateDoubleMatrix(szTimePoints, 1, mxREAL);
      double * timePoints = mxGetPr(plhs[1]);
      timePoints[0] = simInfo.dTimeStart;
      timePoints[szTimePoints-1] = simInfo.dTimeEnd;
      for(size_t i = 1; i < szTimePoints-1; ++i)
        timePoints[i] = timePoints[i-1] + simInfo.dTimeStep;
    }
  }
  else
  {
#ifdef MULTITHREADED
    if(simInfo.bInterruptRequested)
      mexErrMsgIdAndTxt("MATLAB:mexpssa:internal","Simulation cancelled by keyboard interrupt!");
    else
#endif /* MULTITHREADED */
      mexErrMsgIdAndTxt("MATLAB:mexpssa:internal","Simulation failed!");
  }

  // clean up
  simInfo.getModel().free();
}

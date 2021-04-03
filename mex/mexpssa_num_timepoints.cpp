
// Build pSSAlib with MATLAB support:
// CC=$HOME/bin/gcc CXX=$HOME/bin/g++ CFLAGS="-fPIC" CXXFLAGS="-fPIC" ./configure --prefix=$HOME/.local --disable-sbml --disable-optimisation
// make clean; make -j && make install


// To compile:
// mex ./mex/mexpssa_num_timepoints.cpp -I~/.local/include -L~/.local/lib64 -lboost_thread -lboost_chrono -lboost_system -lut -l:libpssa.a -l:libgsl.a -l:libgslcblas.a -v -DHAVE_CONFIG_H -g

#include "mex.h"

#include <pssalib/util/Timing.h>

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
  if (nlhs > 1) {
      mexErrMsgIdAndTxt( "MATLAB:mexpssa:maxlhs",
              "Too many output arguments.");
  }

  /* Check if the inputs are of proper types */
  /* final time of the simulation */
  if(!mxIsNumeric(prhs[0]) || // not numeric
      mxIsComplex(prhs[0]) || // or complex
     !mxIsScalar (prhs[0])) { // or not scalar
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "First argument has to be a real scalar.");
  }
  double timeEnd = mxGetScalar(prhs[0]);
  /* time step for the simulation */
  if(!mxIsNumeric(prhs[1]) || // not numeric
      mxIsComplex(prhs[1]) || // or complex
     !mxIsScalar (prhs[1])) { // or not scalar
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "Second argument has to be a real scalar.");
  }
  double timeStep = mxGetScalar(prhs[1]);
  /* time step for the simulation */
  if(!mxIsNumeric(prhs[2]) || // not numeric
      mxIsComplex(prhs[2]) || // or complex
     !mxIsScalar (prhs[2])) { // or not scalar
      mexErrMsgIdAndTxt("MATLAB:mexpssa:typeargin",
                        "Third argument has to be a real scalar.");
  }
  double timeStart = 0.0;
  if (nrhs > 2) {
    timeStart = mxGetScalar(prhs[2]);
  }

  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double * dest = mxGetPr(plhs[0]);
  *dest = pssalib::timing::getNumTimePoints(timeStart, timeEnd, timeStep);
}

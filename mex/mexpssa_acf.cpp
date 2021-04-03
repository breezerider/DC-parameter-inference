
// To compile:
// mex ./mex/mexpssa_acf.cpp -v -g

#include "mex.h"
#include <cmath>
#include <iostream>

#ifndef MX_API_VER
#  define mxIsScalar(x) (1 >= mxGetN(x)*mxGetM(x))
#endif

void 
mexFunction(
    int nlhs,
    mxArray *plhs[],
    int nrhs,
    const mxArray *prhs[])
{
    /* Check for proper number of input and output arguments */    
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:mexpssa_acf:maxrhs",
                "Two input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:mexpssa_acf:maxlhs",
                "Too many output arguments.");
    }
    /* 
     *Check if the input is of proper type
     */
    /* input dataset */
    if(!mxIsNumeric(prhs[0]) || // not numeric
        mxIsComplex(prhs[0]) || // or complex
        mxIsScalar(prhs[0])) {  // or scalar
        mexErrMsgIdAndTxt("MATLAB:mexpssa_acf:typeargin",
                          "First argument has to be integer scalar.");
    }
    double * arData = mxGetPr(prhs[0]);
    mwSize szDataRows = mxGetM(prhs[0]),
           szDataCols = mxGetN(prhs[0]);
    /*  */
    if(!mxIsNumeric(prhs[1]) || // not numeric
        mxIsComplex(prhs[1]) || // or complex
       ((mxGetM(prhs[1]) > 1)&&
       (mxGetN(prhs[1]) > 1))){ // or wrong size
        mexErrMsgIdAndTxt("MATLAB:mexpssa_acf:typeargin",
                          "Second argument has to be an integer row or column vector.");
    }
    double * arLag = mxGetPr(prhs[1]);
    mwSize szLagLen = mxGetM(prhs[1]) * mxGetN(prhs[1]);

    /* allocate output buffer */
    plhs[0] = mxCreateDoubleMatrix(szLagLen, szDataCols, mxREAL);
    double * arDest = mxGetPr(plhs[0]);

    /* calculate the autocorrelation function */
    for(mwSize lagIdx = 0; lagIdx < szLagLen; ++lagIdx)
    {
      size_t lag = std::floor(arLag[lagIdx]);

      for(mwSize colIdx = 0; colIdx < szDataCols; ++colIdx)
      {
        double mu = 0.0, sigma = 0.0, acf = 0.0;

        if(0 == lag)
        {
          arDest[lagIdx+colIdx*szLagLen] = 1.0;
          continue;
        }

        for(mwSize rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
          mu += arData[rowIdx+colIdx*szDataRows];
        mu /= (double)szDataRows;

        for(mwSize rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
          sigma += std::pow(arData[rowIdx+colIdx*szDataRows] - mu, 2.0);
//         sigma = sigma / (double)szDataRows;

        if(0.0 != sigma)
        {
          for(mwSize rowIdx = 0; rowIdx < (szDataRows - lag); ++rowIdx)
            acf += (arData[rowIdx+colIdx*szDataRows] - mu) *
              (arData[rowIdx+lag+colIdx*szDataRows] - mu);
          acf /= sigma;
        }
        else
          acf = 0.0;

        arDest[lagIdx+colIdx*szLagLen] = acf;
      }
    }

}

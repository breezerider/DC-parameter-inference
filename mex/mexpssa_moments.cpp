
// To compile:
// mex ./mex/mexpssa_moments.cpp -v -g

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
        mexErrMsgIdAndTxt( "MATLAB:mexpssa_moments:maxrhs",
                "Two input arguments required.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:mexpssa_moments:maxlhs",
                "Too many output arguments.");
    }
    /* 
     *Check if the input is of proper type
     */
    /* input dataset */
    if(!mxIsNumeric(prhs[0]) || // not numeric
        mxIsComplex(prhs[0]) || // or complex
        mxIsScalar(prhs[0])) {  // or scalar
        mexErrMsgIdAndTxt("MATLAB:mexpssa_moments:typeargin",
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
        mexErrMsgIdAndTxt("MATLAB:mexpssa_moments:typeargin",
                          "Second argument has to be an integer row or column vector.");
    }
    double * arOrder = mxGetPr(prhs[1]);
    mwSize szOrderLen = mxGetM(prhs[1]) * mxGetN(prhs[1]);

    /* allocate output buffer */
    plhs[0] = mxCreateDoubleMatrix(szOrderLen, szDataCols, mxREAL);
    double * arDest = mxGetPr(plhs[0]);

    /* calculate the moments */
    for(mwSize orderIdx = 0; orderIdx < szOrderLen; ++orderIdx)
    {
      size_t order = std::floor(arOrder[orderIdx]);

      for(mwSize colIdx = 0; colIdx < szDataCols; ++colIdx)
      {
        double mu = 0.0, sigma = 0.0, mu_i = 0.0;

        if(0 == order)
        {
          arDest[orderIdx+colIdx*szOrderLen] = 0.0;
          continue;
        }

        for(mwSize rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
          mu += arData[rowIdx+colIdx*szDataRows];
        mu /= (double)szDataRows;

        if(1 < order)
        {
          for(mwSize rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
            sigma += std::pow(arData[rowIdx+colIdx*szDataRows] - mu, 2.0);
          sigma = std::sqrt(sigma / (double)szDataRows);

          if(2 < order)
          {
            if(sigma > 0.0)
            {
              for(mwSize rowIdx = 0; rowIdx < szDataRows; ++rowIdx)
                mu_i += std::pow(arData[rowIdx+colIdx*szDataRows] - mu, (double)order);
              mu_i = std::pow(std::abs(mu_i / (double)szDataRows), 1.0 / ((double)order)) / sigma;
            }
            else
              mu_i = 0.0;
          }
          else
            mu_i = sigma;
        }
        else
          mu_i = mu;

        arDest[orderIdx+colIdx*szOrderLen] = mu_i;
      }
    }

}

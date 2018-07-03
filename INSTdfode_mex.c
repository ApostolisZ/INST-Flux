/*
 * kfpdxode_mex.c
 * Author: Desmond Lun (dslun@rutgers.edu)
 */

#include "mex.h"
#include "eval_kfpsol.h"

void kfpdfode(double t, double *dX, size_t dX_len, mxArray *dY, double *A, 
        double *B, double *C_diag, double *ddXdt)
{
    int i, j, i1;
    double *mid = mxMalloc(dX_len * sizeof(double));
    size_t mid_len;
    double *Z = mxCalloc(dX_len, sizeof(double));
    
    /* B * dY */
    for (j = 0; j < mxGetM(dY); j++) {
        mid_len = 0;
        
        for (i = 0; i < dX_len; i++) {
            if (B[j * dX_len + i] != 0) {
                if (mid_len == 0)
                    mid_len = eval_kfpsol(mxGetCell(dY, j), t, dX_len, mid, 1);
                
                for (i1 = 0; i1 < mid_len; i1++) {
                    /* mexPrintf("mid[%d] = %f\n", i1, mid[i1]); */
                    Z[i + i1] += B[j * dX_len + i] * mid[i1];
                }
            }
        }
    }
    
    /*
    for (i = 0; i < dX_len; i++)
        mexPrintf("Z[%d] = %f\n", i, Z[i]);
     */
    
    for (i = 0; i < dX_len; i++) {
        ddXdt[i] = Z[i];
        for (j = 0; j < dX_len; j++) {
            ddXdt[i] += A[j * dX_len + i] * dX[j];
        }
        ddXdt[i] /= C_diag[i];
    }
    
    mxFree(mid);
    mxFree(Z);
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double t;
    double *dX;      /* input column vector */
    double *ddXdt;   /* output column vector */
    mxArray *dY;    
    double *A, *B, *C_diag;
    size_t dX_len;
    
    t = mxGetScalar(prhs[0]);
    dX = mxGetPr(prhs[1]);
    dX_len = mxGetM(prhs[1]);
    dY = mxGetField(prhs[2], 0, "dY");
    A = mxGetPr(mxGetField(prhs[2], 0, "A"));
    B = mxGetPr(mxGetField(prhs[2], 0, "B"));
    C_diag = mxGetPr(mxGetField(prhs[2], 0, "C_diag"));
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) dX_len, 1, mxREAL);
    ddXdt = mxGetPr(plhs[0]);
    
    kfpdfode(t, dX, dX_len, dY, A, B, C_diag, ddXdt);
}
/*
 * kfpode_mex.c
 * Author: Desmond Lun (dslun@rutgers.edu)
 */

#include "mex.h"
#include "eval_kfpsol.h"

void kfpode(double t, double *X, size_t X_len, mxArray *Y,
        double *A, double *B, double *C_diag, double *dXdt)
{
    int i, j, i1;
    double *mid = mxMalloc(X_len * sizeof(double));
    size_t mid_len;
    double *Z = mxCalloc(X_len, sizeof(double));
    
    for (j = 0; j < mxGetM(Y); j++) {
        mid_len = 0;
        
        for (i = 0; i < X_len; i++) {
            if (B[j * X_len + i] != 0) {
                if (mid_len == 0)
                    mid_len = eval_kfpsol(mxGetCell(Y, j), t, X_len, mid, 0);
                
                for (i1 = 0; i1 < mid_len; i1++) {
                    /* mexPrintf("mid[%d] = %f\n", i1, mid[i1]); */
                    Z[i + i1] += B[j * X_len + i] * mid[i1];
                }
            }
        }
    }
    
    for (i = 0; i < X_len; i++) {
        dXdt[i] = Z[i];
        for (j = 0; j < X_len; j++) {
            dXdt[i] += A[j * X_len + i] * X[j];
        }
        dXdt[i] /= C_diag[i];
    }
    
    mxFree(mid);
    mxFree(Z);
}

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    double t;
    double *X;      /* input column vector */
    double *dXdt;   /* output column vector */
    mxArray *Y;     /* cell array Y */
    double *A, *B, *C_diag;
    size_t X_len;
    
    t = mxGetScalar(prhs[0]);
    X = mxGetPr(prhs[1]);
    X_len = mxGetM(prhs[1]);
    Y = mxGetField(prhs[2], 0, "Y");
    A = mxGetPr(mxGetField(prhs[2], 0, "A"));
    B = mxGetPr(mxGetField(prhs[2], 0, "B"));
    C_diag = mxGetPr(mxGetField(prhs[2], 0, "C_diag"));
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) X_len, 1, mxREAL);
    dXdt = mxGetPr(plhs[0]);
    
    kfpode(t, X, X_len, Y, A, B, C_diag, dXdt);
}

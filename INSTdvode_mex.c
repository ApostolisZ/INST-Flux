/*
 * kfpdvode_mex.c
 * Author: Desmond Lun (dslun@rutgers.edu)
 */

#include "mex.h"
#include "eval_kfpsol.h"

void kfpdvode(double t, double *dX, size_t dX_len, mxArray *X, mxArray *Y, 
        mxArray *dY, double *A, double *B, double *C_diag, double *dA, 
        double *dB, int *A_inds, double *ddXdt)
{
    int i, j, i1;
    double *mid = mxMalloc(dX_len * sizeof(double));
    size_t mid_len;
    double *Z = mxCalloc(dX_len, sizeof(double));
    
    /* B * dY/dv */
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
    
    /* dA/dV * X */
    for (j = 0; j < mxGetM(X); j++) {
        mid_len = 0;
        
        for (i = 0; i < dX_len; i++) {
            if (dA[A_inds[j] * dX_len + i] != 0) {
                if (mid_len == 0)
                    mid_len = eval_kfpsol(mxGetCell(X, j), t, dX_len, mid, 0);
                
                for (i1 = 0; i1 < mid_len; i1++) {
                    /* mexPrintf("mid[%d] = %f\n", i1, mid[i1]); */
                    Z[i + i1] += dA[A_inds[j] * dX_len + i] * mid[i1];
                }
            }
        }
    }
    
    /*
    for (i = 0; i < dX_len; i++)
        mexPrintf("Z[%d] = %f\n", i, Z[i]);
     */
    
    /* dB/dV * Y */
    for (j = 0; j < mxGetM(Y); j++) {
        mid_len = 0;
        
        for (i = 0; i < dX_len; i++) {
            if (dB[j * dX_len + i] != 0) {
                if (mid_len == 0)
                    mid_len = eval_kfpsol(mxGetCell(Y, j), t, dX_len, mid, 0);
                
                for (i1 = 0; i1 < mid_len; i1++) {
                    /* mexPrintf("mid[%d] = %f\n", i1, mid[i1]); */
                    Z[i + i1] += dB[j * dX_len + i] * mid[i1];
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
    mxArray *X, *Y, *dY;    
    double *A, *B, *dA, *dB, *C_diag;
    double *emu_natoms;
    int *A_inds;
    size_t dX_len;
    int i, nemu;
    
    t = mxGetScalar(prhs[0]);
    dX = mxGetPr(prhs[1]);
    dX_len = mxGetM(prhs[1]);
    X = mxGetField(prhs[2], 0, "X");    
    Y = mxGetField(prhs[2], 0, "Y");
    dY = mxGetField(prhs[2], 0, "dY");
    A = mxGetPr(mxGetField(prhs[2], 0, "A"));
    B = mxGetPr(mxGetField(prhs[2], 0, "B"));
    dA = mxGetPr(mxGetField(prhs[2], 0, "dA"));
    dB = mxGetPr(mxGetField(prhs[2], 0, "dB"));
    C_diag = mxGetPr(mxGetField(prhs[2], 0, "C_diag"));
    emu_natoms = mxGetPr(mxGetField(prhs[2], 0, "emu_natoms"));
    nemu = mxGetM(mxGetField(prhs[2], 0, "emu_natoms"));
    
    A_inds = mxMalloc(nemu * sizeof(int));
    A_inds[0] = 0;
    for (i = 1; i < nemu; i++)
        A_inds[i] = A_inds[i - 1] + (int) emu_natoms[i - 1];
    
    plhs[0] = mxCreateDoubleMatrix((mwSize) dX_len, 1, mxREAL);
    ddXdt = mxGetPr(plhs[0]);
    
    kfpdvode(t, dX, dX_len, X, Y, dY, A, B, C_diag, dA, dB, A_inds, ddXdt);
}
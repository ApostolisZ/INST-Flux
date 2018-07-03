/*
 * eval_kfpsol.c
 * Author: Desmond Lun (dslun@rutgers.edu)
 */

#include "mex.h"

void convolve(const double a[], size_t a_len,
        const double b[], size_t b_len,
        double c[])
{
    size_t n;
    
    for (n = 0; n < a_len + b_len - 1; n++)
    {
        size_t kmin, kmax, k;
        
        c[n] = 0;
        
        kmin = (n >= b_len - 1) ? n - (b_len - 1) : 0;
        kmax = (n < a_len - 1) ? n : a_len - 1;
        
        for (k = kmin; k <= kmax; k++)
        {
            c[n] += a[k] * b[n - k];
        }
    }
}

void deval23t(double *t, double *y, mxArray *idata, size_t t_len,
        size_t m, double tint, double *inds, size_t inds_len, 
        double *sxint)
{
    int bottom, i;
    double s, s2, s3, v1, v2;
    double *z = mxGetPr(mxGetField(idata, 0, "z"));
    double *znew = mxGetPr(mxGetField(idata, 0, "znew"));
    
    for (bottom = 0; bottom < t_len - 1 && t[bottom + 1] < tint; bottom++) {}
        
    if (t[bottom] == tint) {
        for (i = 0; i < inds_len; i++)
            sxint[i] = y[bottom * m + (size_t) inds[i] - 1];
    }
    else if (bottom < t_len - 1) {
        s = (tint - t[bottom]) / (t[bottom + 1] - t[bottom]);
        s2 = s * s;
        s3 = s * s2;
        for (i = 0; i < inds_len; i++) {            
            v1 = y[(bottom + 1) * m + (size_t) inds[i] - 1] -
                    y[bottom * m + (size_t) inds[i] - 1] -
                    z[(bottom + 1) * m + (size_t) inds[i] - 1];
            v2 = znew[(bottom + 1) * m + (size_t) inds[i] - 1] -
                    z[(bottom + 1) * m + (size_t) inds[i] - 1];
            sxint[i] = y[bottom * m + (size_t) inds[i] - 1] +
                    z[(bottom + 1) * m + (size_t) inds[i] - 1] * s +
                    (3 * v1 - v2) * s2 + (v2 - 2 * v1) * s3;
        }
    }
    else
        mexErrMsgIdAndTxt("kfpode_mex", "tint out of range.");
}

void deval23tb(double *t, double *y, mxArray *idata, size_t t_len,
        size_t m, double tint, double *inds, size_t inds_len, 
        double *sxint)
{
    int bottom, i;
    double a1, a2, a3;
    double *t2 = mxGetPr(mxGetField(idata, 0, "t2"));            
    double *y2 = mxGetPr(mxGetField(idata, 0, "y2"));
    
    for (bottom = 0; bottom < t_len - 1 && t[bottom + 1] < tint; bottom++) {}
    
    if (t[bottom] == tint) {
        for (i = 0; i < inds_len; i++)
            sxint[i] = y[bottom * m + (size_t) inds[i] - 1];
    }
    else if (bottom < t_len - 1) {
        a1 = ((tint - t[bottom + 1]) * (tint - t2[bottom + 1])) /
                ((t[bottom] - t[bottom + 1]) * (t[bottom] - t2[bottom + 1]));
        a2 = ((tint - t[bottom]) * (tint - t2[bottom + 1])) /
                ((t[bottom + 1] - t[bottom]) * (t[bottom + 1] - t2[bottom + 1]));
        a3 = ((tint - t[bottom]) * (tint - t[bottom + 1])) /
                ((t2[bottom + 1] - t[bottom]) * (t2[bottom + 1] - t[bottom + 1]));
        for (i = 0; i < inds_len; i++) {
            sxint[i] = a1 * y[bottom * m + (size_t) inds[i] - 1] + 
                    a2 * y[(bottom + 1) * m + (size_t) inds[i] - 1] + 
                    a3 * y2[(bottom + 1) * m + (size_t) inds[i] - 1];
        }
    }
    else
        mexErrMsgIdAndTxt("kfpode_mex", "tint out of range.");
}

size_t eval_kfpsol(mxArray *kfpsol, double t, size_t max_len, double *mid,
        int isderiv)
{
    mxArray *de_sol = mxGetCell(kfpsol, 0);
    mxArray *inds = mxGetCell(kfpsol, 1);
    int i, j, i1;
    size_t m = mxGetM(de_sol);
    size_t n = mxGetN(de_sol);
    double *tmp_mid = mxMalloc((max_len + 1) * sizeof(double));
    double *tmp_mid2 = mxMalloc((max_len + 1) * sizeof(double));
    double *sxint = mxMalloc((max_len + 1) * sizeof(double));
    mwIndex index;
    mxArray *idata, *x, *y, *inds_array;
    size_t tmp_mid_len, inds_array_len;
    
    for (i = 0; i < max_len; i++)
        mid[i] = 0.0;
    
    for (i = 0; i < m; i++) {
        tmp_mid_len = 1;
        tmp_mid[0] = 1.0;
        for (j = 0; j < n; j++) {
            index = (mwIndex) (j * m + i);
            idata = mxGetField(de_sol, index, "idata");
            x = mxGetField(de_sol, index, "x");
            y = mxGetField(de_sol, index, "y");
            inds_array = mxGetCell(inds, index);
            inds_array_len = mxGetN(inds_array);
            
            deval23t(mxGetPr(x), mxGetPr(y), idata, mxGetN(x), mxGetM(y),
                    t, mxGetPr(inds_array), inds_array_len, sxint + 1);
            
            if (isderiv && j == 0)
                sxint[0] = 0.0;
            else
                sxint[0] = 1.0;
            for (i1 = 0; i1 < inds_array_len; i1++)
                sxint[0] -= sxint[1 + i1];

            /*
            for (i1 = 0; i1 < inds_array_len + 1; i1++)
                mexPrintf("sxint[%d] = %f\n", i1, sxint[i1]);
             */
            convolve(sxint, inds_array_len + 1, tmp_mid, tmp_mid_len,
                    tmp_mid2);            
            
            tmp_mid_len += inds_array_len;            
            for (i1 = 0; i1 < tmp_mid_len; i1++) {
                tmp_mid[i1] = tmp_mid2[i1];
            }
        }
        
        for (i1 = 0; i1 < tmp_mid_len - 1; i1++)
            mid[i1] += tmp_mid[i1 + 1];
    }
    
    mxFree(tmp_mid);
    mxFree(tmp_mid2);
    mxFree(sxint);
    
    return tmp_mid_len - 1;
}
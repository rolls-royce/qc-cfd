/*************************************************************************************************/
/*                                                                                               */
/* Copyright (c) 2022 Rolls-Royce plc                                                            */
/*                                                                                               */
/* Redistribution and use in source and binary forms, with or without modification, are          */
/* permitted provided that the following conditions are met:                                     */
/*                                                                                               */
/* 1. Redistributions of source code must retain the above copyright notice, this list of        */
/*    conditions and the following disclaimer.                                                   */
/* 2. Redistributions in binary form must reproduce the above copyright notice, this list of     */
/*    conditions and the following disclaimer in the documentation and/or other materials        */
/*    provided with the distribution.                                                            */
/* 3. Neither the name of the copyright holder nor the names of its contributors may be used to  */
/*    endorse or promote products derived from this software without specific prior written      */
/*    permission.                                                                                */
/*                                                                                               */
/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS   */
/* OR IMPLIED  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF              */
/* MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE    */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE */
/* GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED    */
/* AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING     */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED  */
/* OF THE POSSIBILITY OF SUCH DAMAGE.                                                            */
/*                                                                                               */
/*************************************************************************************************/
/**
 * \file    matsol.c
 * \brief   Routines for solving matrix equations Ax = b
 * \details These are locally coded solvers to avoid the need to link an
 *          external library such as GSL (Gnu Scientific Library).
 *          The routines convert sparse matrices to full matrices which
 *          should be fine for the 1D nozzle: a full coupled matrix for
 *          a compressible nozzle with 64 stations requires 300KB of
 *          memory.
 *          If memory is an issue, the Python interface can be used which
 *          invokes scipy.sparse.linalg.sparse to solve the matrix equations.
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

/******************************************************************************
 *  Include files                                                             *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>               // defines DBL_MAX

#include <nozzle.h>
#include <proto.h>

/******************************************************************************
 *  Inline function for addressing matrix entries                             *
 ******************************************************************************/
static inline int mij(int i,int j, int n) {return n*i + j;}

/******************************************************************************
 *  matsol_LU - solve matrix equation using LU solver                         *
 ******************************************************************************/
/**
 * \brief   Solve matrix equation using LU solver
 * \details Converts incoming sparse matrix to a full matrix for LU decomposition.
 *          Requires storage for 3 full matrices < 1MB for a compressible nozzle 
 *          with 64 stations.
 *          Uses Doolittle algorithm for LU decomposition. Assumes all principal
 *          minors of A are non-singular.
 *
 * \param[in]
 *          A         \ref SMATRIX structure holding matrix and right hand
 *                    side vector, b.
 * \param[out]
 *          x         solution Ax = b, memory allocated here and must be 
 *                    managed by the calling routine. The size of x is set as
 *                    number of columns in A.
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */


int  matsol_LU(SMATRIX *A, double **x)
{
     char routine[] = "matsol_LU";

     int32_t nr = A->nr, nc = A->nc;
     int32_t *col = A->col;
     int32_t *rst = A->rst;
     double  *val = A->val;
     double  *b   = A->b;

     double *M, *L, *U, *y, *xp, *z, sum, min, max;
     int     i, j, k, n, n2;
     bool    checkLU = false, checkXY = false;

/*   check matrix is square
     ---------------------- */
     if(nr != nc){
       error_message(routine, __FILE__, __LINE__, "Matrix S[%d, %d] is not square", nr, nc);
       return 1;
     }
     n = nr;
     n2 = nr*nc;

/*   copy sparse matrix to full matrix
     --------------------------------- */
     M = (double*) malloc(n2 * sizeof(double));
     if(!M){
       error_message(routine, __FILE__, __LINE__, 
                     "Sparse to full memory allocation error for S[%d, %d]", nr, nc);
       return 2;
     }
     memset(M, 0, n2 * sizeof(double));

     for(i=0; i<nr; i++)
     {
       for(k=rst[i]; k<rst[i+1]; k++)
       {
         j = col[k];
         M[mij(i,j,n)] = val[k];
       }
     }

/*   allocate full memory for L and U matrices
     ----------------------------------------- */
     L = (double*) malloc(n2 * sizeof(double));
     U = (double*) malloc(n2 * sizeof(double));
     if(!L || !U){
       error_message(routine, __FILE__, __LINE__, 
                     "LU memory allocation error for S[%d, %d]", nr, nc);
       free(M);
       if(L) free(L);
       if(U) free(U);
       return 3;
     }
     memset(L, 0, n2 * sizeof(double));
     memset(U, 0, n2 * sizeof(double));

/*   perform LU (Doolittle Algorithm)
     -------------------------------- */
     for(i=0; i<n; i++)
     {
       for(k=i; k<n; k++)
       {
         sum = 0.0;
         for (j=0; j<i; j++) sum += L[mij(i,j,n)] * U[mij(j,k,n)];
         U[mij(i,k,n)] = M[mij(i,k,n)] - sum;
       }
 
       for (k=i; k<n; k++)
       {
         if (i == k)
           L[mij(i,i,n)] = 1;
         else
         {
           sum = 0.0;
           for (j=0; j<i; j++) sum += L[mij(k,j,n)] * U[mij(j,i,n)];
           L[mij(k,i,n)] = (M[mij(k,i,n)] - sum)/U[mij(i,i,n)];
         }
       }
     }

/*   check LU = A, actually LU - M = 0
     --------------------------------- */
     if(checkLU){
       min =  DBL_MAX;
       max = -DBL_MAX;

       for(i=0; i<n; i++) {
         for(j=0; j<n; j++) {
           M[mij(i,j,n)] = -M[mij(i,j,n)];
           for(k=0; k<n; k++) {
             M[mij(i,j,n)] += L[mij(i,k,n)] * U[mij(k,j,n)];
           }
           min = MIN(min, M[mij(i,j,n)]);
           max = MAX(max, M[mij(i,j,n)]);
         }
       }
       printf("LU - A: min_ij = % le,\t max_ij =  % le\n", min, max);
     }

/*   free M - ensures memory for x and y
     ----------------------------------- */
     free(M);

/*   solve Ly = b
     ------------ */
     y = (double*) malloc(n * sizeof(double));
     y[0] = b[0]/L[mij(0,0,n)];
     for(i=1; i<n; i++)
     {
       y[i] = b[i];
       for(j=0; j<i; j++) y[i] -= L[mij(i,j,n)] * y[j];
       y[i] /= L[mij(i,i,n)];
     }

/*   solve Ux = y
     ------------ */
     *x = (double*) malloc(n * sizeof(double));
     xp = *x;

     xp[n-1] = y[n-1]/U[mij(n-1,n-1,n)];
     for(i=n-2; i>=0; i--)
     {
       xp[i] = y[i];
       for(j=i+1; j<n; j++) xp[i] -= U[mij(i,j,n)] * xp[j];
       xp[i] /= U[mij(i,i,n)];
     }

/*   check Ly=b and Ux = y
     --------------------- */
     if(checkXY){
       z = (double*) malloc(n * sizeof(double));

       min =  DBL_MAX;
       max = -DBL_MAX;
       for(i=0; i<n; i++) {
         z[i] = 0;
         for(j=0; j<n; j++) {
           z[i] += L[mij(i,j,n)] * y[j];
         }
         min = MIN(min, b[i] - z[i]);
         max = MAX(max, b[i] - z[i]);
       }
       printf("b - Ly: min_ij = % le,\t max_ij =  % le\n", min, max);

       min =  DBL_MAX;
       max = -DBL_MAX;
       for(i=0; i<n; i++) {
         z[i] = 0;
         for(j=0; j<n; j++) {
           z[i] += U[mij(i,j,n)] * xp[j];
         }
         min = MIN(min, y[i] - z[i]);
         max = MAX(max, y[i] - z[i]);
       }
       printf("y - Ux: min_ij = % le,\t max_ij =  % le\n", min, max);

       free(z);
     }

/*   debug
     ----- */
//   for(i=0; i<n; i++)printf("b[%3d] = % le\tx[%3d] = % le\n", i, b[i], i, xp[i]);

/*   free local arrays
     ----------------- */
     free(L);
     free(U);
     free(y);

     return 0;
}


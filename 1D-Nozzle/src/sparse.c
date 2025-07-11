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
 * \file   sparse.c
 * \brief  Routines for sparse matrix operations
 * \author Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date   20th Dec 2022
 */

/******************************************************************************
 *  Include files                                                             *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <nozzle.h>
#include <proto.h>

/******************************************************************************
 *  sparse_free - free a sparse matrix                                        *
 ******************************************************************************/
/**
 * \brief   Free a sparse matrix
 * \details Free a sparse matrix
 *
 * \param[in]
 *          S         \ref SMATRIX structure holding matrix
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    20th Dec 2022
 */

void sparse_free(SMATRIX *S)
{
     if(S->val) free(S->val);
     if(S->col) free(S->col);
     if(S->rst) free(S->rst);
     if(S->b)   free(S->b);

     S->val = NULL;
     S->col = NULL;
     S->rst = NULL;
     S->b   = NULL;
     S->nr  = 0;
     S->nc  = 0;
     S->nnz = 0;
}

/******************************************************************************
 *  sparse_print - print a sparse matrix                                      *
 ******************************************************************************/
/**
 * \brief   Print a sparse matrix
 * \details Print a sparse matrix using precsion based on number of columns
 *
 * \param[in]
 *          name      name to be printed with the matrix, NULL if no name
 *                     required
 * \param[in]
 *          S         \ref SMATRIX structure holding matrix
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    20th Dec 2022
 */

void sparse_print(char *name, SMATRIX S)
{
     int32_t nr = S.nr, nc = S.nc;
     int32_t *col = S.col;
     int32_t *rst = S.rst;
     double  *val = S.val;
     int32_t i, j, jc, n;
     int     sf;

/*   set precision
     ------------- */
     if      (nc <= 4)  sf = 6;
     else if (nc <= 8)  sf = 4;
     else if (nc <= 16) sf = 3;
     else               sf = 2;

/*   name
     ---- */
     if(name) printf("%s\n", name);

/*   column indices
     -------------- */
     printf("    ");
     for(j=0; j<nc; j++) printf(" % *d", sf+4, j);
     printf("\n");

/*   matrix including row indices
     ---------------------------- */
     for(i=0; i<nr; i++)
     {
       printf("%4d",i);
       jc = 0;
       for(n=rst[i]; n<rst[i+1]; n++)
       {
         while(jc != col[n] && jc < nc)     // fill gaps with zeros
         {
           printf("  % .*lf", sf, 0.0);
           jc += 1;
         }
         printf("  % .*lf", sf, val[n]);
         jc += 1;
       }

       for(j=jc; j<nc; j++) printf("  % .*lf", sf, 0.0);     // fill to end of row
       printf("\n");
     }

}

/******************************************************************************
 *  sparse_matvec - post multiply a sparse matrix with a vector               *
 ******************************************************************************/
/**
 * \brief   Calculate y = Sx where S is a sparse matrix and x,y are vectors
 * \details Calling routine must ensure vectors have the correct size.
 *
 * \param[in]
 *          S         \ref SMATRIX structure holding matrix
 * \param[in]
 *          x         vector to be multiplied by S, must have the same size as the
 *                    number of columns in S
 * \param[in]
 *          y         result y = Sx, memory for b must be allocated by the
 *                    calling routine. This must have the same size as the
 *                    number of rows in S
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    29th Dec 2022
 */

void sparse_matvec(SMATRIX *S, double *x, double *y)
{    
     int32_t nr = S->nr, nc = S->nc;
     int32_t *col = S->col;
     int32_t *rst = S->rst;
     double  *val = S->val;
     int32_t i, j, n;

/*   perform multiplication
     ---------------------- */
     memset(y, 0, nr*sizeof(double));

     for(i=0; i<nr; i++)
     {
       for(n=rst[i]; n<rst[i+1]; n++)
       {
         j = col[n];
         y[i] += val[n] * x[j];
       }
     }
}

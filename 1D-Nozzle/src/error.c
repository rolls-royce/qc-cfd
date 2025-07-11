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

#include <string.h>
#include <stdarg.h>

#include <proto.h>

/******************************************************************************
 *  nan_check - check double vector for a NaN                                 *
 ******************************************************************************/
/**
 * \brief   Check for a NaN in double vector
 * \details A NaN is detected when d != d.
 *
 * \param[in]
 *          name      name of vector
 * \param[in]
 *          n         length of vector
 * \param[in]
 *          d         vector 
 *
 * \return  error indicator for traceback
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    15th Nov 2022
 */

int  nan_check(const char *name, int n, double *d)
{
     int i;

     for(i=0; i<n; i++){
       if(d[i] != d[i]){
         error_message("nan_check", __FILE__, __LINE__,
                       "NaN found in array %s[%ld]\n", name, i);
         return 1;
       }
     }
     return 0;
}

/******************************************************************************
 *  error_message - print error message with context                          *
 ******************************************************************************/
/**
 * \brief   Print error message with context
 * \details Includes routine name, filename and line number with error message
 *
 * \param[in]
 *          routine name of routine, this should be explicitly set by the calling 
 *                  module
 * \param[in]
 *          file    name of file, calling routine should use macro __FILE__
 * \param[in]
 *          line    line number, calling routine should use macro __LINE__
 * \param[in]
 *          fmt     format statement for va_list
 * \param[in]
 *          ...     list of variables for va_list
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
 */

void error_message(const char *routine, const char *file, const int line, char *fmt, ...)
{
     const int maxlen = 104;
     char      message[maxlen];
     int       msglen;
     bool      truncate = false;
     va_list   list;

/*   get message from va_list and check for truncation
     ------------------------------------------------- */
     va_start(list, fmt);
     msglen = vsnprintf(message, maxlen, fmt, list);
     if(msglen > maxlen) truncate = true;
     va_end(list);

/*   if message truncated, replace last 3 characters with ellipsis
     ------------------------------------------------------------- */
     if(truncate){
       message[maxlen-4] = '\0';
       strncat(message, "...", 3);
       message[maxlen-1] = '\0';
     }

/*   print error message
     ------------------- */
     printf("Error in %s (%s:%d)\n\t%s\n", routine, file, line, message);
}

/******************************************************************************
 *  error_traceback - print traceback context                                 *
 ******************************************************************************/
/**
 * \brief   Print traceback context
 * \details Includes routine name, filename and line number with error message
 *
 * \param[in]
 *          routine name of routine, this should be explicitly set by the calling 
 *                  module
 * \param[in]
 *          file    name of file, calling routine should use macro __FILE__
 * \param[in]
 *          line    line number, calling routine should use macro __LINE__
 *
 * \author  Leigh Lapworth <leigh.lapworth@rolls-royce.com>
 * \date    1st Dec 2022
 */

void error_traceback(const char *routine, const char *file, const int line)
{

/*   print error traceback message
     ----------------------------- */
     printf("Traceback error in %s (%s:%d)\n", routine, file, line);
}

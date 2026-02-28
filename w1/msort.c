/*  File: msort.c
 *  Author: Jean Thierry-Mieg (mieg@ncbi.nlm.nih.gov)
 *  Copyright (C) J Thierry-Mieg, 2002
 * -------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 * -------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@sanger.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *      This file implements the classic sort/merge algorithm.
 *      
 *      This file is adapted from the GNU msort 
 *      written by Mike Haertel, September 1988.
 *      I got from him the idea of copying only
 *      (n - n2) bytes back in b, and added the
 *      already sorted case short cut.
 *      The rest of the code is standard.
 *	
 * Exported functions:
 *       mSort () with the same prototype as unix qsort()
 * HISTORY:
 * Last edited:
 * Created: Aug 6 2002 (mieg)
 *-------------------------------------------------------------------
 */

/* $Id: msort.c,v 1.9 2014/09/08 04:23:24 mieg Exp $ */

#include "regular.h"


/* time gain is less than 1/200 on 32bytes struct  */
static void insertionSort (char *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb))
{
  mysize_t i, j ;
  char buf[s] ;
  for (i = 1 ; i < n ; i++)
    {
      j = i - 1 ;
      if ((*cmp) (b + i*s, b + j*s) >= 0)
	continue ;
      memcpy (buf, b + i*s, s) ;
      memcpy (b + i*s, b + j*s, s) ;
      while (j > 0 && (*cmp) (b + (j-1)*s, buf) > 0)
	{
	  memcpy (b + j*s, b + (j-1)*s, s) ;
	  j-- ;
	}
      memcpy (b + j*s, buf, s) ;
    }
} /* insertionSort */


static void insertionSortInt (int *b, mysize_t n, int (*cmp)(const void *va, const void *vb))
{
   mysize_t i, j ;
   int tmp ;
   for (i = 1 ; i < n ; i++)
    {
      j = i - 1 ;
      if (cmp (b + i, b+j) >= 0)
	continue ;
      tmp = b[i] ;
      b[i] = b[j] ;
      while (j > 0 && cmp (b + j - 1, &tmp) > 0)
	{
	  b[j] = b[j-1] ;
	  j-- ;
	}
      b[j] = tmp ;
    }
} /* insertionSortInt */

static void mDoSortInt (int *b, mysize_t n, int (*cmp)(const void *va, const void *vb), int *buf)
{
  int *b0 = buf, *b1, *b2 ;
  mysize_t  n1, n2 ;

  /* 2025_05_24: 32 is the sweet spot butthe difference is minimal */
  if (n <= 32) { insertionSortInt(b, n, cmp) ; return ; }
  
  n1 = n / 2 ; n2 = n - n1 ;
  b1 = b ; b2 = b + n1 ;

  if (n1 > 1) mDoSortInt (b1, n1, cmp, buf);
  if (n2 > 1) mDoSortInt (b2, n2, cmp, buf);

  if ((*cmp) (b2 - 1, b2) <= 0) /* already sorted case */
    return ;

  /* 2010_09_05
   * this could be optimized (here and identically in the next function)
   * as long as b1<b2, we do not need to copy to buf and back again to b
   * vice versa as soon as b1 is done, we do not need to copy any more
   * the current code does a little of that but not completely 
   */

  while (n1 > 0 && n2 > 0)
    if ((*cmp) (b1, b2) <= 0)
      { n1-- ; *b0++  = *b1++ ; }
    else
      { n2-- ; *b0++  = *b2++ ; }

  if (n1 > 0)
    memcpy (b0, b1, n1 * sizeof (int)) ;
  memcpy (b, buf, (n - n2) * sizeof (int)) ;
} /* mDoSortInt */

/***********************************************************/

static void mDoSortVoid (void **b, mysize_t n, int (*cmp)(const void *va, const void *vb), void **buf)
{
  void **b0 = buf, **b1, **b2 ;
  mysize_t  n1, n2 ;

  n1 = n / 2 ; n2 = n - n1 ;
  b1 = b ; b2 = b + n1 ;

  if (n1 > 1) mDoSortVoid (b1, n1, cmp, buf);
  if (n2 > 1) mDoSortVoid (b2, n2, cmp, buf);

  if ((*cmp) (b2 - 1, b2) <= 0) /* already sorted case */
    return ;

  /* 2010_09_05
   * this could be optimized (here and identically in the next function)
   * as long as b1<b2, we do not need to copy to buf and back again to b
   * vice versa as soon as b1 is done, we do not need to copy any more
   * the current code does a little of that but not completely 
   */

  while (n1 > 0 && n2 > 0)
    if ((*cmp) (b1, b2) <= 0)
      { n1-- ; *b0++  = *b1++ ; }
    else
      { n2-- ; *b0++  = *b2++ ; }

  if (n1 > 0)
    memcpy (b0, b1, n1 * sizeof (void*)) ;
  memcpy (b, buf, (n - n2) * sizeof (void*)) ;
} /* mDoSortVoid */

/***********************************************************/

static void mDoSort (char *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb), char *buf)
{
  char *b0 = buf, *b1, *b2 ;
  mysize_t  n1, n2 ;

  /* gain is around 2% on 32bytes struct  */
  if (0 && n <= 16) { insertionSort (b, n, s, cmp) ; return ; }
  
  n1 = n >> 1 ; n2 = n - n1 ;
  b1 = b ; b2 = b + (n1 * s);

  if (n1 > 1) mDoSort (b1, n1, s, cmp, buf);
  if (n2 > 1) mDoSort (b2, n2, s, cmp, buf);

  if ((*cmp) (b2 - s, b2) <= 0) /* already sorted case */
    return ;

  while (n1 > 0 && n2 > 0)
    { 
      if ((*cmp) (b1, b2) <= 0)
	{
	  memcpy (b0, b1, s) ;
	  n1-- ; b1 += s ; b0 += s ;
	}
      else
	{
	  memcpy (b0, b2, s) ;
	  n2-- ; b2 += s ; b0 += s ;
	}
    }

  if (n1 > 0)
    memcpy (b0, b1, n1 * s) ;
  memcpy (b, buf, (n - n2) * s) ;
} /* mDoSort */

/***********************************************************/

void mSort (void *b, mysize_t n, int s, int (*cmp)(const void *va, const void *vb))
{ 
  mysize_t size = n * s ;
  int s_of_i = sizeof (int) ;
  int s_of_v = sizeof (void*) ;
  char  sBuf [4*1024] ; 
  char *mBuf = 0, *buf ;
  /* char bufBig [size] ;   2022_01_28, always alloc big on the stack */

  if (n < 2)
    return ;

  /* If size is small, use the stack, else malloc */
  /* optimize ? may be we should use malloc rather than messalloc
   * which initializes to zero
   * also we should have a size_of_void* dedicated code
   */

  if (size < 4*1024)
    buf = sBuf ;
  else
    buf = mBuf = (char *) malloc (size) ;
  
    /* buf = bufBig ; */

  if (s == s_of_i &&   /* aligned integers.  Use direct int pointers.  */
      (((char *)b - (char *) 0) % s_of_i) == 0)
    mDoSortInt ((int*)b, n, cmp, (int*)buf) ;
  else if (s == s_of_v &&   /* aligned void*.  Use direct word pointers.  */
      (((char *)b - (char *) 0) % s_of_v) == 0)
    mDoSortVoid ((void**)b, n, cmp, (void**)buf) ;
  else
    mDoSort ((char *)b, n, s, cmp, buf) ;

  free (mBuf) ;
} /* mSort */

/***********************************************************/
/***********************************************************/


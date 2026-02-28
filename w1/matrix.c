/*  File: matrix.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmba.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg yand R Durbin, 1991
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@crbm.cnrs-mop.fr
 *
 * Description:
 *              Extension of Array, to define typed matrix and tensors
 *              Bool, Int, Float, Complex
 *              Bool is the same as Int modulo 2
 *
 *    matrices are used in particular in the Deep learning library
 *
 * Exported functions:
 *              See Header file: matrix.h
 * Public read-only subfields (can only be altered at creation time)
 *    rank : 1 = vector, 2 = matrix, 3... = tensor
 *    shape : int[], a list of rank positive integers 
 *    mxAdd, mxTranspose, mxDot, mxPointMultiply
 *    
 * HISTORY:
 * Last edited:
 * Created: Feb 3, 2018
 *-------------------------------------------------------------------
 */

#define MAXRANK 64

#include "ac.h"
#include <complex.h>

/* SIYNCHRONIZE VERY CAREFULLY WITH matrix.h matrix.c matrix_slices.c */
/* Type of a matrix */
typedef enum {
  MX_NULL = 0, MX_BOOL, MX_INT, MX_FLOAT, MX_COMPLEX
} MX_TYPE ;

static int MX_MAGIC = 7861663 ;
typedef struct mxStruct *MX ;
typedef struct mxStruct { 
  const char *name ;
  MX_TYPE type ;
  int rank ;
  int shapes[MAXRANK] ;
  int slice0, slice1 ;

  void *magic ;
  complex float *zc ;
  float *zf ;
  int *zi ;
  int size ;
  int typeSize ;
  int deltas[MAXRANK] ; 
} MM ;
static const char *typeName[] = { "NULL", "BOOL", "INT", "FLOAT", "COMPLEX", "overflow"} ;
#define MATRIX_IS_DEFINED
#include "matrix.h"

/**********************************************************************/

BOOL mxCheck (MX a, const char *nam)
{
  if (!a)
    messcrash ("%s received a null matrix", nam) ;
  if (a->magic != &MX_MAGIC)
    messcrash ("%s received a bad pointer matrix", nam) ;
  
  return TRUE ;
} /* mxCheck */

/**********************************************************************/

static BOOL mxCheckSameRank (MX a, MX b, const char *nam)
{
  if (a->rank != b->rank)
    messcrash ("Rank mismatch %d != %d in %s :: %s %s"
	       , a->rank, b->rank, nam, a->name, b->name) ;
  return TRUE ;
} /* mxCheckSameRank */

/**********************************************************************/

static BOOL mxCheckRankN (MX a, int n, const char *nam)
{
  if (a->rank != n)
    messcrash ("Rank mismatch %d != %d in %s :: %s"
	       , a->rank, n, nam, a->name) ;
  return TRUE ;
} /* mxCheckRankN */

/**********************************************************************/

static BOOL mxCheckSameType (MX a, MX b, const char *nam)
{
  if (a->type != b->type)
    messcrash ("Matrix type mismatch in %s :: ((%s)->type = %s) != ((%s)->type = %s)"
	       ,nam, a->name, typeName[a->type], b->name, typeName[b->type]
	       ) ;
  return TRUE ;
} /* mxCheckSameType */

/**********************************************************************/

static BOOL mxCheckSameShapes (MX a, MX b, const char *nam)
{
  int i ;
  mxCheckSameRank (a, b, nam) ;
  for (i = 0 ; i < a->rank ; i++)
    if (a->shapes[i] != b->shapes[i])
      messcrash ("Shapes[%d] mismatch %d != %d in %s :: %s %s"
		 , i, a->shapes[i], b->shapes[i], nam, a->name, b->name) ;
  
  return TRUE ;
} /* mxCheckSameShape */


void uMxDestroy (void *va)
{
  MX a = (MX) va ;
  if (a)
    {
      if (a->magic !=  (&MX_MAGIC)) 
	messcrash ("mxDestroy received corrupt array->magic");	
      a->magic = 0 ;
      messfree (a->zi) ;
      messfree (a->zf) ;
      messfree (a->zc) ;
     }
  return ; 
} /* uMxDestroy  */

/**********************************************************************/
/* Matrix creation */
MX mxCreate (AC_HANDLE h,const char *name, MX_TYPE type, ...)
{
  MX a = (MX) handleAlloc (uMxDestroy, h, sizeof (MM)) ; /* a is zero initialized */
  int ii, delta ;
  int shapes[MAXRANK] ;
  va_list args ;
  
  a->magic = &MX_MAGIC ;
  a->type = type ;
  a->rank = 0 ;
  a->name = name ? strnew (name, h) : "no_name" ;

  va_start (args, type) ;
  ii = 0 ;

  while (1)
    {
      int k = va_arg(args, int );  

      if (k == -999)
	{
	  int *s = va_arg(args, int*) ;
	  for (ii = 0 ; ii< MAXRANK ; ii++)
	    {
	      if (s[ii] > 0)
		{
		  a->rank++ ;
		  shapes[ii] = s[ii] ;
		}
	      else if (s[ii] == 0)
		break ;
	      else
		messcrash ("mxCreate %s, called with a negative dim[%d] = %d < 0 in shapes[]", a->name, ii, s[ii]) ;
	    }
	  break ;
	}
      if (k < 0)
	  messcrash ("mxCreate %s, called with a dim %d = %d < 0", a->name, a->rank, k) ;
      if (k == 0)
	break ;
      a->rank++ ;
      if (a->rank > MAXRANK) 
	messcrash ("mxCreate %s, called with rank = %d > %d", a->name, a->rank, MAXRANK) ;
      shapes[ii++] = k ;
    }
  va_end(args);

  /* we start with the fastest moving index */
  for (ii = 0 ; ii < a->rank ; ii++)
    a->shapes[ii] = shapes[ii] ;

  a->size = 1 ;
  for (ii = 0 ; ii < MAXRANK ; ii++)
    if (a->shapes[ii] < 0)
      a->shapes[ii] = 0 ;
  for (ii = 0 ; ii < a->rank ; ii++)
    {
      a->size *= a->shapes[ii] ;
      if (ii == 0)
	delta = 1 ;
      else
	delta *=  a->shapes[ii -1] ;
      a->deltas[ii] = delta ;

      if (a->size > (1 << 30))
	messcrash ("mxCreate %s, rank %d called with more than 1G cells", a->name, a->rank) ;
    }
  
  switch (type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      a->typeSize = sizeof (int) ;
      a->zi = halloc (a->typeSize * a->size, 0) ;
       break ;
    case MX_FLOAT:
      a->typeSize = sizeof (float) ;
      a->zf = halloc (a->typeSize * a->size, 0) ;
       break ;
    case MX_COMPLEX:
      a->typeSize = sizeof (complex float) ;
      a->zc = halloc (a->typeSize * a->size, 0) ;
      break ;
    }

  return a ;
} /* mxCreate */

/**********************************************************************/

MX mxCopy (MX a, MX b, AC_HANDLE h)
{
  mxCheck (b, "mxCopy(...,b,...)") ;
  if (!a)
    {
      int i, kk[MAXRANK] ;
      int rank = b->rank ;

      for (i = 0 ; i < rank ; i++)
	kk[i] = b->shapes[i] ;
      kk[i] = 0 ;
      
      a = mxCreate (h
		    , b->name
		    , b->type 
		    , kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6], kk[7], kk[8], kk[9], kk[10], kk[11], kk[12], kk[13], kk[14], kk[15], kk[16], kk[17], kk[18], kk[19], kk[20], kk[21], kk[22], kk[23], kk[24], kk[25], kk[26], kk[27], kk[28], kk[29], kk[30], kk[31], kk[32], kk[33], kk[34], kk[35], kk[36], kk[37], kk[38], kk[39], kk[40], kk[41], kk[42], kk[43], kk[44], kk[45], kk[46], kk[47], kk[48], kk[49], kk[50], kk[51], kk[52], kk[53], kk[54], kk[55], kk[56], kk[57], kk[58], kk[59], kk[60], kk[61], kk[62], kk[63]
		    ) ;
    }
  mxCheckSameRank (a, b, "mxCopy") ;
  mxCheckSameType (a, b, "mxCopy") ;
  mxCheckSameShapes (a, b, "mxCopy") ;

  switch (b->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      a->typeSize = sizeof (int) ;
      a->zi = halloc (a->typeSize * a->size, 0) ;
      memcpy (a->zi, b->zi, a->typeSize * a->size)  ;
      break ;
    case MX_FLOAT:
      a->typeSize = sizeof (float) ;
      a->zf = halloc (a->typeSize * a->size, 0) ;
      memcpy (a->zf, b->zf, a->typeSize * a->size)  ;
      break ;
    case MX_COMPLEX:
      a->typeSize = sizeof (complex float) ;
      a->zc = halloc (a->typeSize * a->size, 0) ;
      memcpy (a->zc, b->zc, a->typeSize * a->size)  ;
      break ;
    }
  return a ;
} /* mxCopy */

/**********************************************************************/
/* create an empty A matrix with same type and shapes as the template B */ 
static MX mxDuplicate (const char *name, MX a, MX b, MX c, AC_HANDLE h)
{
  mxCheck (b, "mxDuplicate(...,b,...)") ;
  if (!a)
    {
      int i, kk[MAXRANK] ;
      int rank = b->rank ;
      MX_TYPE type, tb, tc ;
      tb = b->type ;
      tc = c ? c->type : 0 ;
      type = tb >= tc ? tb : tc ;

      for (i = 0 ; i < rank ; i++)
	kk[i] = b->shapes[i] ;
      kk[i] = 0 ;
      
      a = mxCreate (h
		    , name
		    , type 
		    , kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6], kk[7], kk[8], kk[9], kk[10], kk[11], kk[12], kk[13], kk[14], kk[15], kk[16], kk[17], kk[18], kk[19], kk[20], kk[21], kk[22], kk[23], kk[24], kk[25], kk[26], kk[27], kk[28], kk[29], kk[30], kk[31], kk[32], kk[33], kk[34], kk[35], kk[36], kk[37], kk[38], kk[39], kk[40], kk[41], kk[42], kk[43], kk[44], kk[45], kk[46], kk[47], kk[48], kk[49], kk[50], kk[51], kk[52], kk[53], kk[54], kk[55], kk[56], kk[57], kk[58], kk[59], kk[60], kk[61], kk[62], kk[63]
		    ) ;

    }
 return a ;
} /* mxDuplicate */

/**********************************************************************/

MX mxScalarMatrix (const char *name, float x, AC_HANDLE h)
{
  MX Un = mxCreate (h, name, MX_FLOAT, 1, 0) ;
  float un[] = { x } ;
  mxSet (Un, un) ;
  return Un ;
} /* mxScalarMatrix  */

/**********************************************************************/

void  mxValues (MX a, const int **zip, const float **zfp, const complex float **zcp) 
{
  mxCheck (a, "mxValues (a)") ;
  if (zip) *zip = a->zi ;
  if (zfp) *zfp = a->zf ;
  if (zcp) *zcp = a->zc ;
  
  return ;
} /* mxValues */

/**********************************************************************/
/* transform a matrix of size 0 into a complex */
float complex mxTopCorner (MX a)
{
  float complex zc = 0 ;
  mxCheck (a, "mxTopCorner (a)") ;

  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      zc = a->zi[0] ;
      break ;
    case MX_FLOAT:
      zc = a->zf[0] ;
      break ;
    case MX_COMPLEX:
      zc = a->zf[0] ;
      break ;
    default:
      break ;
    }
  return zc ;
} /* mxTopCorner */


/**********************************************************************/

MX mxAdd (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxLinearCombine (a, 1.0, b, 1.0, c, h) ;
}

/**********************************************************************/

MX mxMultiSum (AC_HANDLE h, MX *mms)
{
  int ii ;
  MX m = mms[0] ;

  if (! m)
    return 0 ;
  for (ii = 1 ; mms[ii] ; ii++)
    m = mxAdd (m, m, mms[ii], h) ;

  return m ;
} /* mxMultiSum */

/**********************************************************************/

MX mxSubstract (MX b, MX c, AC_HANDLE h)
{
  MX a = 0 ;
  return mxLinearCombine (a, 1.0, b, - 1.0, c, h) ;
}

/**********************************************************************/
/**********************************************************************/
/* Matrix multiplication */
/*
  a multiplicatiom where we contract the firts n indices
  leads to a very simple block repetitive algorithm
*/

static MX mxContractFirstOrLastIndex (MX a, MX b, MX c, int nFirst, int nLast, AC_HANDLE h)
{
  int ii, kb, kc, kbMax = 0, kcMax = 0, g, gMax = 0, dgb = 0, dgc = 0 ;

  int zi, * Restrict zia, * Restrict zib, * Restrict zic ;
  float zf, * Restrict zfa, * Restrict zfb, * Restrict zfc ;
  complex float zc, * Restrict zca, * Restrict zcb, * Restrict zcc ;

  if (nFirst * nLast > 0 || nFirst + nLast == 0)
    messcrash ("mxContractFirstOrLastIndex nFirst = %d > AND nLast = %d > 0, only one of them should be non null"
	       , nFirst, nLast
	       ) ;
  
  mxCheck (b, " mxContractOrLastIndex(...,b,...)") ;
  mxCheck (c, " mxContractFirstOrLastIndex(...,c)") ;
  if (nFirst > b->rank)
    messcrash ("mxContractFirstOrLastIndex nFirst = %d > b->rank = %d ::  b = %s"
	       , nFirst, b->rank, b->name
	       ) ;
  if (nFirst > c->rank)
    messcrash ("mxContractFirstOrLastIndex nFirst = %d > c->rank = %d ::  c = %s"
	       , nFirst, c->rank, c->name
	       ) ;
  if (nLast > b->rank)
    messcrash ("mxContractFirstOrLastIndex nLast = %d > b->rank = %d ::  b = %s"
	       , nLast, b->rank, b->name
	       ) ;
  if (nLast > c->rank)
    messcrash ("mxContractFirstOrLastIndex nLast = %d > c->rank = %d ::  c = %s"
	       , nLast, c->rank, c->name
	       ) ;

  for (ii = 0 ; ii < nFirst ; ii++)
    if (b->shapes[ii] != c->shapes[ii])
    messcrash ("mxContractFirstNIndices (N=%d) %s %s, shape mismatch, the first %d indices of the 2 matrices should have the same range:  (b->shapes[%d]) = %d != (c->shapes[%d] = %d)"
	       , nFirst
	       , b->name
	       , c->name
	       , nFirst
	       , ii, b->shapes[ii]
	       , ii, c->shapes[ii]
	       ) ;

  for (ii = 0 ; ii < nLast ; ii++)
    if (b->shapes[b->rank - ii - 1] != c->shapes[c->rank - ii - 1])
      messcrash ("mxContractLasttIndex %s %s, shape mismatch, the last index of the 2 matrices should have the same range:  (a->shapes[last]) = %d != (b->shapes[last] = %d)"
		 , b->name
		 , c->name
		 , b->shapes[b->rank - 1]
		 , c->shapes[c->rank - 1]
		 ) ;
  if (!a)
    {
      int i, j ;
      int kk[MAXRANK] ;
      MX_TYPE type = b->type > c->type ? b->type : c->type ;

      memset (kk, 0, sizeof (kk)) ;
      
      if (nFirst)
	{
	  for (i = nFirst, j = 0 ; i < b->rank ; i++)
	    kk[j++] = b->shapes[i] ;
	  for (i = nFirst ; i < c->rank ; i++)
	    kk[j++] = c->shapes[i] ;
	}
      else if (nLast)
	{
	  for (i = 0, j = 0 ; i < b->rank - nLast ; i++)
	    kk[j++] = b->shapes[i] ;
	  for (i = 0 ; i < c->rank - nLast ; i++)
	    kk[j++] = c->shapes[i] ;
	}
	
      a = mxCreate (h
		    , hprintf (h, "%s * %s", b->name, c->name)
		    , type 
		    , kk[0], kk[1], kk[2], kk[3], kk[4], kk[5], kk[6], kk[7], kk[8], kk[9], kk[10], kk[11], kk[12], kk[13], kk[14], kk[15], kk[16], kk[17], kk[18], kk[19], kk[20], kk[21], kk[22], kk[23], kk[24], kk[25], kk[26], kk[27], kk[28], kk[29], kk[30], kk[31], kk[32], kk[33], kk[34], kk[35], kk[36], kk[37], kk[38], kk[39], kk[40], kk[41], kk[42], kk[43], kk[44], kk[45], kk[46], kk[47], kk[48], kk[49], kk[50], kk[51], kk[52], kk[53], kk[54], kk[55], kk[56], kk[57], kk[58], kk[59], kk[60], kk[61], kk[62], kk[63]
		    ) ;
    }
  else
    mxClear (a) ;

  mxCheck (a, " mxContractFirstOrLastIndex (a,...)") ;
  
  if (a == b)
      messcrash ("mxContractFirstOrLastIndex A = <B | C>, A and B ae identical, they cannot be"
	       , a->name
	       , b->name
	       ) ;
 
  if (a == c)
      messcrash ("mxContractFirstOrLastIndex A = <B | C>, A and C ae identical, they cannot be"
	       , a->name
	       , c->name
	       ) ;
 
  if (a->type < b->type)
    messcrash ("mxContractFirstOrLastIndex %s %s %s, type mismatch, a->type must be at least as precise as b->type"
	       , a->name
	       , b->name
	       , c->name
	       ) ;
  if (a->type < c->type)
    messcrash ("mxContractFirstOrLastIndex %s %s %s, type mismatch, a->type must be at least as precise as c->type"
	       , a->name
	       , b->name
	       , c->name
	       ) ;

  
  if (nFirst)
    {
      for (ii = 0, gMax = 1 ; ii < nFirst ; ii++)
	gMax *= b->shapes[ii] ;
      dgb = 1 ; /* g is the fastest index of tensor b */
      dgc = 1 ; /* g is the fastest index of tensor c */
      kbMax = b->size / gMax ;
      kcMax = c->size / gMax ;
    }
  else if (nLast)
    {
      for (ii = 0, gMax = 1 ; ii < nLast ; ii++)
	gMax *= b->shapes[b->rank - ii - 1] ;
      dgb =  b->deltas[b->rank - nLast] ; /* g is the slowest index of tensor b */
      dgc =  c->deltas[c->rank - nLast] ; /* g is the slowest index of tensor c*/
      kbMax = b->size / gMax ;
      kcMax = c->size / gMax ;
    }

  if (a->size * gMax * gMax != b->size * c->size)
    messcrash ("mxContractFirstOrLastIndex %s %s %s, size mismatch, a->shapes should match (b->shapes  union c-<shapes, excluding the zeroth index of b and c"
	       , a->name
	       , b->name
	       , c->name
	       ) ;
  
  for (kb = 0 ; kb < kbMax ; kb++)
    {
      int kb1 = nFirst ? kb * gMax : kb ; /* required offeset in tensor b */
      for (kc = 0 ; kc < kcMax ; kc++)
	{
	  int kc1 = nFirst ? kc * gMax : kc ; /* required offeset in tensor c */
	  int ka1 = kb + kc * kbMax ; /* required offeset in tensor a */
	  switch (a->type)
	    {
	    case MX_NULL:
	    case MX_BOOL:
	    case MX_INT:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zi = 0 ;
		      zia = a->zi + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [ g * dgb] * zic [g * dgc] ;
		      *zia = zi ;
		      break ;
		    case MX_FLOAT:
		      break ;
		    case MX_COMPLEX:
		      break ;
		    }
		  break ;
		default:
		  break ;
		}
	      break ;

	    case MX_FLOAT:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zi = 0 ;
		      zfa = a->zf + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [ g * dgb] * zic [g * dgc] ;
		      *zfa = zi ;
		      break ;
		    case MX_FLOAT:
		      zf = 0 ;
		      zfa = a->zf + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zfc = c->zf + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zib [g] * zfc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zib [ g * dgb] * zfc [g * dgc] ;
		      *zfa = zf ;
		      break ;
		    default:
		      break ;
		    }
		  
		  break ;
		case MX_FLOAT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zf = 0 ;
		      zfa = a->zf + ka1 ;
		      zfb = b->zf + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [ g * dgb] * zic [g * dgc] ;
		      *zfa = zf ;
		      break ;
		    case MX_FLOAT:
		      zf = 0 ;
		      zfa = a->zf + ka1 ;
		      zfb = b->zf + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zfc = c->zf + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [g] * zfc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [ g * dgb] * zfc [g * dgc] ;
		      *zfa = zf ;
		      break ;
		    default:
		      break ;
		    }
		  break ;
		default:
		  break ;
		}
	      break ;
	    case MX_COMPLEX:
	      switch (b->type)
		{
		case MX_NULL:
		case MX_BOOL:
		case MX_INT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zi = 0 ;
		      zca = a->zc + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zi = zi + zib [ g * dgb] * zic [g * dgc] ;
		      *zca = zi ;
		      break ;
		    case MX_FLOAT:
		      zf = 0 ;
		      zca = a->zc + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zfc = c->zf + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zib [g] * zfc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zib [ g * dgb] * zfc [g * dgc] ;
		      *zca = zf ;
		      break ;
		    case MX_COMPLEX:
		      zc = 0 ;
		      zca = a->zc + ka1 ;
		      zib = b->zi + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zcc = c->zc + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zib [g] * zcc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zib [ g * dgb] * zcc [g * dgc] ;
		      *zca = zc ;
		      break ;
		    }
		  
		  break ;
		case MX_FLOAT:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zf = 0 ;
		      zca = a->zc + kb1 + kc * kbMax ;
		      zfb = b->zf + kb ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [ g * dgb] * zic [g * dgc] ;
		      *zca = zf ;
		      break ;
		    case MX_FLOAT:
		      zf = 0 ;
		      zca = a->zc + ka1 ;
		      zfb = b->zf + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zfc = c->zf + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [g] * zfc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zf = zf + zfb [ g * dgb] * zfc [g * dgc] ;
		      *zca = zf ;
		      break ;
		    case MX_COMPLEX:
		      zc = 0 ;
		      zca = a->zc + ka1 ;
		      zfb = b->zf + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zcc = c->zc + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zfb [g] * zcc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zfb [ g * dgb] * zcc [g * dgc] ;
		      *zca = zc ;
		      break ;
		    }
		  break ;
		case MX_COMPLEX:
		  switch (c->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      zc = 0 ;
		      zca = a->zc + ka1 ;
		      zcb = b->zc + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zic = c->zi + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [g] * zic [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [ g * dgb] * zic [g * dgc] ;
		      *zca = zc ;
		      break ;
		    case MX_FLOAT:
		      zc = 0 ;
		      zca = a->zc + ka1 ;
		      zcb = b->zc + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zfc = c->zf + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [g] * zfc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [ g * dgb] * zfc [g * dgc] ;
		      *zca = zc ;
		      break ;
		    case MX_COMPLEX:
		      zc = 0 ;
		      zca = a->zc + ka1 ;
		      zcb = b->zc + kb1 ;  /* we are contracting the zeroth fast index ob */
		      zcc = c->zc + kc1 ;  /* we are contracting the zeroth fast index ob */
		      if (dgb * dgc == 1)
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [g] * zcc [g] ;
		      else
			for (g = 0 ; g < gMax ; g++)
			  zc = zc + zcb [ g * dgb] * zcc [g * dgc] ;
		      *zca = zc ;
		      break ;
		    }
		  break ;
		}
	      break ;
	    }
	}
    }
  
  return a ;
} /* mxContractFirstOrLastIndex */

/**********************************************************************/

MX mxContractFirstIndex (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstOrLastIndex (a, b, c, 1, 0, h) ;
} /* mxContractFirstIndex */

/**********************************************************************/

MX mxContractLastIndex (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstOrLastIndex (a, b, c, 0, 1, h) ;
} /* mxContractLastIndex */

/**********************************************************************/

MX mxContractFirstNIndices (int n, MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstOrLastIndex (a, b, c, n, 0, h) ;
} /* mxContractFirstTwoIndices */

/**********************************************************************/
/* mxDot: corresponds to the fastest tensor product
 *   we trace the zeroth (fastest moving) index of both tensors
 * for tensors A of rank 2, M and N M.N would be the usual matrix product of Mt N
 
 * mxDot is just a convenient short hand for mxContractFirstIndex
*/

MX mxDot (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstIndex (a, b, c, h) ;
} /* mxDot */

/**********************************************************************/
/* mxMatMult corresponds to the standard matrix product
 *   we trace index 1 of b with index 0 of c
*/

MX mxMatMult (MX b, MX c, AC_HANDLE h)
{
  if (!b)
    messcrash ("matrix C is NULL in maMatMult") ;
  if (!c)
    messcrash ("matrix B is NULL in maMatMult") ;
  if (b->rank != 2)
    messcrash ("matrix B not of rank 2 in maMatMult") ;
  if (c->rank != 2)
    messcrash ("matrix C not of rank 2 in maMatMult") ;
  if (b->shapes[1] != c->shapes[0])
    messcrash ("mxMatMult: matrix %s shape(%d,%d) cannot multiply matrix %s shape(%d,%d)\n"
		   , b->name
		   , b->shapes[0]
		   , b->shapes[1]
		   , c->name
		   , c->shapes[0]
		   , c->shapes[1]
		   ) ;



  return mxContract (0, b, c, 1, 0, h) ;
} /* mxMatMult */

/**********************************************************************/

MX mxMatMultiProduct (AC_HANDLE h, MX *mms)
{
  int ii ;
  MX m = mms[0] ;

  if (! m)
    return 0 ;
  for (ii = 1 ; mms[ii] ; ii++)
    m = mxMatMult(m, mms[ii], h) ;

  return m ;
} /* mxMatMultiProduct */

/**********************************************************************/
/* mxMatListMult
 *   multiply a zero terminated list of rank-4 matrices
 * example:   m = mxMaltListMult (h, m1,m2,m3...m77, 0) 
 *    returns m = m1*m2*m3*...*m77, allocated on handle h
*/

MX mxMatListMult (AC_HANDLE h0, MX *ap)
{
  MX m ;
  int i, n ;
  AC_HANDLE h = ac_new_handle () ;

  if (!ap)
    messcrash ("Null list received by mxMatListMult") ;
  for (n = 0 ; n < 100 ; n++)
    {
      m = ap[n] ;
      if (!m)
	break ;
      if (m->rank != 2)
	messcrash ("matrix number %d: %s in list is of rank %d != 2 in mxMatListMult", n, m->name, m->rank) ;
      if (n && m->shapes[0] != ap[n-1]->shapes[1])
	messcrash ("mxMatListMult list of matrices: matrix %d:%s shape(%d,%d) cannot multiply matrix %d:%s shape(%d,%d)\n"
		   , n -1		   
		   , ap[n-1]->name
		   , ap[n-1]->shapes[0]
		   , ap[n-1]->shapes[1]
		   , n
		   , ap[n]->name
		   , ap[n]->shapes[0]
		   , ap[n]->shapes[1]
		   ) ;
    }
  if (n >= 100)
    messcrash ("More than 100 matrices in mxMatListMult, probably the list is not zero terminated") ;

  if (n == 0)
    messcrash ("Zero matrix in mxMatListMult list") ;
  else if (n == 1)
    {
      m = mxCopy (0, ap[0], h0) ;
    }
  else 
    {
      m = 0 ;
      for (i = 1 ; i < n ; i++)
	{
	  MX m1 = i == 1 ? ap[0] : m ;
	  MX m2 = ap[i] ;
	  m = mxMatMult (m1, m2, i == n-1 ? h0 : h) ;
	}
    }
  ac_free (h) ;
  return m ;
} /* mxMatListMult */

/**********************************************************************/
/* full trace of a single square matrix */
float complex  mxMatTrace (MX a)
{
  int k, kMax = 0, dk = 0 ;
  int zi = 0, * Restrict zia ;
  float zf = 0, * Restrict zfa ;
  float complex z = 0, zc = 0, * Restrict zca ;


  if (!a)
    messcrash ("NULL matrix in mxMatTrace") ;
  if (a->rank != 2)
    messcrash ("matrix %s (rank %d) not of rank 2 in mxMatTrace", a->name, a->rank) ;
  if (a->shapes[0] != a->shapes[1])
    messcrash ("Non square matrix %s shape=(%d %d)in mxMatTrace"
	       ,  a->name, a->shapes[0], a->shapes[1]
	       ) ;

  kMax = a->shapes[0] ;
  dk = kMax + 1 ;
  kMax *= kMax ;
  switch (a->type)
    {
    case MX_NULL:
      break ;
    case MX_BOOL:
    case MX_INT:
      zia = a->zi ;
      for (k = 0 ; k < kMax ; k += dk)
	zi = zi + zia [k] ;
      z = zi ;
      break ;
    case MX_FLOAT:
      zfa = a->zf ;
      for (k = 0 ; k < kMax ; k += dk)
	zf = zf + zfa [k] ;	
      z = zf ;
      break ;
    case MX_COMPLEX:
      zca = a->zc ;
      for (k = 0 ; k < kMax ; k += dk)
	zc = zc + zca [k] ;
      z = zc ;
      break ;
    }

  return z ;
} /* mxMatTrace */

/**********************************************************************/
/**********************************************************************/
/* numerical cofactor of a buffer of integers */
static BOOL mxIntCofactor (int *cf, int *aa, int ii, int jj, int rank)
{
  int i ;
  int r = rank -1 ;
  
  if (r <= 0)
    messcrash ("mxIntCofactor received rank=%d < 2", rank) ;

  for (i = 0 ; i < r ; i++)
    {    /* construct the cofactor */
      int j, di = (i >= ii ? 1 : 0) ;
      for (j = 0 ; j < r ; j++)
	{
	  int dj = (j >= jj ? 1 : 0) ;
	  cf [r * i + j] = aa [rank * (i+di) + dj + j] ;
	}
    }
  return TRUE ;
} /* mxIntCofactor */

/**********************************************************************/
/* numerical determinant or a buffer of integers */
int mxIntDeterminant (int *aa, int rank)
{
  int xx = 0, ii, sign = -1 ; 
  if (rank == 1)
    return aa[0] ;
  for (ii = 0 ; ii < rank ; ii++)
    {
      int r = rank -1 ;
      int r2 = r * r ;
      int bb[r2] ;
      sign = -sign ; /* alternate the signs */
      mxIntCofactor (bb, aa, ii, 0, rank) ;
      xx += sign * aa[rank * ii + 0] * mxIntDeterminant (bb, r) ;
    }
  return xx ;
} /* mxIntDeterminant */

/**********************************************************************/
/* returns the determinant and the NOT SCALED numerical inverse or a buffer of integers */
int mxIntInverse (int *ai, int *aa, int rank)
{
  int ii, jj ;
  int r = rank -1 ;
  int r2 = r * r ;
  int det = mxIntDeterminant (aa, rank) ;
  int s = det >= 0 ? 1 : -1 ;

  if (rank == 1)
    {
      ai[0] = s ;
    }
  else
    {
      for (ii = 0 ; ii < rank ; ii++)
	for (jj = 0 ; jj < rank ; jj++)
	  {
	    /* inverse = transposed determinant of cofactors with alterbated signs */
	    int bb[r2] ;
	    mxIntCofactor (bb, aa, ii, jj, rank) ;
	    ai [rank*jj + ii] = s * (1 - 2 * ((ii+jj) % 2)) * mxIntDeterminant (bb, r) ;
	  }
    }
  return s * det ;
} /* mxIntInverse */

/**********************************************************************/

BOOL mxMatDeterminant (MX a, int *ip, float *fp, complex float *cp) 
{
  BOOL ok = TRUE ;

  if (!a)
    messcrash ("NULL matrix in mxMatDeterminant") ;
  if (a->rank != 2)
    messcrash ("matrix %s (rank %d) not of rank 2 in mxMatTrace", a->name, a->rank) ;
  if (a->shapes[0] != a->shapes[1])
    messcrash ("Non square matrix %s shape=(%d %d)in mxMatTrace"
	       ,  a->name, a->shapes[0], a->shapes[1]
	       ) ;

  switch (a->type)
    {
    case MX_NULL:
      ok = FALSE ;
      break ;
    case MX_BOOL:
    case MX_INT:
      if (! *ip)
	messcrash ("mxMatDeterminant of Int matrix, ip not provided") ;
      *ip = mxIntDeterminant (a->zi, a->shapes[0]) ;
      break ;
    default:
      ok = FALSE ;
      break ;
    }

  return ok ;
} /* mxMatDeterminant */

/**********************************************************************/
/* return NULL or the inverse matrix if available, if INT, the matrix is NOT scaled */
MX mxMatInverse (MX a, AC_HANDLE h) 
{
  MX ai = 0 ;
  if (!a)
    messcrash ("NULL matrix in mxMatTrace") ;
  if (a->rank != 2)
    messcrash ("matrix %s (rank %d) not of rank 2 in mxMatTrace", a->name, a->rank) ;
  if (a->shapes[0] != a->shapes[1])
    messcrash ("Non square matrix %s shape=(%d %d)in mxMatTrace"
	       ,  a->name, a->shapes[0], a->shapes[1]
	       ) ;
  ai = mxCopy (0, a, h) ;
  switch (a->type)
    {
    case MX_NULL:
      break ;
    case MX_BOOL:
    case MX_INT:
      {
	int r = a->shapes[0] ;
	mxIntInverse (ai->zi, a->zi, r) ;
      }
      break ;
    default:
      messcrash ("Inverse of non Int matrix not yet programmed ") ;
      break ;
    }
  return ai ;
} /* mxMatInverse */

/**********************************************************************/
/* no verif A = B*C  where A,B,C are rank*rank size square int buffers */
void  mxIntMult (int *a, int *b, int *c, int rank)
{
  int i, j, k ;
  if (! a || ! b || !c || a==b || a==c)
    messcrash ("bad call to low level mxIntMult") ;
  for (i = 0 ; i < rank ; i++)
    for (j = 0 ; j < rank ; j++)
      {
	int x = 0 ;
	for (k = 0 ; k < rank ; k++)
	  x += b[rank * i + k] * c[rank*k + j] ;
	a[rank*i+j] = x ;
      }
}

/**********************************************************************/
/**********************************************************************/
/**********************************************************************/
/* mxD2Dot: corresponds to a frequent 2D matrix product
 *   we trace the zeroth and first (fastest moving) indices of both tensors
 * In Einstein's notation
      A_kl..mn.. = B_ijkl..  C_ijmn...
 * This is useful for 2D convolution
 * It is equivalent to a trace for a single matrix
*/

MX mxD2Dot (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstNIndices (2, a, b, c, h) ;
} /* mxD2Dot */

/**********************************************************************/
/* mxD2Dot: corresponds to a frequent 2D matrix product
 *   we trace the zeroth and first (fastest moving) indices of both tensors
 * In Einstein's notation
      A_kl..mn.. = B_ijkl..  C_ijmn...
 * This is useful for 3D convolution
*/

MX mxD3Dot (MX a, MX b, MX c, AC_HANDLE h)
{
  return mxContractFirstNIndices (3, a, b, c, h) ;
} /* mxD3Dot */

/**********************************************************************/
/**********************************************************************/

MX mxSumLastIndex (MX a, MX b, AC_HANDLE h)
{  
  AC_HANDLE h1 = ac_new_handle () ;
  int i, k = b->shapes[b->rank - 1] ;  /* size of last index */
  float z = 1.0, *zfp ;

  MX c = mxCreate (h1, "unit", MX_FLOAT, k, 0) ;

  for (i = 0, zfp = c->zf ; i < k ; i++, zfp++)
    *zfp = z ;

  a =  mxContractLastIndex (a, b, c, h) ;
  ac_free (h1) ;
  return a ;
} /* msAverageLastIndex */

/**********************************************************************/

MX mxSumFirstIndex (MX a, MX b, AC_HANDLE h)
{
  int i, k = b->shapes[0] ;  /* size of last index */
  float z = 1.0, *zfp ;

  MX c = mxCreate (0, "unit", MX_FLOAT, k, 0) ;

  for (i = 0, zfp = c->zf ; i < k ; i++, zfp++)
    *zfp = z ;

  a =  mxContractFirstIndex (a, b, c, h) ;
  ac_free (c) ;
  return a ;
} /* msAverageFirstIndex */

/**********************************************************************/
/* full trace of 2 matrices with identical shapes */
float complex  mxFullContraction (MX a, MX b, int *nGoodp, int *nBadp)
{
  int i, g, gMax ; 
  float complex zc = 0 ;
  double zf = 0 ;
  int * Restrict zia, * Restrict zib ;
  float * Restrict zfa, * Restrict zfb ;
  complex float * Restrict zca, * Restrict zcb ;

  mxCheck (a, " mxFullContraction (a, ...)") ;
  mxCheck (b, " mxFullContraction (..., b)") ;

  if (a->rank != b->rank)
    messcrash ("mxFullContraction The 2 tensors should have equal rank\n  (%s)->rank = %d != (%s)->rank = %d"
	       , a->name
	       , a->rank
	       , b->name
	       , b->rank
	       ) ;
  for (i = 0 ; i < b->rank ; i++)
    if (a->shapes[i] != b->shapes[i])
      messcrash ("mxFullContraction (%s)->shapes[%d] = %d != (%s)->shapes[%d] = %d"
		 , a->name
		 , i
		 , a->shapes[i]
		 , b->name
		 , i
		 , b->shapes[i]
		 ) ;
      
  gMax = a->size ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      switch (b->type)
	{
	case MX_NULL:
	case MX_BOOL:
	case MX_INT:
	  zia = a->zi ;
	  zib = b->zi ;	 
	  for (g = 0 ; g < gMax ; g++)
	    zf = zf + zia[g] * zib[g] ;
	  break ;
	case MX_FLOAT:
	  zia = a->zi ;
	  zfb = b->zf ;
	  for (g = 0 ; g < gMax ; g++)
	    zf = zf + zia[g] * zfb[g] ;
	  break ;
	case MX_COMPLEX:
	  zia = a->zi ;
	  zcb = b->zc ;
	  for (g = 0 ; g < gMax ; g++)
	    zc = zc + zia[g] * zcb[g] ; 
	  break ;
	}
      break ;
    case MX_FLOAT:
      switch (b->type)
	{
	case MX_NULL:
	case MX_BOOL:
	case MX_INT:
	  zfa = a->zf ;
	  zib = b->zi ;
	  for (g = 0 ; g < gMax ; g++)
	    zf = zf + zfa[g] * zib[g] ;
	  break ;
	case MX_FLOAT:
	  zfa = a->zf ;
	  zfb = b->zf ;
	  for (g = 0 ; g < gMax ; g++)
	    zf = zf + zfa[g] * zfb[g] ;
	  if (nGoodp)
	    {
	      int ng = 0, nb = 0 ;
	      int i, iMax = a->shapes[a->rank - 1] ;
	      int j, jMax = a->size/iMax ;

	      zfa = a->zf ;
	      zfb = b->zf ;
	      for (i = 0 ; i < iMax ; i++)
		{
		  float aBest = zfa[i * jMax] ;
		  float bBest = zfb[i * jMax] ;
		  int ja = 0, jb = 0 ;
		  for (j = 1 ; j < jMax ; j++)
		    {
		      if (zfa[j + i * jMax] > aBest)
			{ aBest = zfa[j + i * jMax] ; ja = j ; }
		      if (zfb[j + i * jMax] > bBest)
			{ bBest = zfb[j + i * jMax] ; jb = j ; }
		    }
		  if (ja == jb)
		    ng++ ;
		  else
		    nb++ ;
		}
	      *nGoodp = ng ;
	      *nBadp = nb ;
	    }
	  break ;
	case MX_COMPLEX:
	  zfa = a->zf ;
	  zcb = b->zc ;
	  for (g = 0 ; g < gMax ; g++)
	    zc = zc + zfa[g] * zcb[g] ;
	  break ;
	}
      break ;
    case MX_COMPLEX:
      switch (b->type)
	{
	case MX_NULL:
	case MX_BOOL:
	case MX_INT:
	  zca = a->zc ;
	  zib = b->zi ;
	  for (g = 0 ; g < gMax ; g++)
	    zc = zc + zca[g] * zib[g] ;
	  break ;
	case MX_FLOAT:
	  zca = a->zc ;
	  zfb = b->zf ;
	  for (g = 0 ; g < gMax ; g++)
	    zc = zc + zca[g] * zfb[g] ;	
	  break ;
	case MX_COMPLEX:
	  zca = a->zc ;
	  zcb = b->zc ;
	  for (g = 0 ; g < gMax ; g++)
	    zc = zc + zca[g] * zcb[g] ;
	  break ;
	}
      break ;
    }

  if (zf) zc = zf ;
  return zc ;
  } /* mxFullContraction */

/**********************************************************************/
/* Froebenius norm, sum of squares of all cells */
double mxFNorm (MX a)
{
  int g, gMax ; 
  double z = 0 ;
 
  mxCheck (a, " mxFNorm (a)") ;
  gMax = a->size ;
  switch (a->type)
    {
    case MX_NULL:
      break ;
    case MX_BOOL:
    case MX_INT:
      {
	int * Restrict zia = a->zi ;
	for (g = 0 ; g < gMax ; g++)
	  z = z + zia[g] * zia[g] ;
      }
      break ;
    case MX_FLOAT:
      {
	float * Restrict zfa = a->zf ;
	for (g = 0 ; g < gMax ; g++)
	  z = z + zfa[g] * zfa[g] ;
      }
      break ;
    case MX_COMPLEX:
      {
	complex float * Restrict zca = a->zc ;
	for (g = 0 ; g < gMax ; g++)
	  z = z + creal (zca[g] * conj (zca[g])) ;
      }
      break ;
    }
  
  return z ;
  } /* mxFNorm */

/**********************************************************************/
/**********************************************************************/

void mxShow (MX a)
{
  int ii, i0, i1, i2, i3 ;
  int i3Max = (a->shapes[3] ? a->shapes[3] : 1) ;
  int i2Max = (a->shapes[2] ? a->shapes[2] : 1) ;
  int i1Max = (a->shapes[1] ? a->shapes[1] : 1) ;
  int i0Max = (a->shapes[0] ? a->shapes[0] : 1) ;

  if (!a || a->magic != &MX_MAGIC)
    {
      fprintf (stderr, "bad matrix\n") ;
      return ;
    }
  
  for (i3 = 0 ; i3 < i3Max ; i3++)
    {
      fprintf (stderr, "\n") ;
      for (i2 = 0 ; i2 < i2Max ; i2++)
	{
	  fprintf (stderr, "\n") ;
	  for (i0 = 0 ; i0 < i0Max ; i0++)
	    {
	      fprintf (stderr, "\n %d:%d:.:%d", i3, i2, i0) ;
	      ii = i3 * a->deltas[3] + i2 * a->deltas[2] + i0 ;
	      for (i1 = 0 ; i1 < i1Max ; i1++)
		switch (a->type)
		  {
		  case MX_NULL:
		  case MX_BOOL:
		  case MX_INT:
		    fprintf (stderr, "\t%d", a->zi[ii + i1* a->deltas[1]]) ;
		    break ;
		  case MX_FLOAT:
		    fprintf (stderr, "\t%f", a->zf[ii + i1* a->deltas[1]]) ;
		    break ;
		  case MX_COMPLEX:
		    fprintf (stderr, "\t%f + %f*i", creal(a->zc[ii+  i1* a->deltas[1]]), cimag(a->zc[ii+  i1* a->deltas[1]])) ;
		    break ;
		  }
	    }
	}
    } 
  fprintf (stderr, "\n") ;
  fprintf (stderr, "%s\n", a->name) ;
  return ;
} /* mxShow */

/**********************************************************************/

void mxSet (MX a, const void *x)
{  
  mxCheck (a, "mxSet") ;

  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
      break ;
    case MX_INT:
      memcpy (a->zi, x, a->size * a->typeSize) ;
      break ;
    case MX_FLOAT:
      memcpy (a->zf, x, a->size * a->typeSize) ;
      break ;
    case MX_COMPLEX:
      memcpy (a->zc, x, a->size * a->typeSize) ;
      break ;
    }
  return ;
} /* mxSet */

/**********************************************************************/
/* reset to zero */
void mxClear (MX a)
{
  mxCheck (a, "mxClear") ;
   
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
      break ;
    case MX_INT:
      memset (a->zi, 0, a->typeSize * a->size) ;
      break ;
    case MX_FLOAT:
      memset (a->zf, 0, a->typeSize * a->size) ;
      break ;
    case MX_COMPLEX:
      memset (a->zc, 0, a->typeSize * a->size) ;
      break ;
    }
  return ;
} /* mxClear */

/**********************************************************************/

void mxRandomInit (MX a)
{
  int i ;
  float complex *cp ;
  float *fp ;
  int *ip ;

  mxCheck (a, "mxRandomInit") ;

  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
      break ;
    case MX_INT:
      for (i = 0 , ip = a->zi ; i < a->size ; i++)
	ip[i] = randint () % 100 ;
      break ;
    case MX_FLOAT:
      for (i = 0 , fp = a->zf ; i < a->size ; i++)
	fp[i] = randfloat () / 1000 ;
      break ;
    case MX_COMPLEX:
      for (i = 0 , cp = a->zc ; i < a->size ; i++)
	cp[i] = randfloat() + randfloat () * _Complex_I / 1000 ;
      break ;
    }
  return ;
}

/**********************************************************************/

MX mxRetype (MX a, MX b, MX_TYPE type, AC_HANDLE h)
{
  mxCheck (b, "mxRetype (...,b,...)") ;
  if (!a)					
    {
      MM mm ;
      mm.type = type ;
      a = mxDuplicate (hprintf (h, "retype %s", b->name), a, b, &mm, h) ;
    }
  else
    {
      if (type == a->type)
	{
	  switch (a->type)
	    {
	    case MX_NULL:
	    case MX_BOOL:
	    case MX_INT:
	      memset (a->zi, 0, a->size * a->typeSize) ;
	      break ;
	    case MX_FLOAT:
	      memset (a->zf, 0, a->size * a->typeSize) ;
	      break ;
	    case MX_COMPLEX:
	      memset (a->zc, 0, a->size * a->typeSize) ;
	      break ;
	    }
	}
      else
	{
	  switch (a->type)
	    {
	    case MX_NULL:
	    case MX_BOOL:
	    case MX_INT:
	      ac_free (a->zi) ;
	      break ;
	    case MX_FLOAT:
	      ac_free (a->zf) ;
	      break ;
	    case MX_COMPLEX:
	      ac_free (a->zc) ;
	      break ;
	    }
	  switch (type)
	    {
	    case MX_NULL:
	    case MX_BOOL:
	    case MX_INT:
	      a->typeSize = sizeof (int) ;
	      a->zi = halloc (a->typeSize * a->size, 0) ;
	      break ;
	    case MX_FLOAT:
	      a->typeSize = sizeof (float) ;
	      a->zf = halloc (a->typeSize * a->size, 0) ;
	      break ;
	    case MX_COMPLEX:
	      a->typeSize = sizeof (complex float) ;
	      a->zc = halloc (a->typeSize * a->size, 0) ;
	      break ;
	    }
	}
      a->type = type ;
    }
  mxAdd (a, a, b, h) ;
  
  return a ;
} /* mxRetype */

/**********************************************************************/

MX mxElementWiseMultiplication (MX a, MX b, MX c, AC_HANDLE h)
{
  int i, iMax ;

  mxCheck (b, "mxElementWiseMultiplicationmxRetype (...,b,...)") ;
  mxCheck (c, "mxElementWiseMultiplicationmxRetype (...,c,...)") ;
  if (!a)
    a = mxDuplicate ("mxElementWiseMultiplication", a , b, c, h) ;

  if (0 && a == b)
      messcrash ("mxElementWiseMultiplication A = B * C, A and B are identical, they cannot be"
	       , a->name
	       , b->name
	       ) ;

  mxCheckSameType (a, b, "mxElementWiseMultiplication") ;
  mxCheckSameType (a, c, "mxElementWiseMultiplication") ;
  mxCheckSameRank (a, b, "mxElementWiseMultiplication") ;
  mxCheckSameRank (a, c, "mxElementWiseMultiplication") ;
  mxCheckSameShapes (a, b, "mxElementWiseMultiplication") ;
  mxCheckSameShapes (a, c, "mxElementWiseMultiplication") ;

  iMax = a->size ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      {
	int * Restrict zia = a->zi ;
	int * Restrict zib = b->zi ;
	int * Restrict zic = c->zi ;

	for (i = 0 ; i < iMax ; i++)
	  zia[i] = zib[i] * zic[i] ;
      }
      break ;
    case MX_FLOAT:
      {
	float * Restrict zfa = a->zf ;
	float * Restrict zfb = b->zf ;
	float * Restrict zfc = c->zf ;

	for (i = 0 ; i < iMax ; i++)
	  zfa[i] = zfb[i] * zfc[i] ;
      }

      break ;
    case MX_COMPLEX:
      {
	complex float * Restrict zca = a->zc ;
	complex float * Restrict zcb = b->zc ;
	complex float * Restrict zcc = c->zc ;

	for (i = 0 ; i < iMax ; i++)
	  zca[i] = zcb[i] * zcc[i] ;
      }
      
      break ;
    }
  
  return a ;
} /* mxElementWiseMultiplication */

/**********************************************************************/
/* Matrix conjugate */
MX  mxConjugate (MX ac, MX a, AC_HANDLE h)
{ 
  mxCheck (a, "mxConjugate (.., a,...)") ;
  if (! ac)
    ac = mxDuplicate ("conjugated", ac, a, 0, h) ;
  mxCheck (a, "mxConjugate (.., a,...)") ;
  mxCheckSameRank (a, ac, "mxConjugate") ;
  mxCheckSameType (a, ac, "mxConjugate") ;
  mxCheckSameShapes (a, ac, "mxConjugate") ;

  if (a->type == MX_COMPLEX)
    {
      complex float *Restrict zc = ac->zc ;
      int i ;
      for (i = 0 ; i < a->size ; i++)
	zc[i] = conjf (zc[i]) ;
    }
  return ac ;
} /* mxConjugate */

/**********************************************************************/
/* Matrix conjugate */
MX  mxMatTranspose (MX at, MX a, AC_HANDLE h)
{ 
  mxCheck (a, "mxMatTranspose (.., a,...)") ;
  if (! at)
    at = mxDuplicate ("transposed", at, a, 0, h) ;
  mxCheck (a, "mxTranspose (.., a,...)") ;
  mxCheckSameRank (a, at, "mxTranspose") ;
  mxCheckSameType (a, at, "mxTranspose") ;
  mxCheckSameShapes (a, at, "mxTranspose") ;
  mxCheckRankN (a, 2, "mxTranspose") ;
  at =  mxTranspose (at, a, 0, 1, h) ;

  return at ; 
} /* mxMatTranspose */

/**********************************************************************/
/* Matrix tranposition */
MX  mxTranspose (MX at, MX a, int g1, int g2, AC_HANDLE h)
{ 
  int *da ;
  int *db = a->deltas ;
  int da1, da2, da3, delta ;
  int db0, db1, db2, db3 ; 
  int kk[4], ssa[4] ;
  int ii, jj ;
  int i0, i1, i2, i3 ;
  MX_TYPE type = a->type ;

  kk[0] = 0 ; kk[1] = 1 ; kk[2] = 2 ; kk[3] = 3 ;
  kk[g1] = g2 ; kk[g2] = g1 ;
  mxCheck (a, "mxTranspose (.., a,...)") ;
  if (! at)
    at = mxDuplicate ("transposed", at, a, 0, h) ;
  at->shapes[g1] = a->shapes[g1] ;
  at->shapes[g2] = a->shapes[g2] ;
  mxCheck (a, "mxTranspose (.., a,...)") ;
  mxCheckSameRank (a, at, "mxTranspose") ;
  mxCheckSameType (a, at, "mxTranspose") ;
  mxCheckSameShapes (a, at, "mxTranspose") ;
  at->shapes[g1] = a->shapes[g2] ;
  at->shapes[g2] = a->shapes[g1] ;
  if (a == at)
    messcrash ("A == B %s = %sin mxTranspose (a, b)"
	       , at->name
	       , a->name
	       ) ;

  da = at->deltas ;
  mxCheck (at, "mxTrnspose (at,...)") ;

  if (a->rank > 4)
    messcrash ("rank = %d > 4 in mxTranspose", a->rank) ;

  for (i3 = 0 ; i3 < 4 ; i3++)
    {
      ssa[i3] = at->shapes[i3] ? at->shapes[i3] : 1 ; 
    }
  for (ii = 0 ; ii < at->rank ; ii++)
    {
      if (ii == 0)
	delta = 1 ;
      else
	delta *=  at->shapes[ii -1] ;
      at->deltas[ii] = delta ;
    }

  if (g1 < 0 || g1 > 3)
    messcrash ("bad g1 = %d received by mxTranspose", g1) ;
  if (g2 < 0 || g2 > 3)
    messcrash ("bad g2 = %d received by mxTranspose", g2) ;
  
  da1 = da[1] ;
  da2 = da[2] ;
  da3 = da[3] ;

  db0 = db[kk[0]] ;
  db1 = db[kk[1]] ;
  db2 = db[kk[2]] ;
  db3 = db[kk[3]] ;

  for (i3 = 0 ; i3 < ssa[3] ; i3++) 
    {
      for (i2 = 0 ; i2 < ssa[2] ; i2++)
	{
	  for (i1 = 0 ; i1 < ssa[1] ; i1++)
	    {
	      ii =  i3 * da3 + i2 * da2 + i1 * da1 ;
	      jj =  i3 * db3 + i2 * db2 + i1 * db1 ;
	      for (i0 = 0 ; i0 < ssa[0] ; i0++, ii++, jj += db0)
		{
		   switch (type)
		     {
		     case MX_NULL:
		     case MX_BOOL:
		     case MX_INT:
		       at->zi[ii] = a->zi[jj] ; 
		       break ;
		     case MX_FLOAT:
		       at->zf[ii] = a->zf[jj] ; 
		       break ;
		     case MX_COMPLEX:
		       at->zc[ii] = a->zc[jj] ; 
		       break ;
		     }
		}
	    }
	}
    }
  return at ;
} /* mxTranspose */

/**********************************************************************/
/* Matrix contraction */
MX mxContract (MX a, MX b, MX c, int g1, int g2, AC_HANDLE h)
{
  int i, j, k, g, dgb, dgc ;
  int gMax ;  /* size of the contacted dimemsion */
  int kk[4] ; /* correspondance from the a shapes to b or c shape, if kk[] % 4 = 0 or 1  */
  int ss[4] ; /* shapes of the a matrix */
 
  int i0, i1, i2, i3 ;
  int ib[4] ;
  int ic[4] ;

  int ja0, ja1, ja2, ja3 ;
  int jb0, jb1, jb2, jb3 ;
  int jc0, jc1, jc2, jc3 ;

  int jja, jjb, jjc ;

  int *sa ;
  int *da, *db, *dc ;

  int zi, *zia, *zib, *zic ;
  int *zzia, *zzib, *zzic ;
  float zf, *zfa, *zfb, *zfc ;
  float *zzfa, *zzfb, *zzfc ;
  complex float zc, *zca, *zcb, *zcc ;
  complex float *zzca, *zzcb, *zzcc ;

  
  mxCheck (b, "mxContract (..,b,..)") ;
  mxCheck (c, "mxContract (..,c,..)") ; 

  if (a && a->rank > 4)
    messcrash ("A rank = %d > 4 in mxContract", a->rank) ;
  if (b->rank > 4)
    messcrash ("B rank = %d > 4 in mxContract", b->rank) ;
  if (c->rank > 4)
    messcrash ("C rank = %d > 4 in mxContract", c->rank) ;

  if (a == b)
    messcrash ("A and B must differ in mxContract (a,b,c..)");
  if (a == c)
    messcrash ("A and C must differ in mxContract (a,b,c..)");
  memset (kk, 0, sizeof (kk)) ;
  memset (ss, 0, sizeof (ss)) ;
  if (1)
    {
      int rank = b->rank + c->rank - 2 ;
      MX_TYPE type = b->type > c->type ? b->type : c->type ;
      
      for (i = j = 0 ; i < b->rank ; i++)
	if (i != g1)
	  { ss[j] = b->shapes[i] ;kk[j++] = i ; }
      for (i = 0 ; i < c->rank ; i++)
	if (i != g2)
	  { ss[j] = c->shapes[i] ; kk[j++] = 4 + i ; }
      
      if (! a)
	a = mxCreate (h, "mxContract", type, ss[0], ss[1], ss[2], ss[3], 0) ;
      else
	{
	  if (a->type < type)
	    messcrash ("type mismatch  in mxContract a->type = %d b->type=%d c->type = %d"
		       , a->type
		       , b->type
		       , c->type
		       ) ;

	  if (a->rank != rank)
	    messcrash ("rank mismatch in mxContract (a->rank =%d) != (b->rank = %d) + (c->rank = %d) - 2"
		       , a->rank, b->rank, c->rank) ;
	  for (i = 0 ; i < a->rank ; i++)
	    if (ss[i] != a->shapes[i])
	      messcrash ("shape mismatch in mxContract (a->shapes[%d]=%d != %c->shapes[%d]=%d"
			 , i, a->shapes[i]
			 , kk[i] > 4 ? 'c' : 'b', kk[i] % 4, ss[i]
			 ) ;

	  if (g1 < 0)
	    messcrash ("contractile shape error g1 = %d < 0", g1) ;
	  if (g2 < 0)
	    messcrash ("contractile shape error g2 = %d < 0", g2) ;
	  if (g1 >= b->rank)
	    messcrash ("contractile shape error g1 = %d >= b->rank = %d", g1, b->rank) ;
	  if (g2 >= c->rank)
	    messcrash ("contractile shape error g1 = %d >= b->rank = %d", g2, c->rank) ;
	  
	  if (b->shapes[g1] != c->shapes[g2])
	    messcrash ("contractile shape mismatch in mxContract (b->shapes[%d]=%d) != (c->shapes[%d]=%d)"
		       , g1, b->shapes[g1], g2, c->shapes[g2]
		       ) ;
	  gMax = b->shapes[g1] ;
	}
    }
  gMax = b->shapes[g1] ;
  
 
  mxCheck (a, "mxContract(a,...)") ;
  sa = a->shapes ;
  
  da = a->deltas ;
  db = b->deltas ;
  dc = c->deltas ;
  
  dgb = db[g1] ;
  dgc = dc[g2] ;

  zia = a->zi ;
  zib = b->zi ;
  zic = c->zi ;
  
  zfa = a->zf ;
  zfb = b->zf ;
  zfc = c->zf ;
  
  zca = a->zc ;
  zcb = b->zc ;
  zcc = c->zc ;
  
  ja0 = ja1 = ja2 = ja3 = 0 ;
  jb0 = jb1 = jb2 = jb3 = 0 ;
  jc0 = jc1 = jc2 = jc3 = 0 ;

  for (i3 = 0 ; i3 < (sa[3] ?  sa[3] : 1) ; i3++)
    {
      ja3 = i3 * da[3] ;
      k = kk[3] ;
      if (k & 4) 
	{ k -= 4 ; ic[k] = i3 ; jc3 = ic[k] * dc[k] ; }
      else 
	{ ib[k] = i3 ; jb3 = ib[k] * db[k] ; }
      
      for (i2 = 0 ; i2 < (sa[2] ? sa[2] : 1) ; i2++)
	{
	  ja2 = i2 * da[2] ;
	  k = kk[2] ;
	  if (k & 4) 
	    { k -= 4 ; ic[k] = i2 ; jc2 = ic[k] * dc[k] ; }
	  else 
	    { ib[k] = i2 ; jb2 = ib[k] * db[k] ; }
	  
	  
	  for (i1 = 0 ; i1 < (sa[1] ? sa[1] : 1); i1++)
	    {
	      ja1 = i1 * da[1] ;
	      k = kk[1] ;
	      if (k & 4) 
		{ k -= 4 ; ic[k] = i1 ; jc1 = ic[k] * dc[k] ; }
	      else 
		{ ib[k] = i1 ; jb1 = ib[k] * db[k] ; }
	      for (i0 = 0 ; i0 < sa[0] ; i0++)
		{
		  ja0 = i0 * da[0] ;
		  k = kk[0] ;
		  if (k & 4) 
		    { k -= 4 ; ic[k] = i0 ; jc0 = ic[k] * dc[k] ; }
		  else 
		    { ib[k] = i0 ; jb0 = ib[k] * db[k] ; }
		  
		  jja = ja3 + ja2 + ja1 + ja0 ; 
		  jjb = jb3 + jb2 + jb1 + jb0 ;
		  jjc = jc3 + jc2 + jc1 + jc0 ;
		  
		  switch (a->type)
		    {
		    case MX_NULL:
		    case MX_BOOL:
		    case MX_INT:
		      switch (b->type)
			{
			case MX_NULL:
			case MX_BOOL:
			case MX_INT:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzia = zia + jja ; 
			      zzib = zib + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zi = 0 ; g < gMax ; g++, zzib+=dgb, zzic+=dgc)
				zi += (*zzib) * (*zzic) ;
			      *zzia = zi ;
			      break ;
			    case MX_FLOAT:
			      break ;
			    case MX_COMPLEX:
			      break ;
			    }
			  break ;
			default:
			  break ;
			}
		      break ;
		    case MX_FLOAT:
		      switch (b->type)
			{
			case MX_NULL:
			case MX_BOOL:
			case MX_INT:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzfa = zfa + jja ; 
			      zzib = zib + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zf = 0 ; g < gMax ; g++, zzib+=dgb, zzic+=dgc)
				zf += (*zzib) * (*zzic) ;
			      *zzfa = zf ;
 			      break ;
			    case MX_FLOAT:
			      zzfa = zfa + jja ; 
			      zzib = zib + jjb ; 
			      zzfc = zfc + jjc ;
			      for (g = 0, zf = 0 ; g < gMax ; g++, zzib+=dgb, zzfc+=dgc)
				zf += (*zzib) * (*zzfc) ;
			      *zzfa = zf ;
			      break ;
			    default:
			      break ;
			    }
			  
			  break ;
			case MX_FLOAT:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzfa = zfa + jja ; 
			      zzfb = zfb + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zf = 0 ; g < gMax ; g++, zzfb+=dgb, zzic+=dgc)
				zf += (*zzfb) * (*zzic) ;
			      *zzfa = zf ;
 			      break ;
			    case MX_FLOAT:
			      zzfa = zfa + jja ; 
			      zzfb = zfb + jjb ; 
			      zzfc = zfc + jjc ;
			      for (g = 0, zf = 0 ; g < gMax ; g++, zzfb+=dgb, zzfc+=dgc)
				zf += (*zzfb) * (*zzfc) ;
			      *zzfa = zf ;
 			      break ;
			    default:
			      break ;
			    }
			  break ;
			default:
			  break ;
			}
		      break ;
		    case MX_COMPLEX:
		      switch (b->type)
			{
			case MX_NULL:
			case MX_BOOL:
			case MX_INT:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzca = zca + jja ; 
			      zzib = zib + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzib+=dgb, zzic+=dgc)
				zc += (*zzib) * (*zzic) ;
			      *zzca = zc ;
			      break ;
			    case MX_FLOAT:
			      zzca = zca + jja ; 
			      zzib = zib + jjb ; 
			      zzfc = zfc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzib+=dgb, zzfc+=dgc) 
				zc += (*zzib) * (*zzfc) ;
			      *zzca = zc ;
			      break ;
			    case MX_COMPLEX:
			      zzca = zca + jja ; 
			      zzib = zib + jjb ; 
			      zzcc = zcc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzib+=dgb, zzcc+=dgc)
				zc += (*zzib) * (*zzcc) ;
			      *zzca = zc ;
			      break ;
			    }
			  
			  break ;
			case MX_FLOAT:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzca = zca + jja ; 
			      zzfb = zfb + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzfb+=dgb, zzic+=dgc)
				zc += (*zzfb) * (*zzic) ;
			      *zzca = zc ;
			      break ;
			    case MX_FLOAT:
			      zzca = zca + jja ; 
			      zzfb = zfb + jjb ; 
			      zzfc = zfc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzfb+=dgb, zzfc+=dgc)
				zc += (*zzfb) * (*zzfc) ;
			      *zzca = zc ;
			      break ;
			    case MX_COMPLEX:
			      zzca = zca + jja ; 
			      zzfb = zfb + jjb ; 
			      zzcc = zcc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzfb+=dgb, zzcc+=dgc)
				zc += (*zzfb) * (*zzcc) ;
			      *zzca = zc ;
			      break ;
			    }
			  break ;
			case MX_COMPLEX:
			  switch (c->type)
			    {
			    case MX_NULL:
			    case MX_BOOL:
			    case MX_INT:
			      zzca = zca + jja ; 
			      zzcb = zcb + jjb ; 
			      zzic = zic + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzcb+=dgb, zzic+=dgc) 
				zc += (*zzcb) * (*zzic) ;
			      *zzca = zc ;
     			      break ;
			    case MX_FLOAT:
			      zzca = zca + jja ; 
			      zzcb = zcb + jjb ; 
			      zzfc = zfc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzcb+=dgb, zzfc+=dgc)
				zc += (*zzcb) * (*zzfc) ;
			      *zzca = zc ;
			      break ;
			    case MX_COMPLEX:
			      zzca = zca + jja ; 
			      zzcb = zcb + jjb ; 
			      zzcc = zcc + jjc ;
			      for (g = 0, zc = 0 ; g < gMax ; g++, zzcb+=dgb, zzcc+=dgc)
				zc += (*zzcb) * (*zzcc) ;
			      *zzca = zc ;
			      break ;
			    }
			  break ;
			}
		      break ;
		    }
		}
	    }
	}
    }

  return a ;
} /* mxContract */

/**********************************************************************/
/**********************************************************************/
/* Returns a = b :: W
 * a and b are matrices with same geometry
 * W, the convolution kernel may be 1, 2 or 3 dimensional
 */
MX mxConvolution (MX a, MX b, MX W, AC_HANDLE h)
{
  int i, iMax ;
  int j, jMax ;
  int k, kMax ;
  int du1, du2, u, uMax ;
  int dv1, dv2, v, vMax ;

  int dk ;
  float zf, *Restrict zfa, *Restrict zfb, *Restrict zfW ;

  messcrash ("mxConvolution pullback cache not yet defined") ;
  mxCheck (a, "mxConvolutionmxValues (a, ..)") ;
  mxCheck (b, "mxConvolution (.., b, ..)") ;
  mxCheck (W, "mxConvolution (.., W, ..)") ;
  mxCheckSameRank (a, b, "mxConvolution (a, b, ...)") ;
  mxCheckSameType (a, b, "mxConvolution (a, b, ...)") ;
  mxCheckSameType (a, W, "mxConvolution (a, .., W)") ;
  mxCheckSameShapes (a, b, "mxConvolution (a, b, ...)") ;

  if (a->rank < 2)
     messcrash (" mxConvolution (a=%s)->rank = %d < 2"
		, a->name, a->rank) ;
  if (b->rank < 2)
     messcrash (" mxConvolution (b=%s)->rank = %d < 2"
		, b->name, b->rank) ;
  if (W->rank != 2)
     messcrash (" mxConvolution (W=%s)->rank = %d != 2"
		, W->name, W->rank) ;

  if (b->type != MX_FLOAT)
    messcrash (" mxConvolution only accepts MX_FLOAT matrices") ;

  iMax = b->shapes[0] ;
  jMax = b->shapes[1] ;
  dk = b->deltas[2] ;
  kMax = b->size / dk ;
  uMax = W->shapes[0] ;
  vMax = W->shapes[1] ;
  du1 = uMax/2 ; du2 = uMax - du1 ; 
  dv1 = vMax/2 ; dv2 = vMax - dv1 ;
  
  for (k = 0 ; k < kMax ; k++)
    for (i = 0 ; i < iMax - uMax + 1 ; i++)
      for (j = 0 ; j < jMax - vMax + 1 ; j++)
	{
	  int u1, u2, v1, v2 ;
	  zfa = a->zf + i + j * iMax + k * dk ;
	  zfb = b->zf + i + j * iMax + k * dk ;
	  zfW = W->zf ;

	  /* u1, u2, v1, v2 create the equivalent of a zero padded frame around the picture */
	  zf = 0 ;
	  u1 = du1 < i ? -du1 : -i ;
	  v1 = dv1 < j ? -dv1 : -j ;
	  u2 = du2 < iMax - i ? du2 : iMax - i ;
	  v2 = dv2 < jMax - j ? dv2 : jMax - j ;

	  for (u = -u1 ; u < u2 ; u++)
	    for (v = -v1 ; v < v2 ; v++)
		zf = zf + zfb[u + v * iMax] * zfW[u + v * uMax] ;
	  zfa[0] = zf ;
	}

  return a ;
} /* mxConvolution */

/**********************************************************************/
/* Returns local max in a little square and squeezes the dimensionality
 * of the first W->rank axis (fastest moving) of b to produce a
 * c is the control matrix, same shapes as b
 * in forward mode a = ly->Z, b = c = ly->X
 * in pullback mode a = ly->dX, b = dZ, c = ly->X
 * W is used only to provide its rank and shape, the W->zf values are ignored
 */
MX mxMaxPooling (MX a, MX b, MX W, AC_HANDLE h)
{ 
  int i, iMax ;
  int j, jMax ;
  int g, gMax ;
  
  int du1, du2, u, uMax ;
  int dv1, dv2, v, vMax ;

  int dga, dgb ;
  float zf, uf, *Restrict zfa, *Restrict zfb ;
  
  mxCheck (a, "mxMaxPooling (a, ..)") ;
  mxCheck (b, "mxMaxPooling (.., b, ..)") ;
  mxCheck (W, "mxMaxPooling (.., b, ..)") ;

  mxCheckSameRank (a, b, "mxMaxPooling (a, b, ...)") ;
  mxCheckSameType (a, b, "mxMaxPooling (a, b, ...)") ;
  
  uMax = W->shapes[0] ;
  vMax = W->rank > 1 ? W->shapes[1] : 1 ;
  /* wMax = W->shapes[2] ;*/
  
  if (a->rank < W->rank)
    messcrash (" mxMaxPooling ((a=%s)->rank) = %d < ((W=%s)->rank = %d)"
	       , a->name, a->rank
	       , W->name, W->rank
	       ) ;
  

  for (i = 0 ; i < a->rank - W->rank ; i++)
    if (a->shapes [W->rank + i] != b->shapes[W->rank + i])
      messcrash ("mxMaxPooling: the spectator shapes of a and b should be equal") ;
  
  
  
  if (b->type != MX_FLOAT)
    messcrash (" mxMaxPooling only accepts MX_FLOAT matrices") ;
  if (W->type != MX_FLOAT)
    messcrash (" mxMaxPooling only accepts MX_FLOAT W") ;
  if (W->rank > 2)
    messcrash ("mxMaxPooling called with rank (W=%s) = %d > 2", W->name, W->rank) ;
  
  if (a->shapes[0] * uMax != b->shapes[0])
    messcrash (" mxMaxPooling ((a=%s)shapes[0]=%d)* (widht0=%d) != ((b=%s)shapes[0]=%d)" 
	       , a->name, a->shapes[0], uMax
	       , b->name, b->shapes[0]
	       ) ;
  if (vMax && a->shapes[1] * vMax != b->shapes[1])
    messcrash (" mxMaxPooling ((a=%s)shapes[1]=%d)* (widht1=%d) != ((b=%s)shapes[1]=%d)" 
	       , a->name, a->shapes[1], vMax
	       , b->name, b->shapes[1]
	       ) ;
  iMax = a->shapes[0] ;
  jMax = vMax ? a->shapes[1] : 1 ;

  dga = W->rank < a->rank ? a->deltas[W->rank] : a->size ;
  dgb = W->rank < b->rank ? b->deltas[W->rank] : b->size ;

  gMax = b->size / dgb ;
  du1 = 0 * uMax/2 ; du2 = uMax - du1 ; 
  dv1 = 0 * vMax/2 ; dv2 = vMax - dv1 ;

  
  for (g = 0 ; g < gMax ; g++)
    for (i = 0 ; i < iMax ; i++)
      for (j = 0 ; j < jMax ; j++)
	{
	  int uv, u1, u2, v1, v2 ;
	  zfa = a->zf + i + j * iMax + g * dga ;
	  zfb = b->zf + i * uMax + j * vMax * iMax * uMax + g * dgb ;
	  
	  zf = zfb[0] ;
	  u1 = du1 < i ? -du1 : -i ;
	  v1 = dv1 < j ? -dv1 : -j ;
	  u2 = du2 < iMax * uMax - i ? du2 : iMax * uMax - i ;
	  v2 = dv2 < jMax * vMax- j ? dv2 : jMax * vMax - j ;

	  for (u = -u1 ; u < u2 ; u++)
	    for (v = -v1 ; v < v2 ; v++)
	      {
		uv = u + v * iMax * uMax ;
		uf = zfb[uv] ;
		if (zf < uf)
		  { zf = uf ; }
	      }
	  zfa[0] = zf ;
	}
  
  return a ;
} /* mxMaxPooling */

/**********************************************************************/
/* Returns local max in a little square and squeezes the dimensionality
 * of the first W->rank axis (fastest moving) of b to produce a
 * c is the control matrix, same shapes as b
 * in forward mode a = ly->Z, b = c = ly->X
 * in pullback mode a = ly->dX, b = dZ, c = ly->X
 * W is used only to provide its rank and shape, the W->zf values are ignored
 */
MX mxMaxPoolingBack (MX a, MX b, MX c, MX W, AC_HANDLE h)
{ 
  int i, iMax ;
  int j, jMax ;
  int g, gMax ;
  
  int du1, du2, u, uMax ;
  int dv1, dv2, v, vMax ;

  int dga, dgb ;
  float zf, uf, *Restrict zfa, *Restrict zfb, *Restrict zfc ;
  
  mxCheck (a, "mxMaxPooling (a, ..)") ;
  mxCheck (b, "mxMaxPooling (.., b, ..)") ;
  mxCheck (c, "mxMaxPooling (.., c, ..)") ;
  mxCheck (W, "mxMaxPooling (.., W)") ;

  mxCheckSameRank (a, b, "mxMaxPooling (a, b, ...)") ;
  mxCheckSameType (a, b, "mxMaxPooling (a, b, ...)") ;
  mxCheckSameRank (a, c, "mxMaxPooling (a, b, ...)") ;
  mxCheckSameType (a, c, "mxMaxPooling (a, .., c, ...)") ;
  mxCheckSameShapes (a, c, "mxMaxPooling (a, .., c, ...)") ;
  
  uMax = W->shapes[0] ;
  vMax = W->rank > 1 ? W->shapes[1] : 1 ;
  /* wMax = W->shapes[2] ;*/
  
  if (a->rank < W->rank)
    messcrash (" mxMaxPooling ((a=%s)->rank) = %d < ((W=%s)->rank = %d)"
	       , a->name, a->rank
	       , W->name, W->rank
	       ) ;
  

  for (i = 0 ; i < a->rank - W->rank ; i++)
    if (a->shapes [W->rank + i] != b->shapes[W->rank + i])
      messcrash ("mxMaxPooling: the spectator shapes of a and b should be equal") ;
  
  
  
  if (b->type != MX_FLOAT)
    messcrash (" mxMaxPooling only accepts MX_FLOAT matrices") ;
  if (W->type != MX_FLOAT)
    messcrash (" mxMaxPooling only accepts MX_FLOAT W") ;
  if (W->rank > 2)
    messcrash ("mxMaxPooling called with rank (W=%s) = %d > 2", W->name, W->rank) ;
  
  if (b->shapes[0] * uMax != a->shapes[0])
    messcrash (" mxMaxPooling ((b=%s)shapes[0]=%d)* (widht0=%d) != ((a=%s)shapes[0]=%d)" 
	       , b->name, b->shapes[0], uMax
	       , a->name, a->shapes[0]
	       ) ;
  if (vMax && b->shapes[1] * vMax != a->shapes[1])
    messcrash (" mxMaxPooling ((b=%s)shapes[1]=%d)* (widht1=%d) != ((a=%s)shapes[1]=%d)" 
	       , b->name, b->shapes[1], vMax
	       , a->name, a->shapes[1]
	       ) ;
  iMax = b->shapes[0] ;
  jMax = vMax ? b->shapes[1] : 1 ;

  dga = W->rank < a->rank ? a->deltas[W->rank] : a->size ;
  dgb = W->rank < b->rank ? b->deltas[W->rank] : b->size ;
  
  gMax = b->size / dgb ;
  du1 = 0 * uMax/2 ; du2 = uMax - du1 ; 
  dv1 = 0 * vMax/2 ; dv2 = vMax - dv1 ;

  
  for (g = 0 ; g < gMax ; g++)
    for (i = 0 ; i < iMax ; i++)
      for (j = 0 ; j < jMax ; j++)
	{
	  int u1, u2, v1, v2 ;
	  int uv, bestUV = 0 ;
	  zfb = b->zf + i + j * iMax + g * dgb ;
	  zfa = a->zf + i * uMax + j * vMax * iMax * uMax + g * dga ;
	  zfc = c->zf + i * uMax + j * vMax * iMax * uMax + g * dga ;
	  
	  zf = zfc[0] ; 
	  u1 = du1 < i ? -du1 : -i ;
	  v1 = dv1 < j ? -dv1 : -j ;
	  u2 = du2 < iMax * uMax - i ? du2 : iMax * uMax - i ;
	  v2 = dv2 < jMax * vMax- j ? dv2 : jMax * vMax - j ;

	  for (u = -u1 ; u < u2 ; u++)
	    for (v = -v1 ; v < v2 ; v++)
	      {
		uv = u + v * iMax * uMax ;
		uf = zfc[uv] ;
		if (zf < uf)
		  { zf = uf ; bestUV = uv ; }
	      }
	  zfa[bestUV] = zfb[0] ;
	}
  
  return a ;
} /* mxMaxPoolingBack */

/**********************************************************************/
/* Returns local max in a little square and sequuezes the dimensionality
 * of the first 2 axis (fastest moving)
 */
MX mxMaxPoolingWithCache (MX a, MX b, MX W, MX cache, AC_HANDLE h)
{ 
  int i, iMax ;
  int j, jMax ;
  int g, gMax ;
  
  int u, uMax ;
  int v, vMax ;

  int dc1, dc2, dc3 ;
  int dga, dgb, dgc ;
  float zf, uf, *Restrict zfa, *Restrict zfb ;
  int *Restrict zic ;
  
  mxCheck (a, "mxMaxPooling (a, ..)") ;
  mxCheck (b, "mxMaxPooling (.., b, ..)") ;
  mxCheck (W, "mxMaxPooling (.., b, ..)") ;
  mxCheck (cache, "mxMaxPooling (.., cache)") ;
  mxCheckSameRank (a, b, "mxMaxPooling (a, b, ...)") ;
  mxCheckSameRank (a, cache, "mxMaxPooling (a, .., cache)") ;
  mxCheckSameType (a, b, "mxMaxPooling (a, b, ...)") ;
  
  uMax = W->shapes[0] ;
  vMax = W->rank > 1 ? W->shapes[1] : 1 ;
  /* wMax = W->shapes[2] ;*/
  
  if (a->rank < W->rank)
    messcrash (" mxMaxPooling ((a=%s)->rank) = %d < ((W=%s)->rank = %d)"
	       , a->name, a->rank
	       , W->name, W->rank
	       ) ;
  
  if (cache->rank != a->rank+W->rank)
    messcrash ("rank error") ;
  
  if (W->rank == 1)
    {
      if (cache->shapes[0] != a->shapes[0])
	messcrash ("mxMaxPooling since (W=%s)->rank = 1, but ((cache=%s)->shape[0]=%d) != ((a=%s)->shape[0]=%d)"
		   , cache->name, cache->shapes[0], a->name, a->shapes[0]
		   ) ;
      if (cache->shapes[1] != b->shapes[0])
	messcrash ("mxMaxPooling since (W=%s)->rank = 1, but ((cache=%s)->shape[1]=%d) != ((b=%s)->shape[0]=%d)"
		   , cache->name, cache->shapes[1], b->name, b->shapes[0]
		   ) ;
    }
  else if (W->rank == 2)
    {
      if (cache->shapes[0] != a->shapes[0])
	messcrash ("mxMaxPooling since (W=%s)->rank = 2, but ((cache=%s)->shape[0]=%d) != ((a=%s)->shape[0]=%d)"
		   , cache->name, cache->shapes[0], a->name, a->shapes[0]
		   ) ;
      if (cache->shapes[1] != a->shapes[1])
	messcrash ("mxMaxPooling since (W=%s)->rank = 2, but ((cache=%s)->shape[1]=%d) != ((a=%s)->shape[1]=%d)"
		   , cache->name, cache->shapes[1], a->name, a->shapes[1]
		   ) ;
      if (cache->shapes[2] != b->shapes[0])
	messcrash ("mxMaxPooling since (W=%s)->rank = 2, but ((cache=%s)->shape[2]=%d) != ((b=%s)->shape[0]=%d)"
		   , cache->name, cache->shapes[2], b->name, a->shapes[0]
		   ) ;
      if (cache->shapes[3] != b->shapes[1])
	messcrash ("mxMaxPooling since (W=%s)->rank = 2, but ((cache=%s)->shape[3]=%d) != ((b=%s)->shape[1]=%d)"
		   , cache->name, cache->shapes[3], b->name, b->shapes[1]
		   ) ;
    }
  for (i = 0 ; i < a->rank - W->rank ; i++)
    if (a->shapes [W->rank + i] != b->shapes[2*W->rank + i])
      messcrash ("mxMaxPooling: the spectator shapes of a and b and W should be equal") ;
  
  
  
  if (b->type != MX_FLOAT)
    messcrash (" mxMaxPooling only accepts MX_FLOAT matrices") ;
  if (W->type != MX_BOOL)
    messcrash (" mxMaxPooling only accepts MX_FLOAT W") ;
  if (cache->type != MX_BOOL)
    messcrash (" mxMaxPooling only accepts MX_BOOL cache") ;
  if (W->rank > 2)
    messcrash ("mxMaxPooling called with rank (W=%s) = %d > 2", W->name, W->rank) ;
  
  if (a->shapes[0] * uMax != b->shapes[0])
    messcrash (" mxMaxPooling ((a=%s)shapes[0]=%d)* (widht0=%d) != ((b=%s)shapes[0]=%d)" 
	       , a->name, a->shapes[0], uMax
	       , b->name, b->shapes[0]
	       ) ;
  if (vMax && a->shapes[1] * vMax != b->shapes[1])
    messcrash (" mxMaxPooling ((a=%s)shapes[1]=%d)* (widht1=%d) != ((b=%s)shapes[1]=%d)" 
	       , a->name, a->shapes[1], vMax
	       , b->name, b->shapes[1]
	       ) ;
  iMax = a->shapes[0] ;
  jMax = vMax ? a->shapes[1] : 1 ;
  dc1 = cache->deltas[1] ;
  dc2 = cache->deltas[2] ;
  dc3 = cache->deltas[3] ;

  dga = a->deltas[W->rank] ;
  dgb = b->deltas[W->rank] ;
  dgc = cache->deltas[2 * W->rank] ;
  gMax = b->size / dgb ;
 
  
  for (g = 0 ; g < gMax ; g++)
    for (i = 0 ; i < iMax - uMax + 1 ; i++)
      for (j = 0 ; j < jMax - vMax + 1 ; j++)
	{
	  int bestUV = 0 ;
	  
	  zfa = a->zf + i + j * iMax + g * dga ;
	  zfb = b->zf + i * uMax + j * vMax * iMax + g * dgb ;
	  zic = cache->zi + i + j * dc1 + i * uMax * dc2 + j * vMax * dc3 + g * dgc ;
	  
	  zf = zfb[0] ;
	  
	  for (u = 0 ; u < uMax ; u++)
	    for (v = 0 ; v < vMax ; v++)
	      {
		int uv = u * dc2 + v * dc3 ;
		zic[uv] = 0 ;
		uf = zfb[u + v * iMax] ;
		if (zf < uf)
		  { zf = uf ; bestUV = uv ; }
		zfa[0] = zf ;
		zic[bestUV] = 1 ;
	      }
	}
  
  return a ;
} /* mxMaxPoolingWithCache */

/**********************************************************************/

MX  mxMaxPoolingPullback (MX dX, MX dZ, MX cache)
{
  return dX ;
} /* mxMaxPoolingPullback */

/**********************************************************************/
/**********************************************************************/

static void iRotate (MX a, MX b, MX cache)
{
  int i ;
  float complex *zc, myI = _Complex_I ;
  a = mxRetype (a, b, MX_COMPLEX, 0) ;
  zc = a->zc ;
  for (i = 0 ; i < a->size ; i++, zc++)
    *zc = myI * (*zc) ;
} /* iRotate */

/**********************************************************************/
/* neural network activation functions */
/* Python finction
def sigmoid (Z):
    A = 1 / (1 + np.exp(-Z))
    cache = A * (1-A)
    return A, cache


def complex_sigmoid (Z):
    ZR = Z.real
    K = (ZR > -30) * (ZR < 30)
    U = np.exp(-(ZR))
    A = 1.0 / (1 + U)
    # print ("Sigmoid X=", X, "\nz=", z)
    cache = (A) * (1.0 - A) /2.0
    return A, cache

*/

static void mxSigmoidDo (MX af, MX a, MX cache)
{
  int i, iMax = a->size ;
  int *Restrict ipa, *Restrict ipb, *Restrict ipc ;
  float *Restrict fpa, *Restrict fpb, *Restrict fpc ;
  float complex *Restrict cpa, *Restrict cpb, *Restrict cpc ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      ipa = af->zi ;
      ipb = a->zi ;
      ipc = cache->zi ;
      for (i = 0 ; i < iMax ; i++)
	{
	  ipa[i] = 1 / (1 + exp(-(ipb[i]))) ; 
	  ipc[i] = ipa[i] * ( 1 - ipa[i]) ;
	}
      break ; 
    case MX_FLOAT: 
      fpa = af->zf ;
      fpb = a->zf ;
      fpc = cache->zf ;
      for (i = 0 ; i < iMax ; i++)
	{
	  fpa[i] = 1 / (1 + exp(-(fpb[i]))) ; 
	  fpc[i] = fpa[i] * ( 1 - fpa[i]) ;
	}
       break ;
    case MX_COMPLEX:
      cpa = af->zc ;
      cpb = a->zc ;
      cpc = cache->zc ;
      for (i = 0 ; i < iMax ; i++)
	{
	  double z = creal (cpb[i]) ;
	  cpa[i] = 1 / (1 + exp(-z)) ;
	  cpc[i] = cpa[i] * ( 1 - cpa[i])/2.0 ;
	}

       break ;
    }	
} /* mxSigmoidDo */

/**********************************************************************/
/* neural network activation functions */

static void mxTanhDo (MX af, MX a, MX cache)
{
  int i, iMax = a->size ;
  int *Restrict ipa, *Restrict ipb, *Restrict ipc ;
  float *Restrict fpa, *Restrict fpb, *Restrict fpc ;
  float complex *Restrict cpa, *Restrict cpb, *Restrict cpc ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      ipa = af->zi ;
      ipb = a->zi ;
      ipc = cache->zi ;
      for (i = 0 ; i < iMax ; i++)
	{
	  double z = ipb[i] ;
	  if (z > 0)
	    {
	      double z1 = exp (-z) ;
	      double z2 = z1 * z1 ;
	      double z3 = (1 - z2) / (1 + z2) ;
	      ipa[i] = z3 ; 
	      ipc[i] = 1 - z3 * z3 ;
	    }
	  else
	    {
	      double z1 = exp (z) ;
	      double z2 = z1 * z1 ;
	      double z3 = (z2 - 1) / (1 + z2) ;
	      ipa[i] = z3 ;
	      ipc[i] = 1 - z3 * z3 ;
	    }
	}
      break ; 
    case MX_FLOAT: 
      fpa = af->zf ;
      fpb = a->zf ;
      fpc = cache->zf ;
      for (i = 0 ; i < iMax ; i++)
	{
	  double z = fpb[i] ;
	  if (z > 0)
	    {
	      double z1 = exp (-z) ;
	      double z2 = z1 * z1 ;
	      double z3 = (1 - z2) / (1 + z2) ;
	      fpa[i] = z3 ; 
	      fpc[i] = 1 - z3 * z3 ;
	    }
	  else
	    {
	      double z1 = exp (z) ;
	      double z2 = z1 * z1 ;
	      double z3 = (z2 - 1) / (1 + z2) ;
	      fpa[i] = z3 ;
	      fpc[i] = 1 - z3 * z3 ;
	    }
	}
       break ;
    case MX_COMPLEX:
      cpa = af->zc ;
      cpb = a->zc ;
      cpc = cache->zc ;
      for (i = 0 ; i < iMax ; i++)
	{
	  double z = creal (cpb[i]) ;
	  double z3 = (exp(z) - exp(-z)) / (exp(z) + exp(-z)) ;
	  cpa[i] = z3 ; 
	  cpc[i] = 1 - z3 * z3 ;
	}

       break ;
    }	
} /* mxTanhDo */

/**********************************************************************/

static void mxSoftMaxDo (MX a, MX b, MX c)
{
  int i, iMax, j, jMax, *Restrict ipb ;
  float fmax, fsum, *Restrict fpa, *Restrict fpb , *Restrict fpc ;
  float complex csum, *Restrict cpa, *Restrict cpb, *Restrict cpc ;

  iMax = b->shapes[0] ;
  jMax = b->size / iMax ;

  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      for (j = 0 ; j < jMax ; j++)
	{   
	  fpa = a->zf + j * iMax ;
	  ipb = b->zi + j * iMax ;
	  fpc = c->zf + j * iMax ;

	  fmax = ipb[0] ;
	  for (i = 1 ; i < iMax ; i++)
	    if (fmax < ipb[i])
	      fmax = ipb[i] ;
	  for (fsum = 0, i = 0 ; i < iMax ; i++)
	    {
	      fpa[i] = exp (ipb[i] - fmax) ;
	      fsum += fpa[i] ;
	    }
	  for (i = 0 ; i < iMax ; i++)
	    {
	      fpa[i] /= fsum ;
	      fpc[i] = fpa[i] * ( 1 - (fpa[i])) ;
	    }
	}
       break ;
    case MX_FLOAT: 
      for (j = 0 ; j < jMax ; j++)
	{   
	  fpa = a->zf + j * iMax ;
	  fpb = b->zf + j * iMax  ;
	  fpc = c->zf + j * iMax ;
	  fmax = fpb[0] ;
	  for (i = 1 ; i < iMax ; i++)
	    if (fmax < fpb[i])
	  fmax = fpb[i] ;
	  for (fsum = 0, i = 0 ; i < iMax ; i++)
	    {
	      fpa[i] = exp (fpb[i] - fmax) ;
	      fsum += fpa[i] ;
	    }
	  for (i = 0 ; i < iMax ; i++)
	    {
	      fpa[i] /= fsum ;
	      fpc[i] = fpa[i] * ( 1 - fpa[i]) ;
	    }
	}
      break ;
    case MX_COMPLEX:
       for (j = 0 ; j < jMax ; j++)
	{   
	  cpa = a->zc + j * iMax  ;
	  cpb = b->zc + j * iMax   ;
	  cpc = c->zc + j * iMax   ;
	  fmax = creal (cpb[0]) ;
	  for (i = 1 ; i < iMax ; i++)
	    if (fmax < creal (cpb[i]))
	      fmax = creal (cpb[i]) ;
	  for (csum = 0, i = 0 ; i < iMax ; i++)
	    {
	      cpa[i] = exp (cpb[i] - fmax) ;
	      csum += cpa[i] ;
	    }
	  for (i = 0 ; i < iMax ; i++)
	    {
	      cpa[i] /= csum ;
	      cpc[i] = cpa[i] * ( 1 - (cpa[i])) ;
	    }
	}
       break ;
    }	
} /* mxSoftMaxDo */

/**********************************************************************/
/* RELU */
/*
#@jit
def relu (X):
    cache = (X>0)
    Z = X * cache

    return Z, cache

def complex_relu (X):
    cache = (X.real > 0)
    Z = X * cache

    return Z, cache
*/
static void mxReluDo (MX af, MX a, MX cache)
{
  int i = a->size, *ipa, *ipb, *ipc ;
  float *fpa, *fpb, *fpc ;
  float complex *cpa, *cpb, *cpc ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      ipa = af->zi ;
      ipb = a->zi ;
      ipc = cache->zi ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (ipb[i] > 0)
	    {
	      ipa[i] = ipb[i] ;
	      ipc[i] = 1.0 ;
	    }
	  else
	    {
	      ipa[i] = ipc[i] = 0 ;
	    }
	}
      break ; 
    case MX_FLOAT: 
      fpa = af->zf  ;
      fpb = a->zf  ;
      fpc = cache->zf ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (fpb[i] > 0)
	    {
	      fpa[i] = fpb[i] ;
	      fpc[i] = 1.0 ;
	    }
	  else
	    {
	      fpa[i] = fpc[i] = 0 ;
	    }
	}
       break ;
    case MX_COMPLEX:
      cpa = af->zc ;
      cpb = a->zc ;
      cpc = cache->zc ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (creal (cpb[i]) > 0)
	    {
	      cpa[i] = cpb[i] ;
	      cpc[i] = 1.0 ;
	    }
	  else
	    {
	      cpa[i] = cpc[i] = 0 ;
	    }
	}
      break ; 
     }	
} /*  nnReluDo */

/**********************************************************************/
/* ELU */

static void mxEluDo (MX af, MX a, MX cache)
{
  int i = a->size ; 
  float *fpa, *fpb, *fpc ;
  float complex *cpa, *cpb, *cpc ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      break ; 
    case MX_FLOAT: 
      fpa = af->zf  ;
      fpb = a->zf  ;
      fpc = cache->zf ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (fpb[i] > 0)
	    {
	      fpa[i] = fpb[i] ;
	      fpc[i] = 1.0 ;
	    }
	  else
	    {
	      fpc[i] = exp (fpb[i]) ;
	      fpa[i] = fpc[i] - 1 ;
	    }
	}
       break ;
    case MX_COMPLEX:
      cpa = af->zc ;
      cpb = a->zc ;
      cpc = cache->zc ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (creal (cpb[i]) > 0)
	    {
	      cpa[i] = cpb[i] ;
	      cpc[i] = 1.0 ;
	    }
	  else
	    {
	      cpc[i] = exp (cpb[i]) ;
	      cpa[i] = cpc[i] - 1 ;
	    }
	}
      break ; 
     }	
} /*  nnEluDo */

/**********************************************************************/
/* SOFTRELU */

static void mxSoftReluDo (MX af, MX a, MX cache)
{
  int i = a->size, *ipa, *ipb, *ipc ;
  float *fpa, *fpb, *fpc ;
  float complex *cpa, *cpb, *cpc ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      ipa = af->zi ;
      ipb = a->zi ;
      ipc = cache->zi ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (ipb[i] > 0)
	    {
	      ipa[i] = ipb[i] ;
	      ipc[i] = 1.0 ;
	    }
	  else
	    {
	      ipa[i] = 0 ; /* ipa[i] = .10 * ipb[i] ; but */
	      ipc[i] = 0 ; /* this makes no sense, but .10 is not an int */
	    }
	}
      break ; 
    case MX_FLOAT: 
      fpa = af->zf  ;
      fpb = a->zf  ;
      fpc = cache->zf ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (fpb[i] > 10)
	    {
	      fpa[i] = 10 ;
	      fpc[i] = 0 ;
	    }
	  else if (fpb[i] > 0)
	    {
	      fpa[i] = fpb[i] ;
	      fpc[i] = 1.0 ;
	    }
	  else
	    {
	      fpa[i] = 0 * .10 * fpb[i] ;
	      fpc[i] = 0 * .10 ;
	    }
	}
       break ;
    case MX_COMPLEX:
      cpa = af->zc ;
      cpb = a->zc ;
      cpc = cache->zc ;
      for (i = 0 ; i < a->size ; i++)
	{
	  if (creal (cpb[i]) > 0)
	    {
	      cpa[i] = cpb[i] ;
	      cpc[i] = 1.0 ;
	    }
	  else
	    {
	      cpa[i] = .10 * cpb[i] ;
	      cpc[i] = 0.10 ;
	    }
	}
      break ; 
     }	
} /*  nnSoftReluDo */

/**********************************************************************/

static void mxLogDo (MX af, MX a, MX cache)
{
  int i = a->size, *ipb ;
  float *fpa, *fpb ;
  float complex *cpb ;
  double zf ;
  double min = 1E-12 ;
  double minLog = log (min) ;
  int sign = 1 ;

  if (cache == (MX)2)
    sign = -1 ;
  switch (a->type)
    {
    case MX_NULL:
    case MX_BOOL:
    case MX_INT:
      fpa = af->zf ;
      ipb = a->zi ;
      for (i = 0 ; i < a->size ; i++)
	{
	  zf = ipb[i] ;
	  fpa[i] = sign * (zf > min ? log (zf) : minLog) ;
	}
      break ; 
    case MX_FLOAT: 
      fpa = af->zf ;
      fpb = a->zf ;
      for (i = 0 ; i < a->size ; i++)
	{ 
	  zf = fpb[i] ;
	  fpa[i] =  sign * (zf > min ? log (zf) : minLog) ;
	}
       break ;
    case MX_COMPLEX:
      fpa = af->zf ;
      cpb = a->zc ;
      for (i = 0 ; i < a->size ; i++)
	{
	  zf = creal (cpb[i]) ;
	  fpa[i] =  sign * (zf > min ? log (zf) : minLog) ;
	}
       break ;
    }	
} /* mxLogDo */

/**********************************************************************/
/* Map a function f on every element of a matrix */
MX mxMap (MX af, MX a, MX cache, void (f)(MX af, MX a, MX cache), AC_HANDLE h) 
{
   mxCheck (a, "mxMap (..,a,..)") ;
   if (cache && cache != ((MX) 1) && cache != ((MX) 2)) mxCheck (cache, "mxMap (..,cache,..)") ;
   if (! af)
     af = mxDuplicate ("map", af, a, 0, h) ;
   if (h && ! cache)
     cache = mxDuplicate ("cache", cache, a, 0, h) ;
   mxCheck (af, "mxMap (af,..)") ;
   mxCheckSameShapes (af, a, "mxMap") ;
   if (cache && cache != ((MX) 1) && cache != ((MX) 2)) mxCheckSameShapes (af, cache, "mxMap") ;
   
   f (af, a, cache) ;
   return af ;
} /* mxFilter */

/**********************************************************************/
/**********************************************************************/
/* API to the activation functions */

MX mxSoftMax (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxSoftMaxDo, h) ;
  return af ;
} /* mxSoftMax */

MX mxSigmoid (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxSigmoidDo, h) ;
  return af ;
} /* mxSigmoid */

MX mxRelu (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxReluDo, h) ;
  return af ;
} /* mxRelu */

MX mxElu (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxEluDo, h) ;
  return af ;
} /* mxElu */

MX mxL3elu (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxEluDo, h) ;
  return af ;
} /* mxL3elu */

MX mxTanh (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxTanhDo, h) ;
  return af ;
} /* mxTanh */

MX mxSoftRelu (MX af, MX a, MX cache, AC_HANDLE h)
{
  mxMap (af, a, cache, mxSoftReluDo, h) ;
  return af ;
} /* mxRelu */

MX  mxLog (MX af, MX a, BOOL isMinus, AC_HANDLE h)
{
  mxCheck (af, "mxLog (af,..)") ;
  if (af->type != MX_FLOAT)
    messcrash ("mxLog(af : %s) af->type should be MX_FLOAT, a type can be chosen freely"
	       , af->name) ;
  mxMap (af, a, (MX)(isMinus ? (void*)2 : (void*)1), mxLogDo, h) ;
  return af ;
} /* mxLog */

/**********************************************************************/
/**********************************************************************/


/* test the different matrix operations */
BOOL mxTest (void)
{
  AC_HANDLE h = ac_new_handle () ;
  MX a, b, c, v, w, s ;
  MX af, bf, cf, ac, bc, cc, cc1 ;
  int ii[] = {1,2,3,4,5,6,7,8,9} ;
  int ii10[] = {10,20,30,40,50,60,70,80,90} ;
  /* int id[] = {1,0,0,0,1,0,0,0,1} ; */
  int jj[] = {1,1,1,1,2,2,2,2,3} ;
  /*
    freeinit () ; 
    messErrorInit (argv[0]) ;
  */

  /* check creation */
  a = mxCreate (h, "a", MX_INT, 2, 3, 3, 0) ;
  b = mxCreate (h, "b", MX_INT, 3, 3, 0) ;
  c = mxCreate (h, "c", MX_INT, 3, 3, 0) ;
  v = mxCreate (h, "v", MX_INT, 3, 0) ;
  s = mxScalarMatrix ("Deux", 2, h) ;
  mxCheck (a, "mxTestCreation") ;
  mxCheck (b, "mxTestCreation") ;
  mxCheck (c, "mxTestCreation") ;

  mxShow (s) ;
  mxShow (v) ;
  w = 0 ;
  w = mxAdd (w, v, s, h) ;
  mxShow (w) ;

  /* set b, c */
  mxSet (b, ii) ;
  mxSet (c, jj) ;

  mxShow (b) ;
  mxShow (c) ;
  /* check add */
  fprintf (stderr, "\nAdd")  ;
  mxAdd (a, b, c, h) ;
  mxShow (a) ;
  w = mxAdd (w, v, s, 0) ;
  mxShow (w) ;
  mxAdd (a, b, v, 0) ;
  mxShow (a) ;
  a = 0 ;
  a = mxAdd (a, b, s, 0) ;
  mxShow (a) ;
  mxAdd (a, b, s, 0) ;
  mxShow (a) ;
 
  /* chexk multiply */ 
  fprintf (stderr, "\nScalar product") ;
  
  a = 0 ;
  mxShow (w) ;
  a = mxDot (a, w, w, h) ;
  mxShow (a) ;


  fprintf (stderr, "\nMultiply")  ;
  mxSet (b, ii) ;
  mxSet (c, ii10) ;

  a = 0 ;
  a = mxDot (a, b, c, h) ; 
  mxShow (b) ;
  mxShow (c) ;
  mxShow (a) ;

  fprintf (stderr, "\nFloat Multiply")  ;
  af = bf = cf = 0 ;
  bf = mxRetype (bf, b, MX_FLOAT, h) ;
  cf = mxRetype (cf, c, MX_FLOAT, h) ;

  af = mxDot (af, bf, cf, h) ; 
  mxShow (af) ;

  fprintf (stderr, "\nComplex Multiply")  ;
  ac = bc = cc = 0 ;
  bc = mxRetype (bc, b, MX_FLOAT, h) ;
  cc = mxRetype (cc, c, MX_COMPLEX, h) ;

  ac = mxDot (ac, bc, cc, h) ; 
  mxShow (ac) ;

  fprintf (stderr, "\niRotate c")  ;
  cc1 = 0 ;
  cc1 = mxMap (cc1, cc, 0, iRotate, h) ;

 
  /* transpose index 0 and 1 */
  fprintf (stderr, "\nTranspose")  ;
  a = 0 ;
  a = mxTranspose (a, b, 0, 1, h) ;
  mxShow (a) ;
  return TRUE ;
} /* mxTest */

/**********************************************************************/
/**********************************************************************/

/***********************************************************************************************************************************************/
/*********************************************************** Display utilitie ******************************************************************/
/***********************************************************************************************************************************************/

static double  nicePrintFraction (const char *prefix, double x, const char *suffix);

float complex nicePrint (const char *prefix, float complex z)
{
  double a = creal (z) ;
  double b = cimag (z) ;
  if (b == 0)
    a = nicePrintFraction (prefix, a, "") ;
  else if (a == 0)
    b = nicePrintFraction (prefix, b, "i") ;
  else
    {
      printf("(") ;
      a = nicePrintFraction (prefix, a, "") ;
      b = nicePrintFraction (" + ", b, "i") ;
      printf(")") ;
    }
  return a + I * b ;
} /* nicePrint */

/*******************************************************************************************/
static int niceInt (float x)
{
  int k = x >= 0 ? x + .01 : x - .01 ;
  return k ;
} /* niceInt */

/*******************************************************************************************/

static double nicePrintFraction (const char *prefix, double x, const char *suffix)
{
  int i = 1000000 ;
  int s = x >= 0 ? 1 : -1 ;
  
  x *= s ;
  if (niceInt (1000000*x) == 0)
    { x = 0 ; s = 1 ; }
  printf ("%s", prefix) ;
  if (niceInt (1000000*x) == 1000000 * niceInt (x))
    printf ("%.0f", s*x) ;
  else
    for (i = 2 ; i <= 2048 ; i++)
      if (niceInt (10000*i*x) == 10000 * niceInt (i*x))
	{
	  int y = i * x + .1 ;
	  printf ("%d/%d", s*y, i) ;
	  i = 1000000 ;
	}
  if (i < 1000000)
    printf ("%.5f", s*x) ;
  if (x != 0)
    printf ("%s", suffix) ;

  return s * x ;
} /* nicePrintFraction */

/*******************************************************************************************/

static void niceComplexShow (MX a)
{
  int i,j, iMax = a->shapes[0], jMax = a->shapes[1] ;
  const complex float *zc ;

  if (a->rank == 1) { jMax = iMax ; iMax = 1 ;}
  mxValues (a, 0, 0, &zc) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < jMax ; j++)
	nicePrint ("\t", zc[iMax*j + i]) ;
    }
  printf ("\n") ;
  return ;
} /* niceShow */

/*******************************************************************************************/

static void niceFloatShow (MX a)
{
  int i,j, iMax = a->shapes[0], jMax = a->shapes[1] ;
  const float *zf ;

  if (a->rank == 1) { jMax = iMax ; iMax = 1 ;}
  mxValues (a, 0, &zf, 0) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < jMax ; j++)
	nicePrintFraction ("\t", zf[iMax*j + i],"") ;
    }
  printf ("\n") ;
  return ;
} /* niceFloatShow */

/*******************************************************************************************/
static void niceIntShow (MX a)
{
  int i,j, iMax = a->shapes[0], jMax = a->shapes[1] ;
  const int *zi ;

  if (a->rank == 1) { jMax = iMax ; iMax = 1 ;}
  mxValues (a, &zi, 0, 0) ;
  printf ("## %s", a->name) ;
  for (i = 0 ; i < iMax ; i++)
    {
      printf ("\n") ;
      for (j = 0 ; j < jMax ; j++)
	printf ("\t%d", zi[iMax*j + i]) ;
    }
  printf ("\n") ;
  return ;
} /* niceIntShow */

/*******************************************************************************************/

void mxNiceShow (MX a)
{
  switch (a->type)
    {
    case MX_NULL: break ;
    case MX_BOOL:
    case MX_INT: niceIntShow (a) ; break ;
    case MX_FLOAT: niceFloatShow (a) ; break ;
    case MX_COMPLEX: niceComplexShow (a) ; break ;
    }
  return ;
} /* mxNiceShow */

/**********************************************************************/
/************************  end of file ********************************/
/**********************************************************************/
 
 
 
 

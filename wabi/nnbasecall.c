/*  File: nnbasecall.c
 *  Author: Jean Thierry-Mieg (mieg@kaa.crbm.cnrs-mop.fr)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1994
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmba.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Nov 20 14:01 1998 (fw)
 * Created: Fri Dec 23 16:43:01 1994 (mieg)
 *-------------------------------------------------------------------
 */

/* @(#)nnbasecall.c	1.15 5/23/97 */

/*
#define ARRAY_CHECK
#define MALLOC_CHECK
*/

#include "acedb.h"
#include "acembly.h"
#include "basecall.h"
/* use the Net to evaluate a single base */

#define  FENETRE 36
     /* window for the neural net code */

static void   callNeuralEvaluator (float *yy, float *rr)
{ 
  extern int myNetwork(float *in, float *out, int init) ;
  myNetwork (yy, rr, 1) ;
}

char nnBaseCall (LANE *lane, int pos)
{ float rr, nnReply[4] ;
  int x1, xx, i, j, k ;
  Array zz ;
  Read* seq ;
  float yMax ;
  char cc ;
  TRACE *bp[4], y, max ;

  if (!baseCallGetSeq (lane))
    return 0 ;

  xx = arr (lane->basePos, pos, short) ;
  x1 = xx - FENETRE / 2 ; 
  seq = lane->seq ;
                 /* A C G T special order for neural net */
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceC + x1 ;
  bp[2] = seq->traceG + x1 ;
  bp[3] = seq->traceT + x1 ;

  max = 0 ;  /* normalize in window */
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  if (y > max) max = y ;
	}
    }
  if (max <= 0) return 0 ;

  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceC + x1 ;
  bp[2] = seq->traceG + x1 ;
  bp[3] = seq->traceT + x1 ;

  yMax = max ; /* cast to float */
  zz = arrayCreate (4*FENETRE, float) ;
  i = FENETRE ; k = 0 ; /* prepare output */
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  array (zz, k++, float) = y / yMax ;
	}
    }

  
  callNeuralEvaluator (arrp(zz, 0, float), nnReply) ;
  arrayDestroy (zz) ;
  i = 4 ; yMax = 0 ;
  while (i--)
    { rr = nnReply [i] ;
      if (rr > yMax)
	{ yMax = rr ; k = i ; }
    }
  cc = 0 ;
  if (yMax > .35)
    switch (k)
      {
      case 0:
	cc = T_ ;
	break ;
      case 1:
	cc = G_ ;
	break ;
      case 2:
	cc = C_ ;
	break ;
      case 3:
	cc = A_ ;
	break ;
      }
  return cc ;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

#define  FENETRE 36
/*********************************************************/
/*********************************************************/

static void   callNeuralTrainer (FILE *fil, float *yy, float *rr)
{ int i = 60 , j ;
  static int nTraining = 0 ;

  fprintf (fil, "\n#%d \n", ++nTraining) ;
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--) fprintf (fil, " %g ", *yy++) ;
      fprintf (fil, "\n") ;
    }
  fprintf (fil, "#output \n") ;
  i = 4 ;
  while (i--)
    fprintf (fil, " %g ", *rr++) ;
  fprintf (fil, "\n") ;
}

void nnBaseTrain (FILE *fil, LANE *lane, int pos, char cc)
{ float rr[4] ;
  int x1, xx, i, j, k ;
  static Array zz = 0 ;
  Read* seq ;
  float yMax ;
  TRACE *bp[4], y, max ;

  if (!baseCallGetSeq (lane))
    return ;
  if (pos < 0 || pos >= arrayMax (lane->basePos))
    return ;
  xx = arr (lane->basePos, pos, short) ;
  x1 = xx - FENETRE / 2 ; 
  seq = lane->seq ;
                 /* A C G T special order for neural net */
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceC + x1 ;
  bp[2] = seq->traceG + x1 ;
  bp[3] = seq->traceT + x1 ;

  max = 0 ;  /* normalize in window */
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  if (y > max) max = y ;
	}
    }
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceC + x1 ;
  bp[2] = seq->traceG + x1 ;
  bp[3] = seq->traceT + x1 ;

  yMax = max ; /* cast to float */
  zz = arrayReCreate (zz, 4*FENETRE, float) ;
  i = FENETRE ; k = 0 ; /* prepare output */
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  array (zz, k++, float) = y / yMax ;
	}
    }
  
  i = 4 ;
  while (i--)
    rr[i] = 0 ;
  switch (cc)
    {
    case A_:
      rr [3] = 1 ;
      break ;
    case T_:
      rr [0] = 1 ;
      break ;
    case G_:
      rr [1] = 1 ;
      break ;
    case C_:
      rr [2] = 1 ;
      break ;
    default:
      return ;
    }
      
  callNeuralTrainer (fil, arrp(zz, 0, float), rr) ;
}

static FILE *nnTrainInit (int nErr)
{ FILE *fil ;
  int ii = 4 * FENETRE ;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;

  if (!(fil = filqueryopen (dirName, fileName, "pat", "w",
		     "Where do you wish to export the nn training data ?")))
	return 0 ;

  fprintf(fil, "%s\n%s\n\n\n%s",
	  "SNNS pattern definition file V1.4",
	  "generated at Fri Jul 14 12:40:40 1995",
	  "No. of patterns     : ") ;

  if (nErr > -1)
    fprintf (fil, "%d\n", nErr) ;
  else
    fprintf (fil, "  ###### \n") ;

  fprintf(fil, "%s%d\n%s\n\n",
	  "No. of input units  : ", ii,
	  "No. of output units : 4") ;

  return fil ;
}
/********************/

void nnLaneDoTrain (FILE *fil, Array dna, LANE *lane) 
{ Array errArray = 0 ;
  A_ERR *ep ;
  char cc ;
  int i ;

  if (!baseCallGetSeq (lane))
    return ;

  errArray = lane->errArray ;
  
  if (!(arrayMax (errArray)))
    return ;

  i = arrayMax (errArray) ;
  ep = arrp (errArray, 0, A_ERR) - 1 ;

  while (ep++, i--)
    { 
      switch(ep->type)
	{
	case AMBIGUE: 
/* 	case ERREUR:  */
	  if (ep->iLong < 0 || ep->iLong >= arrayMax (dna))
	    break ;
	  cc = arr(dna, ep->iLong, char) & 0x0F ;
	  if (lane->upSequence)
	    cc = complementBase(cc) ;
	  nnBaseTrain (fil, lane, ep->iShort, cc) ;
	  break ;
	default:
	  break ;
	}
    }
}

void nnLaneTrain (LOOK look, LANE *lane) 
{ FILE *fil ;
  Array errArray = 0 ;
  
  errArray = lane->errArray ;
  
  if (!(arrayMax (errArray)))
    return ;
 
  fil = nnTrainInit (-1) ;
  if (!fil) 
    return ;

  nnLaneDoTrain (fil, look->dna, lane) ;
  filclose (fil) ;
}

void nnAssemblyTrain (LOOK look)
{ LANE *lane ;
  int i ;
  FILE *fil ;

  fil = nnTrainInit (-1) ;
  if (!fil) 
    return ;

  if (look->lanes &&
      (i = arrayMax(look->lanes)))
    while (i--)
      { lane = arrp(look->lanes, i, LANE) ;
	if (lane->scf < 1)
	  traceGetLane (look, lane) ;
	if (lane->scf < 2 ||
	    !findBaseCall (lane) ||
	    ! baseCallGet (lane))
	  continue ;
	nnLaneDoTrain (fil, look->dna, lane) ;
      }
  filclose (fil) ;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/

#ifdef JUNK
/* Method 2: Try to export any sucession of letters and count them */

static void   callNeuralTrainer2 (FILE *fil, float *yy, float *rr)
{ int i = 60 , j ;
  static int nTraining = 0 ;

  fprintf (fil, "\n#%d \n", ++nTraining) ;
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--) fprintf (fil, " %g ", *yy++) ;
      fprintf (fil, "\n") ;
    }
  fprintf (fil, "#output \n") ;
  i = 4 ;
  while (i--)
    fprintf (fil, " %g ", *rr++) ;
  fprintf (fil, "\n") ;
}

void nnBaseTrain2 (FILE *fil, LANE *lane, int pos)
{ float rr[4] ;
  int x1, xx, i, j, k, ibc, t ;
  static Array zz = 0 ;
  Read* seq ;
  float yMax ; 
  TRACE *bp[4], y, max ;
  BASECALL *bb ;
  Array bc ;
  char cc ;

  if (!baseCallGetSeq (lane))
    return ;
  if (pos < 0 || pos >= arrayMax (lane->basePos))
    return ;
  xx = arr (lane->basePos, pos, short) ;
 
  bc = lane->baseCall ;
  for (ibc = 0 ;
       ibc < arrayMax(bc) - 1 ; ibc++)  /* will stop sur avant dernier */
    { bb = arrp (bc, ibc, BASECALL) ;
      if ((bb+1)->x > xx)
	break ;
    }
  if (xx - bb->x > (bb + 1)->x - xx)
    bb++ ;

  cc = array (lane->base, pos, char) ;
  switch (cc)
    { 
    case A_: t = 0 ; break ;
    case G_: t = 1 ; break ;
    case C_: t = 2 ; break ;
    case T_: t = 3 ; break ;
    default: return ;
    }


  if (t != bb->t)
    return ;
  xx = bb->x ;

  x1 = xx - 6 ;  /* center at position 6 */
  seq = lane->seq ;
                 /* A C G T special order for neural net */
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceG + x1 ;
  bp[2] = seq->traceC + x1 ;
  bp[3] = seq->traceT + x1 ;

  max = 0 ;  /* normalize in window */
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  if (y > max) max = y ;
	}
    }
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceG + x1 ;
  bp[2] = seq->traceC + x1 ;
  bp[3] = seq->traceT + x1 ;

  yMax = max ; /* cast to float */
  zz = arrayReCreate (zz, 4*FENETRE, float) ;
  i = FENETRE ; k = 0 ; /* prepare output */
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j]++ ;
	  array (zz, k++, float) = y / yMax ;
	}
    }
  
  i = 4 ;
  while (i--)
    rr[i] = 0 ;
  switch (cc)
    {
    case A_:
      rr [3] = 1 ;
      break ;
    case T_:
      rr [0] = 1 ;
      break ;
    case G_:
      rr [1] = 1 ;
      break ;
    case C_:
      rr [2] = 1 ;
      break ;
    default:
      return ;
    }
      
  callNeuralTrainer (fil, arrp(zz, 0, float), rr) ;
}

static FILE *nnTrainInit2 (int nErr)
{ FILE *fil ;
  int ii = 4 * FENETRE ;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;

  if (!(fil = filqueryopen (dirName, fileName, "pat", "w",
		     "Where do you wish to export the nn training data ?")))
	return 0 ;

  fprintf(fil, "%s\n%s\n\n\n%s",
	  "SNNS pattern definition file V1.4",
	  "generated at Fri Jul 14 12:40:40 1995",
	  "No. of patterns     : ") ;

  if (nErr > -1)
    fprintf (fil, "%d\n", nErr) ;
  else
    fprintf (fil, "  ###### \n") ;

  fprintf(fil, "%s%d\n%s\n\n",
	  "No. of input units  : ", ii,
	  "No. of output units : 4") ;

  return fil ;
}
/********************/

void nnLaneDoTrain2 (FILE *fil, Array dna, LANE *lane) 
{ Array errArray = 0 ;
  A_ERR *ep ;
  int i, nErr,  cTop = lane->clipTop, cEnd = lane->clipEnd ;

  if (lane->upSequence) return ;
  if (!baseCallGetSeq (lane))
    return ;
  errArray = lane->errArray ;
  
  if (nErr = arrayMax (errArray))
    { ep = arrp (errArray, 0, A_ERR) ; nErr-- ; }
  else 
    ep = 0 ;

  for (i = cTop ; i < cEnd ; i++)
    { 
      if (ep && i > ep->iShort)
	{  i = ep->iShort + 2 ;
	   if (nErr-- > 0) ep++ ; else ep = 0 ;
	   continue ; 
	 }

      /* else no error I hope */
      nnBaseTrain2 (fil, lane, i) ;
    }
}

#endif

/*********************************/

void nnContigTrain (KEY contig)
{ LANE *lane ;
  static FILE *fil = 0 ;
  OBJ obj ;
  Array aa, dna ;
  KEY dnaKey ;
  int i, x1, x3, ctop, c3 ;

  if (!fil)
    fil = nnTrainInit (-1) ;
  if (!fil) 
    return ;
  if (!contig)
    { filclose (fil) ; fil = 0 ;
      return ;
    }

  obj = bsCreate (contig) ;
  if (!obj ||
      !bsFindTag (obj, _Assembled_from) || 
      !bsGetKey (obj, _DNA, &dnaKey) ||
      !(dna = dnaGet (dnaKey)))
    { bsDestroy (obj) ;
      return ;
    }
  
  messStatus ("TrainNN") ;
  aa = arrayCreate (100, BSunit) ;
  bsFindTag (obj, _Assembled_from) ;
  bsFlatten (obj, 3, aa) ;
  bsDestroy (obj) ;

  lane = (LANE*) messalloc (sizeof (struct LaneStruct)) ;

  for (i = 0 ; i < arrayMax(aa) ; i+= 3)
    { 
      if (messIsInterruptCalled())
	  goto abort ;

      lane->key = arr(aa,i, BSunit).k ;
      lane->x1 = arr(aa,i + 1, BSunit).i - 1;
      lane->x2 = arr(aa,i + 2, BSunit).i ;
      lane->x3 = lane->x2 ;
  
      traceGetLane (0, lane) ;
      if (lane->scf < 2 ||
	  !findBaseCall (lane) ||
	  ! baseCallGet (lane)) ;
      else
	{ x1 = lane->x1 ; ctop = lane->clipTop ;
	  x3 = lane->x3 - 1 ; c3 = lane->clipExtend - 1 ;  
	  x3 = arrayMax(dna) - 1;
	  x3 = lane->x2 ; c3 = lane->clipEnd - 1 ;
	  if (x1 < x3)
	    lane->errArray = dnaAlignCptErreur (dna, lane->dna, 
				 &x1, &x3,
				 &ctop, &c3) ;
	  if (lane->errArray && arrayMax(lane->errArray))
	    nnLaneDoTrain (fil, dna, lane) ;
	  laneDestroy (lane) ;
	}
    }
  
 abort: 
  arrayDestroy (aa) ;
  arrayDestroy (dna) ;
  messfree (lane) ;
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

/* Quality Training */
#define  Q_FENETRE 60

static void   nnExportQualityTrainer (FILE *fil, float *yy, int rate)
{ int i = 60 , j ;
  static int nTraining = 0 ;

  fprintf (fil, "\n#%d \n", ++nTraining) ;
  i = FENETRE ;
  while (i--)
    { j = 4 ;
      while (j--) fprintf (fil, " %g ", *yy++) ;
      fprintf (fil, "\n") ;
    }
  fprintf (fil, "#output \n") ;

  for (i=0; i < rate; i++)
    fprintf (fil, "0 ") ;

  i++ ;
  fprintf (fil, "1 ") ;

  for (; i < 5; i++)
    fprintf (fil, "0 ") ;

  fprintf (fil, "\n\n") ;
}

static void nnDoQualityTrain (FILE* fil, LANE *lane, int pos, int rate) 
{ int x1, xx, i, j, k, step = 2 ;
  static Array zz = 0 ;
  Read* seq ;
  float yMax ;
  TRACE *bp[4], y, max ;

  if (!baseCallGetSeq (lane))
    return ;

  xx = arr (lane->basePos, pos, short) ;
  x1 = xx - Q_FENETRE * step / 2 ; 
  seq = lane->seq ;
                 /* A C G T special order for neural net */
  bp[0] = seq->traceA + x1 ;
  bp[1] = seq->traceC + x1 ;
  bp[2] = seq->traceG + x1 ;
  bp[3] = seq->traceT + x1 ;

  max = lane->laneMax ;  /* normalize in complete trace */
  yMax = max ; /* cast to float */
  zz = arrayReCreate (zz, 4*Q_FENETRE, float) ;
  i = Q_FENETRE ; k = 0 ; /* prepare output */
  while (i--)
    { j = 4 ;
      while (j--)
	{ y = *bp[j] ; bp[j] += step ; 
	  array (zz, k++, float) = y / yMax ;
	}
    }
  
  nnExportQualityTrainer (fil, arrp(zz, 0, float), rate) ;
}

/************/

static FILE *nnQualityTrainInit (int nErr)
{ FILE *fil ;
  int ii = 4 * Q_FENETRE ;
  static char fileName[FIL_BUFFER_SIZE] , dirName[DIR_BUFFER_SIZE] ;

  strcpy (fileName,"quality") ;
  if (!(fil = filqueryopen (dirName, fileName, "pat", "w",
		     "Where do you wish to export the nn training data ?")))
	return 0 ;

  fprintf(fil, "%s\n%s\n\n\n%s",
	  "SNNS pattern definition file V1.4",
	  "generated at Fri Jul 14 12:40:40 1995",
	  "No. of patterns     : ") ;

  if (nErr > -1)
    fprintf (fil, "%d\n", nErr) ;
  else
    fprintf (fil, "  ###### \n") ;

  fprintf(fil, "%s%d\n%s\n\n",
	  "No. of input units  : ", ii,
	  "No. of output units : 5") ;

  return fil ;
}

/***************/

void nnQualityTrain (LANE *lane, int pos, int rate)
{ int i ;
  static FILE *fil ;

  if (!fil)
    fil = nnQualityTrainInit (-1) ;
  if (!fil) 
    return ;

  if (lane->scf < 2 ||
      !findBaseCall (lane) ||
      ! baseCallGet (lane))
    return ;

  i = 10; 
  while (i-- && pos > Q_FENETRE)
    nnDoQualityTrain (fil, lane, pos -= 7, rate) ;
}

/*********************************************************/

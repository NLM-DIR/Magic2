/*  File: randsubs.c
 *  Author: Jean Thierry-Mieg (mieg@mrc-lmb.cam.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1992
 *-------------------------------------------------------------------
 * This file is part of the ACEDB genome database package, written by
 * 	Richard Durbin (MRC LMB, UK) rd@mrc-lmb.cam.ac.uk, and
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * Description: reproducable low quality random numbers
 *
 * Exported functions:
 * HISTORY:
 * Last edited: Dec  4 15:06 1998 (fw)
 * Created: Mon Jun 15 14:44:56 1992 (mieg)
 *-------------------------------------------------------------------
 */

 /* $Id: randsubs.c,v 1.4 2014/05/28 22:29:56 mieg Exp $ */

#ifdef ALPHA
int random(void);		/* in libc.a */
#endif /* ALPHA */

/* despite the statics, i consider this package reentrant
 * if we are mutithreaded, the numbers just become
 * non reproducible
 */

static int xrand = 18721 ;
static int yrand = 37264 ;  
static int zrand = 28737 ;

/*******************************/

double randfloat (void)
{
  double x ;
  
  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323;
  x = xrand/30269.0 + yrand/30307.0 + zrand/30323.0;
  return (x-(int)x);
}

/*******************************/

#if defined(ALPHA) || defined(LINUX) || defined(OPTERON)
long random (void) ;
int randint (void)
{   
  return (int) random() ;
}
#else
int randint (void)
{
  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323; 
  return (zrand) ;
}
#endif

#ifdef JUNK
/* 2024_04: alternative code found on the web, not tested */

KEY randKey (KEY *ks)
{
  KEY x ;
  ks[0] += 1111111111;
  x = ks[0] ^ ((ks[1] << 21) | (ks[1] >> 11));
  ks[1] = ((ks[2] << 13) | (ks[2] >> 19)) + x ;
  ks[2] = ks[0] + x ;
  return ks[2];
}
#endif

/*******************************/

double randgauss (void)
{
  double sum ;
  double fac = 3.0/90899.0 ;

  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323;
  sum = xrand + yrand + zrand ;
  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323;
  sum += xrand + yrand + zrand ;
  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323;
  sum += xrand + yrand + zrand ;
  xrand = 171*xrand % 30269;
  yrand = 172*yrand % 30307;
  zrand = 170*zrand % 30323;
  sum += xrand + yrand + zrand ;
  return (sum*fac - 6.0) ;
 }

/*******************************/

void randsave (int *arr)
{
  arr[0] = xrand ;
  arr[1] = yrand ;
  arr[2] = zrand ;
}

/*********************************/

void randrestore (int *arr)
{
  xrand = arr[0] ;
  yrand = arr[1] ;
  zrand = arr[2] ;
}

/*********************************/
/*********************************/

/*  File: bitset_.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) J Thierry-Mieg and R Durbin, 1997
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
 *	Jean Thierry-Mieg (CRBM du CNRS, France) mieg@kaa.cnrs-mop.fr
 *
 * SCCS: $Id: bitset_.h,v 1.4 2020/05/30 15:03:14 mieg Exp $
 * Description: transformation of array package for bits
 * Exported functions: a pure header - no code
 * HISTORY:
 * * 2015 : this file bitset_.h is a copy of the previous bitset_.h
 *  it is now automatically included by ac.h
 *   codes that explicitley use bit sets, should include bitset_.h 
 * Last edited: Dec  4 14:51 1998 (fw)
 * * Sep 15 09:18 1998 (edgrif): Add a macro to define bitField so that
 *              those that want it can define it (avoids compiler warnings).
 * Created: Sat Feb  1 11:07:23 1997 (rd)
 *-------------------------------------------------------------------
 */

#ifndef BITSET__H
#define BITSET__H

#include "regular.h"
#include "bigarray.h"

/* Some of the below macros operate on a array called 'bitField', to use     */
/* these macros code this macro first at the appropriate level in your       */
/* .c file.                                                                  */

typedef BigArray BitSet ; 

BitSet bitSetCreate (unsigned long int n, AC_HANDLE h) ;
BitSet bitSetReCreate (BitSet bb, unsigned long int n)  ;
void bitExtend (BitSet bb, unsigned long int n) ;
#define bitSetDestroy(_x)	bigArrayDestroy(_x)

#define bitSetCopy(_x, _handle) bigArrayHandleCopy(_x, _handle)

/* the following use array() and so are "safe" */

#define bitSet(_x,_n)	(bigArray((_x), (_n) >> 5, unsigned int) |= \
			 bitField[(_n) & 0x1f])
#define bitUnSet(_x,_n)	(bigArray((_x), (_n) >> 5, unsigned int) &= \
			 ~bitField[(_n) & 0x1f])
#define bitt(_x,_n)	(bigArray((_x), (_n) >> 5, unsigned int) & \
			 bitField[(_n) & 0x1f])

/* bit() uses arr() for optimal performance */

#ifndef MEM_DEBUG
#define bit(_x,_n)	(bigArr((_x), (_n) >> 5, unsigned int) & bitField[(_n) & 0x1f])
#else
#define bit(_x,_n) bitt(_x,_n)
#endif
                  /* Returns the number of set bits */
unsigned long int bitSetCount (BitSet a) ;
      /* return the maximal allocated bit, useful to scan the table */
#define bitSetExists(_x) ((_x) && (_x)->magic == BIG_ARRAY_MAGIC ? (_x)->id : 0 )


           /* performs b1 = b1 OPERATOR b2, returns count (b1) */
unsigned long int bitSetAND (BitSet b1, BitSet b2) ;
unsigned long int bitSetOR (BitSet b1, BitSet b2) ;
unsigned long int bitSetXOR (BitSet b1, BitSet b2) ;
unsigned long int bitSetMINUS (BitSet b1, BitSet b2) ;
unsigned long int bitSetANDcount (BitSet b1, BitSet b2) ;

#endif

/***************** end of file *****************/
 
 

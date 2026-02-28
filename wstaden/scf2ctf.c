#include "Read.h"

#include "regular.h"
#include "freeout.h"

/***************************************************************/
/* missing func */
char *name(KEY key) { return "trace file" ; }
void cDNAAlignInit (void) { return ; }
char* baseCallNameGuess (KEY key, int type) { return 0 ; }
BOOL baseCallGet (void *lane) { return 0 ; } 

/***************************************************************/

int main(int argc, char **argv)
{ 
  Read *rr = 0 ;

  freeinit() ;
  freeOutInit() ;

  if (argc>1)
    goto usage ;

  rr = fread_reading(stdin, "stdin", 0) ; 
  fwrite_reading (stdout, rr, TT_CTF) ; 
  fwrite_reading (stderr, rr, TT_PLN) ; 

  return 0 ;

usage:
  freeOutf("Usage:  scf2ctf < ff.scf > ff.ctf\n") ;
  return 1 ;
}

/***************************************************************/
/***************************************************************/

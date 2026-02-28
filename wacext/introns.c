/*
 * Authors: Jean Thierry-Mieg, NCBI, 
 * Apr 2016
 * Statistical analysis of intron distribution
 *  type, banc de Poisson
*/

#define BITSET_MAKE_BITFIELD   
#include <ac.h>
#include <acedna.h>

typedef struct ppStruct { 
  AC_HANDLE h ; 
  const char *inFileName, *outFileName ;
  BOOL gzi, gzo ;
  BOOL detect, count ;
  const char *run ;
  int minSupport, minOverhang, minFrequency ;
  const char *mrnaRemapFilename ;
  const char *intronListFilename ;
} PP ;

typedef struct ddStruct { int chrom, a1, a2, n ; BOOL isDown ; char type ;} DD ;

/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (char *message)
{
  if (! message)
    fprintf  (stderr,
	      "// intron: Intron analysis\n"
	      "// Authors: Danielle and Jean Thierry-Mieg, NCBI, October 2009, mieg@ncbi.nlm.nih.gov\n"
	      "// Purpose\n"
	      "// Create a set of Introns databases in tmp/Introns_DB/zoner.*\n"
	      "//   Import donor acceptors and introns detected in each run, in genome or transcript\n"
	      "//   Remap to genome and use this as the main identifier\n"
	      "//   Count all detected introns in all runs using the remapping to count the suport only once\n"
	      "//   The idea is to count donors and acceptors, including the non-spliced reference\n"
	      "//   The we can compute a chi-square for each donor/acceptor, giving the raltive usage of the intron in that gene\n"    
	      "// Caveat:\n"
	      "//   Lines starting with '#' are considered as comments and dropped out\n"
	      "//\n"
	      "// File names\n"
	      "//   -i mainInputFileName\n"
	      "//   -o prefix for all output files, ex: dir1/dir2/introns. \n"
	      "//   --gzi : force unzip the input, useful if importing from stdin\n"
	      "//          all input files names .gz are automatically unzipped\n"
	      "//   --gzo : gzip all output files\n"
	      "//\n"
	      "// Usage \n"
	      "//  --detect : detect the introns in the best_hit file\n"
	      "//    -i best_hits_file.gz\n"
	      "//    -r (--run) run (run name, mandatory)\n"
	      "//    --minSupport <int> [default 10] : provide the run intron threshold \n"
	      "//    --minOverhang <int> [default 8] : number of flawless base on both side of donor or acceptor \n"
	      "//    --minFrequency <int> [default 5]: percentage of reads jumping in this way at donor/acceptor sites\n"
	      "//    --mrnaRemap mrnaRemap.gz : position of each intron in the transcript and on the genome\n"
	      "//\n"
	      "//  --count :count the introns and the reference in the best_hit file\n"
	      "//    -i best_hits_file.gz\n"
	      "//    -r run (run name, mandatory)\n"
	      "//    --minSupport <int> [default 10] : provide the run intron threshold \n"
	      "//    --minOverhang <int> [default 8] : number of flawless base on both side of donor or acceptor \n"
	      "//    --minFrequency <int> [default 5]: percentage of reads jumping in this way at donor/acceptor sites\n"
	      "//    --mrnaRemap mrnaRemap.gz : position of each intron in the transcript and on the genome\n"
	      "//    --intronList file : [optional] the default in to count all the introns in the database\n"
	      "//  -h --help : this help\n"
 	      ) ;
  if (message)
    {
      fprintf (stderr, "// %s\nFor more information try:  dna2dna -help\n", message) ;
    }
  exit (1);
  
} /* usage */

/*************************************************************************************/
/*************************************************************************************/

int main (int argc, const char **argv)
{
  PP pp ;
  AC_HANDLE h = 0 ;

  freeinit () ; 
  h = ac_new_handle () ;
  memset (&pp, 0, sizeof (PP)) ;
  pp.h = h ;


  /* optional arguments */

  if (argc == 1 ||
      getCmdLineOption (&argc, argv, "-h", 0) ||
      getCmdLineOption (&argc, argv, "-help", 0) ||
      getCmdLineOption (&argc, argv, "--help", 0)
      )
    usage (0) ;

  getCmdLineOption (&argc, argv, "-i", &pp.inFileName) ;
  getCmdLineOption (&argc, argv, "-o", &pp.outFileName) ;

  pp.gzi = getCmdLineBool (&argc, argv, "-gzi") ;
  pp.gzo = getCmdLineBool (&argc, argv, "-gzo") ;
  pp.gzi |= getCmdLineBool (&argc, argv, "--gzi") ;
  pp.gzo |= getCmdLineBool (&argc, argv, "--gzo") ;

  pp.detect = getCmdLineBool (&argc, argv, "--detect") ;
  pp.count = getCmdLineBool (&argc, argv, "--count") ;

  if (! getCmdLineOption (&argc, argv, "-r", &pp.run))
    getCmdLineOption (&argc, argv, "--run", &pp.run) ;


  getCmdLineOption (&argc, argv, "--mrnaRemap", &pp.mrnaRemapFilename) ;
  getCmdLineOption (&argc, argv, "--intronList", &pp.intronListFilename) ;

  getCmdLineInt (&argc, argv, "--minSupport", &pp.minSupport) ;
  getCmdLineInt (&argc, argv, "--minOverhang", &pp.minOverhang) ;
  getCmdLineInt (&argc, argv, "--minFrequency", &pp.minFrequency) ;

  fprintf (stderr , "Processed %d introns\n", 0) ;
  {
    int mx ;
    messAllocMaxStatus (&mx) ;   
    fprintf (stderr, "// done: %s\tmax memory %d Mb\n", timeShowNow(), mx) ;
  }
  ac_free (h) ;
  if (1) sleep (1) ; /* to prevent a mix up between stdout and stderr */
  return 0 ;
}

/*************************************************************************************/
/*************************************************************************************/
/*************************************************************************************/



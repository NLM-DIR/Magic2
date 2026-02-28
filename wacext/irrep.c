#include "ac.h"
#include "matrix.h"
#include "bitset.h"

#define MALLOC_CHECK
#define ARRAY_CHECK

/* construct the table of representations of the basic classical Lie superalgebras 
 */

#define RMAX 8 

typedef struct wStruct
{
  int x[RMAX] ;
  int yRaw, y, yDenominator ;
  BOOL isBlack ;
  int mult ;
  int crystal ;
  int n21 ;
  int k ;
  int layer ;
  int l2 ;
  BOOL odd, hw ;
  int myOddPart ;
  BitSet oddB ;
} WW ; /* representation weigth vector */

typedef struct saStruct
{
  AC_HANDLE h ;
  
  int method ;
  BOOL table ;
  const char *type ; /* A,B.C,D,F,G */
  const char *DynkinWeights ; /* 1:0:2:.... */
  BOOL hasY ;
  int YY[RMAX], YYd ;
  int m, n, alpha, rank, rank1, rank2, rank3 ;
  Array scale, lowerMetric, upperMetric ;
  Array Cartan, Cartan1, Cartan2, Cartan3 ;
  int detgg ;
  int detGG ;
  int hasAtypic ;
  BOOL hasOdd ;
  BOOL isBlack ;
  int hasExtended ;
  int nEven, nOdd ;
  BOOL odd[RMAX] ;
  BOOL odd1[RMAX] ;
  BOOL odd2[RMAX] ;
  BOOL odd3[RMAX] ;
  BOOL extended[RMAX] ;
  Array Kac ; /* Kac crystal */
  WW evenHw, evenHw1, evenHw2, evenHw3 ;
  WW oddHw ;
  WW extendedHw ;
  Array evenRoots ; /* odd roots */
  Array negativeOddRoots ; /* odd roots */
  Array wws ; /* weights */
  WW hw, rho, rho0, rho1 ;
  DICT *dict ;
  int pass ;
  int D1, D2 ; /* number of even and odd generators of the adjoiunt rep */
  Array dd ;   /* dims of all the submodules */
  Array atypic ;
  Array wwsShifted ;
  MX chi ;
} SA ; /* SuperAlgebra struct */


static int wwScalarProduct (SA *sa, WW *ww1, WW *ww2) ;
static BOOL demazureEvenOdd (SA *sa, Array wws, int r1, BOOL even, int *dimEvenp, int *dimOddp, BOOL show) ;
static int demazure (SA *sa, int *dimEvenp, int *dimOddp, int method, BOOL show) ;

/******************************************************************************************************/

static int wwCreationOrder (const void *a, const void *b)
{
  const WW *up = (WW *)a ;
  const WW *vp = (WW *)b ;
  int n = 0 ;

  n = up->k - vp->k ;
  return n ;  
} /* wwCreationOrder */

/******************************************************************************************************/

static int wwLayerOrder (const void *a, const void *b)
{
  const WW *up = (WW *)a ;
  const WW *vp = (WW *)b ;
  int n = 0 ;

  n = up->layer - vp->layer ; if (n) return n ;
  n = up->k - vp->k ;  if (n) return n ;

  n = up->x[0] - vp->x[0] ;  if (n) return -n ;
  n = up->x[1] - vp->x[1] ;  if (n) return -n ;
  n = up->x[2] - vp->x[2] ;  if (n) return -n ;
  n = up->x[3] - vp->x[3] ;  if (n) return -n ;
  n = up->x[4] - vp->x[4] ;  if (n) return -n ;

  return n ;  
} /* wwLayerOrder */

/******************************************************************************************************/  

static void wwY (SA *sa, WW *ww)
{
  int *YY = sa->YY ;
  int i, y1 = 0, y2 = 0 ;
  
  for (i = 0 ; i < sa->rank ; i++)
    y1 += YY[i] * ww->x[i] ;

  y2 = sa->YYd ;
  if (y2 < 0) {y2 = -y2 ; y1 = -y1 ;}
  ww->yRaw = y1 ;
  if (y1 == 0)
    y2 = 1 ;
  if (y1 == y2)
    y1 = y2 = 1 ;
  if (y1 == -y2)
    { y1 = -1 ; y2 = 1 ; }
  for (i = 2 ; i<= y2 ; i++)
    {
      int z1= y1/i, z2 = y2/i ;
      if (i*z1 == y1 && i*z2 == y2)
	{ y1 = z1 ; y2 = z2 ; i = 1 ; }
    }
  ww->y = y1 ;
  ww->yDenominator = y2 ;
  return  ;
} /* wwY */

/******************************************************************************************************/

static int locateWeight (SA *sa, WW *w, BOOL create)
{
  int kk = 0 ;
  char buf[1024] ;

  /* locate it to construct the multiplicity */
  memset (buf, 0, sizeof (buf)) ;
  sprintf (buf, "%d:%d:%d:%d:%d:%d:%d:%d", w->x[0] , w->x[1] , w->x[2] , w->x[3] , w->x[4] , w->x[5] , w->x[6] , w->x[7] ) ;
  if (create)
    dictAdd (sa->dict, buf, &kk) ;
  else
    dictFind (sa->dict, buf, &kk) ;
  w->l2 = wwScalarProduct (sa, w, w) ;
  wwY (sa, w) ;
  return kk ;
} /* locateWeight */

/******************************************************************************************************/

static void metricRescale (Array metric, int r, int ii, int jj, int scale)
{
  int i, k ;
    for (i = ii ; i <= jj ; i++)
  for (k = 0 ; k < r ; k++)
    {
      int x = arr (metric, r * i + k, int) ;
      arr (metric, r * i + k, int) = scale * x ;
    }
  return ;
} /* metricRescale */

/******************************************************************************************************/
/* merge and complete the Cartan Matrix for a superalgebra */
static void mergeCartan (SA *sa, BOOL show)
{
  int i, j ;
  int r1 = sa->rank1 ;
  int r2 = sa->rank2 ;
  int r3 = sa->rank3 ;
  int r = r1 + r2 + r3 ;
  int rr = r * r ;
  int dx = 0 ;
  BOOL hasY ;
  
  sa->rank = r ;
  sa->Cartan = arrayHandleCreate (rr, int , sa->h) ;

  array (sa->Cartan, rr-1, int) = 0 ;  

  sa->rank = r ;
  hasY = sa->hasY ;

  if (r1 > 0)
    {
      for (i = 0 ; i < r1 ; i++)
	for (j = 0 ; j < r1 ; j++)
	  {
	    array (sa->Cartan, r * i + j, int) = arr (sa->Cartan1, r1 * i + j, int) ;
	  }
      dx = r1 ;
      if (*sa->type == 'A')
	sa->evenHw1.x[dx] = r1 ;
    }
  if (r2 > 0)
    {
      WW oldhw2 = sa->evenHw2 ;


      for (i = 0 ; i < r ; i++)
	sa->evenHw2.x[i] = 0 ;
      for (i = 0 ; i < r2 ; i++)
	{	
	  for (j = 0 ; j < r2 ; j++)
	    {
	      array (sa->Cartan, r * (i+dx) + j + dx, int) = arr (sa->Cartan2, r2 * i + j, int) ;
	    }
	  sa->evenHw2.x[dx+i] = oldhw2.x[i] ;
	}
      dx = r1 + r2 ;
    }
  if (r3 > 0)
    {
      WW oldhw3 = sa->evenHw3 ;

      for (i = 0 ; i < r ; i++)
	sa->evenHw3.x[i] = 0 ;
      for (i = 0 ; i < r3 ; i++)
	{	
	  for (j = 0 ; j < r3 ; j++)
	    {
	      array (sa->Cartan, r * (i+dx) + j + dx, int) = arr (sa->Cartan3, r3 * i + j, int) ;
	    }
	  sa->evenHw3.x[dx+i] = oldhw3.x[i] ;
	}
      if (*sa->type == 'A')
	sa->evenHw3.x[r1] = r3 ;
    }
  if (0 && hasY)
    {
      if (r1) array (sa->Cartan, r * (r1 - 1) + r1, int) = -1 ;
      if (r1) array (sa->Cartan, r * (r1) + r1 - 1, int) = -1 ;
      if (r2) array (sa->Cartan, r * (r1) + r1 + 1, int) = 1 ;
      if (1 && r1 && r2) array (sa->Cartan, r * (r1) + r1 - 1, int) = -1 ;
      if (r2) array (sa->Cartan, r * (r1+1) + r1, int) = -1 ;
      dx++ ;
    }

  if (0 && hasY)
    {
      int i, ro = r1 ;
      int YY[RMAX] ;
      
      memset (YY, 0, sizeof (YY)) ;

      switch ((int)sa->type[0])
	{
	default:
	  messcrash ("Only A and C have a U(1) Y hypercharge") ;
	  break ;
	case 'A':
	case 'C':
	  for (i = 0 ; i < ro ; i++)
	    YY[i] =  -(r2+1)*(i+1) ;
	  for (i = r-1 ; i > ro ; i--)
	    YY[i] =  (r1+1)*(r-i) ;
	  YY[ro]  = (r1+1)*(r2+1) * (r1 <= r2 ? -1 : -1) ;
	  sa->YYd = (r1 == r2 ? 1 : r1 - r2) ;
	  for (i = 0 ; i < r ; i++)
	    sa->YY[i] = YY[i] ;

	  if (show)
	    {
	      printf ("............................................. YYd=%d  YY[]=", sa->YYd) ;
	      for (i = 0 ; i < r ; i++)
		printf (" %d", sa->YY[i]) ;
	      printf ("\n") ;
	    }
	  break ;
	}
    }
} /* mergeCartan */

/******************************************************************************************************/
/* prepare the Cartan Matrix for a simple Lie algebra block */
static void getOneCartan (SA *sa, char type, int r, int Lie, BOOL show)
{
  int i, r2 = r * r ;
  WW hw ;
  Array Cartan = arrayHandleCreate (r2, int, sa->h) ;

  array (Cartan, r2 - 1, int) = 0 ;
  
  memset (&hw, 0, sizeof (hw)) ;
  switch ((int)type)
    {
    case 'A':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      
      hw.x[0] = 1 ;
      hw.x[r-1] += 1 ;
      break ;
      
    case 'B':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r * 0 + 1, int) = -2 ;
      hw.x[r-2] = 1 ;
      break ;
      
    case 'C':           
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*r - r - 1, int) = -2 ;

      hw.x[0] = 2 ;
      break ;
      
    case 'D':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 1) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i> 0 && i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      if (r > 2)
	{
	  array (Cartan, r*0 + 2, int) = -1 ;
	  array (Cartan, r*2 + 0, int) = -1 ;
	}
      if (r == 2)
	{
	  hw.x[0] = 2 ;
	  sa->evenHw2.x[1] = 2 ;
	}
      else if (r == 3)
	{
	  hw.x[0] = 1 ;
	  hw.x[1] = 1 ;
	}
      else
	hw.x[r-2] = 1 ;
      break ;
      
    case 'E':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, (r-2) * r + r-3, int) = 0 ;
      array (Cartan, (r-3) * r + r-2, int) = 0 ;
      array (Cartan, (r-2) * r + r-4, int) = -1 ;
      array (Cartan, (r-4) * r + r-2, int) = -1 ;

      switch (r)
	{
	case 6: hw.x[3] = 1 ; break ;
	case 7: hw.x[6] = 1 ; break ;
	case 8: hw.x[0] = 1 ; break ;
	}
      break ;
      
    case 'F':
      for (i = 0 ; i < r ; i++)
	{
	  array (Cartan, r*i + i, int) = 2 ;
	  if (i > 0) array (Cartan, r*i + i-1, int) = -1 ;
	  if (i < r-1) array (Cartan, r*i + i+1, int) = -1 ;
	}
      array (Cartan, r*2 + 1, int) = -2 ;

      hw.x[0] = 1 ;
      break ;
      
    case 'G':
      array (Cartan, r*0 + 0, int) = 2 ;
      array (Cartan, r*0 + 1, int) = -1 ;
      array (Cartan, r*1 + 1, int) = 2 ;
      array (Cartan, r*1 + 0, int) = -3 ;
      
      hw.x[0] = 1 ;
      break ;
      
    case 'U':
      sa->hasY = TRUE ;
      break ;
    default:
      messcrash ("The Cartan matrix of a Lie algebra of type %c is not yet programmed, sorry.", type) ;
    }


  switch (Lie)
    {
    case 0:
      sa->Cartan = Cartan ;
      sa->rank = r ;
      sa->evenHw = hw ;
      break ;
    case 1:
      sa->Cartan1 = Cartan ;
      sa->rank1 = r ;
      sa->evenHw1 = hw ;
      break ;
    case 2:
      sa->Cartan2 = Cartan ;
      sa->rank2 = r ;
      sa->evenHw2 = hw ;
      break ;
    case 3:
      sa->Cartan3 = Cartan ;
      sa->rank3 = r ;
      sa->evenHw3 = hw ;
      break ;
    }
  
  if (show)
    {
      int i, j ;
      printf ("\n### getOneCartan (lower) Matrix type %c rank = %d\n", type, r) ;
      for (i = 0 ; i < r ; i++)
	{
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (Cartan, r*i + j, int)) ;
	  printf ("\n") ;
	}
    }
  
} /* getOneCartan */

/******************************************************************************************************/

static void getCartan (SA *sa, BOOL show)
{
  int m = sa->m, n = sa->n ;
  
  if (sa->n == -99999)
    {
      n = sa->n = 0 ;
      sa->rank = m ; 
      sa->hasOdd = FALSE ;
    }
  else
    sa->hasOdd = TRUE ;
  
  switch ((int)sa->type[0])
    {
    case 'A':
      if (! sa->hasOdd)
	getOneCartan (sa, 'A', m, 0, TRUE) ;
      else
	{
	  int r = m + n + 1 ;
	  if (m >= 1) getOneCartan (sa, 'A', m, 1, TRUE) ;
	  if (1)      getOneCartan (sa, 'U', 1, 2, TRUE) ;
	  if (n >= 1) getOneCartan (sa, 'A', n, 3, TRUE) ;
	  mergeCartan (sa, show) ;
	  if (m > 0)
	    {
	      array (sa->Cartan, (m-1)*r + m, int) = -1 ;
	      array (sa->Cartan, (m)*r + m-1, int) = -1 ;
	    }
	  if (n > 0)
	    {
	      array (sa->Cartan, (m)*r + m+1, int) = 1 ;
	      array (sa->Cartan, (m+1)*r + m, int) = -1 ;
	    }
	  else if (n == 0)
	    array (sa->Cartan, (m)*r + m-1, int) = 1 ;
	  
	  if (m>0) sa->oddHw.x[m-1] = 1 ;
	  if (n>0) sa->oddHw.x[m+1] = 1 ;
	  sa->odd[m] = 1 ;
	  sa->odd1[m] = TRUE ;
	  if (sa->method >= 10)
	    sa->odd2[m] = TRUE ;
	}
      break ;
      
    case 'B':
      if (! sa->hasOdd)
	{
	  if (m == 1)
	    getOneCartan (sa, 'A', 1, 0, TRUE) ;
	  else if (m == 2)
	    getOneCartan (sa, 'C', 2, 0, TRUE) ;
	  else
	    getOneCartan (sa, 'B', m, 0, TRUE) ;
	}
      else
	{
	  if (m <= 1)
	    messcrash ("Type Lie B(m) and Kac OSp(1/2m) cannot have m=%d < 2\n", m) ;
	  if (n == 0)
	    getOneCartan (sa, 'B', m, 0, TRUE) ;
	  else
	    {
	      sa->n = 0 ;
	      sa->extended[0] = TRUE ;
	      sa->hasOdd = TRUE ;
	      sa->isBlack = TRUE ;
	      getOneCartan (sa, 'C', m, 0, TRUE) ;
	      
	      sa->odd[0] = 1 ;
	      sa->oddHw.x[m-1] = 1 ;
	      sa->oddHw.x[0] = -1 ;
	    }
	}
      break ;
      
    case 'C':
      if (! sa->hasOdd)
	{
	  if (m == 1)
	    messcrash ("C(1)=sp(2)=su(2)=A(1), please request A(1)") ;
	  else
	    getOneCartan (sa, 'C', m, 0, TRUE) ;
	}
      else
	{
	  sa->rank = n + 1 ;
	  sa->n = 0 ;
	  sa->hasOdd = TRUE ;
	  getOneCartan (sa, 'U', 1, 1, TRUE) ;
	  getOneCartan (sa, 'C', n, 2, TRUE) ; 
	  mergeCartan (sa, show) ;
	  sa->oddHw.x[n-1] = 1 ;
	  if (sa->method >= 10)
	    {
	      sa->odd1[0] = TRUE ;
	      sa->odd2[1] = TRUE ;
	    }
	}
      break ;

    case 'D':   /* D(m/n) = OSp(2m/2n):SO(2m)+Sp(2n):D(m)+C(n) */
      if (! sa->hasOdd)
	{
	  if (m == 1)
	    messcrash ("D(1)=so(2)=u(1)is Abelian") ;
	  else if (m == 2)
	    {
	      getOneCartan (sa, 'A', 1, 1, TRUE) ;
	      getOneCartan (sa, 'A', 1, 2, TRUE) ;
	      mergeCartan (sa, show) ;
	    }
	  else
	    getOneCartan (sa, 'D', m, 0, TRUE) ;
	}
      else
	{
	  if (m == 1)
	    messcrash ("Type Kac D(1,n) = OSp(2/2n) = C(1,n), please request C(1,n)") ;
	  else if (m == 2)
	    {
	      getOneCartan (sa, 'A', 1, 1, TRUE) ;
	      getOneCartan (sa, 'A', 1, 2, TRUE) ;
	    }
	  else 
	    getOneCartan (sa, 'D', m, 1, TRUE) ;
	  if (n == 1)
	    getOneCartan (sa, 'A', 1, 3, TRUE) ;
	  else
	    getOneCartan (sa, 'C', n, 3, TRUE) ;

	  mergeCartan  (sa, show) ;
	  sa->extended[m] = TRUE ;
	  sa->hasOdd = TRUE ;
	  sa->odd[2] = 1 ;
	  if (sa->method >= 1000)
	    {
	      sa->odd1[m] = TRUE ;
	      sa->odd2[m] = TRUE ;
	    }
	  if (m == 2)
	    {
	      sa->oddHw.x[0] = 1 ;
	      sa->oddHw.x[1] = 1 ;
	    }
	  else
	    sa->oddHw.x[m-1] = 1 ;
	  sa->oddHw.x[m] = -1 ;
	}
      break ;
      
    case 'E':
      if (sa->hasOdd)
	messcrash ("superalgebra of type E(m,n) is not defined") ;
      switch (m)
	{
	case 6:  /* Lie algebra E6 */
	case 7:     /* Lie algebra E7 */
	case 8:     /* Lie algebra E8 */
	  getOneCartan (sa, 'E', m, 0, TRUE) ;
	  break ;
	default:
	  messcrash ("Type Lie E(m) should be E(6 or 7 or 8) not E(%d)", m) ;
	  break ;
	}
      break ;
      
     case 'F':
       if (!sa->hasOdd)
	 {
	   if (m == 4)
	    getOneCartan (sa, 'F', 4, 0, TRUE) ;
	   else
	    messcrash ("Type Lie F(m) should be F(4) not F(%d)", m) ;
	 }
       else
	 {
	   if (m != 1 || n != 3)
	     messcrash ("Type Kac F(m,n) should F(1,3) not F(%d,%d)", m, n) ;
	   sa->rank = sa->m = 4 ; sa->n = 0 ;
	   sa->extended[3] = TRUE ;
	   sa->hasOdd = TRUE ;
	   sa->oddHw.x[2] = 1 ;
	   sa->oddHw.x[3] = -1 ;
	   getOneCartan (sa, 'B', 3, 1, TRUE) ;
	   getOneCartan (sa, 'A', 1, 2, TRUE) ;
	   mergeCartan  (sa, show) ;
	 }
       break ;
       
    case 'G':
       if (!sa->hasOdd)
	 {
	   if (m == 2)
	    getOneCartan (sa, 'G', 2, 0, TRUE) ;
	   else
	    messcrash ("Type Lie G(m) should be G(2) not G(%d)", m) ;
	 }
       else
	 {	  /* Lie superalgebra G(3) */
	   if (m != 2 || n != 1)
	     messcrash ("Type Kac G(m,n) should G(2,1) not G(%d,%d)", m, n) ;
	   sa->extended[2] = TRUE ;
	   sa->hasOdd = TRUE ;
	   sa->odd[2] = 1 ;
	   if (sa->method >= 10)
	     {
	       sa->odd1[1] = TRUE ;
	       sa->odd2[2] = TRUE ;
	     }
	   sa->oddHw.x[1] = 1 ;
	   sa->oddHw.x[2] = -1 ;
	   getOneCartan (sa, 'G', 2, 1, TRUE) ;
	   getOneCartan (sa, 'A', 1, 2, TRUE) ;
	   mergeCartan  (sa, show) ;
	 }
       break ;
      
    default:
      messcrash ("The Cartan matrix of a superalgebra of type %c is not yet programmed, sorry.", sa->type) ;
    }
  
  if (sa->hasOdd || show)
    {
      int r = sa->rank ;
      printf ("\n### Cartan (lower) Matrix type %s m=%d n=%d rank = %d\n", sa->type, m, n, r) ;
      for (int i = 0 ; i < r ; i++)
	{
	  for (int j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (sa->Cartan, r*i + j, int)) ;
	  printf ("\n") ;
	}
    }
}   /* KK */
/* getCartan */

/******************************************************************************************************/
/* the lower Cartan metric is the symmetrized form of the Cartan matrix
*  the upper metric is the inverse of the lower metric
* the scale is the metric relative to cartan
*/
static void getMetric (SA *sa, BOOL show)
{
  Array aa ;
  int i, j ;
  int m = sa->m, n = sa->n, r = sa->rank ;

  aa = sa->lowerMetric = arrayHandleCopy (sa->Cartan, sa->h) ;
  sa->upperMetric = arrayHandleCreate (r*r, int , sa->h) ;
  sa->scale = arrayHandleCreate (r, int , sa->h) ;
  array (sa->scale, r - 1 , int) = 1 ;

  for (i = 1 ; i < r ; i++)
    {
      int x = - arr (aa, r * (i-1) + i , int) ;
      int y = - arr (aa, r * (i+0) + i - 1, int) ;

      if (x * y < 0 )
	{
	  y = -y ;
	  metricRescale (aa, r, i, r - 1, -1) ;
	}
      if (x > y)
	  metricRescale (aa, r, i, r - 1, x/y) ;
      else if (x < y)
	metricRescale (aa, r, 0, i - 1, y/x) ;
      }

  if (sa->hasOdd)
    switch ((int)sa->type[0])
      {
      case 'D':
	for (i = m ; i < r ; i++)
	  metricRescale (aa, r, r - 1, r - 1, -2) ;
	break ;
      }

  for (i = 0 ; i < r ; i++)
    {
      for (j = 0 ; j < r ; j++)
	{
	   int x = - arr (aa, r * i + j , int) ;
	   int y = - arr (sa->Cartan, r * i + j, int) ;
	   if (x * y != 0)
	     {
	       arr (sa->scale, i , int) = x / y ;
	       break ;
	     }
	}
    }

    
  sa->detgg = mxIntDeterminant (arrp(sa->lowerMetric, 0, int), r) ;
  if (show && sa->lowerMetric)
      {
	
	printf ("\n### Lower Metric type %s m=%d n=%d rank = %d det=%d\n", sa->type, m, n, r, sa->detgg) ;
	for (i = 0 ; i < r ; i++)
	  {
	    for (j = 0 ; j < r ; j++)
	      printf ("\t%d", arr (sa->lowerMetric, r*i + j, int)) ;
	    printf ("\n") ;
	  }
      }
    sa->detGG = mxIntInverse (arrp (sa->upperMetric, 0, int),  arrp(sa->lowerMetric, 0, int), sa->rank) ;
    if (show && sa->upperMetric)
      {

	printf ("\n### Upper Metric type %s m=%d n=%d rank = %d, det=%d\n", sa->type, m, n, r, sa->detGG) ;
	for (i = 0 ; i < r ; i++)
	{
	  for (j = 0 ; j < r ; j++)
	    printf ("\t%d", arr (sa->upperMetric, r*i + j, int)) ;
	  printf ("\n") ;
	}
      }

    if (show && sa->scale)
      {

	printf ("\n### Scale = metric/Cartan\n") ;
	for (i = 0 ; i < r ; i++)
	  printf ("\t%d", arr (sa->scale, i, int)) ;
	printf ("\n") ;
      }
}  /* getMetric */

/******************************************************************************************************/  
#ifdef JUNK
static int wwNaturalProduct (SA *sa, WW *ww1, WW *ww2)
{
  int i, x = 0, rank = sa->rank ;
  
  for (i = 0 ; i < rank ; i++)
    x += ww1->x[i] * ww2->x[i] ;
  
  return x ;
} /* wwNaturalProduct */
#endif
/******************************************************************************************************/  

static int wwScalarProduct (SA *sa, WW *ww1, WW *ww2)
{
  Array scale = sa->scale ;
  Array upperMetric = sa->upperMetric ;
  int i, j, x = 0, rank = sa->rank ;

  for (i = 0 ; i < rank ; i++)
    for (j = 0 ; j < rank ; j++)
      x += ww1->x[i] * arr (upperMetric, rank*i+j,int) * ww2->x[j] * arr (scale, i, int) * arr (scale, j, int) ;
  
  return x ;
} /* wwScalarProduct */

/******************************************************************************************************/  

static void wwsShow (SA *sa, char *title, int type, Array wws, WW *hw)
{
  if (wws)  
    {
      int dimE = 0, dimO = 0 ;
      int ii, iMax = arrayMax (wws) ;
      int r, rank = sa->rank ;

      printf ("\n##################### %s: t=%d : hw", title, type) ;
      for (r = 0 ; hw && r < rank ; r++)
	printf (" %3d", hw->x[r]) ;
      printf ("\n") ;

      for (ii = 0 ; ii < iMax ; ii++)
	{
	  WW *ww = arrp (wws, ii, WW) ;
	  if (ww->mult)
	    {
	      locateWeight (sa, ww, FALSE) ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %3d", ww->x[r]) ;
	      printf ("\tmult=%2d n21=%2d k=%2d l=%3d %5s %s Lsquare=", ww->mult, ww->n21, ww->k, ww->layer, ww->odd ? "Odd" : "", ww->hw ? "*" : "") ;
	      if (ww->layer < 0)
		messcrash ("Negative layer") ; 
	      ww->l2 = wwScalarProduct (sa, ww, ww) ;
	      if (sa->detGG && ww->l2 % sa->detGG == 0)
		printf ("%d", ww->l2/sa->detGG) ;
	      else
		printf ("%d/%d", ww->l2, sa->detGG) ;
	      if (sa->YYd)
		{
		  printf ("\tYraw=%d",ww->yRaw) ;
		  if (ww->yDenominator <= 1) printf ("\tY=%d", ww->y) ;
		  else printf ("\tY=%d/%d", ww->y, ww->yDenominator) ;
		}
	      if (sa->hasOdd && ww->oddB)
		{
		  int i ;
		  char *sep = "\todd=" ;
		  
		  for (i = 0 ; i < bitSetCount (ww->oddB) ; i++)
		    if (bit (ww->oddB, i))
		      {
			printf ("%s%d", sep, i+1) ;
			sep = "," ;
		      }
		}
	      printf ("\n") ;
	      if (ww->odd)
		dimO += ww->mult ;
	      else
		dimE += ww->mult ;
	    }
	}
      printf ("## %s dimEven=%d dimOodd=%d dim=%d sdim=%d\n", title, dimE, dimO, dimE+dimO, dimE - dimO) ;
    }
  return ;
} /* wwsShow */

/******************************************************************************************************/

static WW getHighestWeight (SA *sa, int type, BOOL create, BOOL show)
{
  WW *hw ;
  Array wws ;
  int k = 0 ;
  int rank = sa->rank ;
  int r ;

  /* reinitialize a single hw.w */
  sa->pass = 0 ;

  wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
  hw = arrayp (wws, 1, WW) ;
  hw[-1].layer = -1000 ;
  
  switch (type)
    {
    case 0: /* trivial h.w. */
      break ;
    case -1: /* use as h.w. the lowering simple root */
      hw->odd = TRUE ;
      for (r = 0 ; r < rank ; r++)
	hw->x[r] = sa->oddHw.x[r] ;
      break ;
    case -2: /* construct from the parameters */
      if (sa->wws)
	k = sscanf (sa->DynkinWeights, "%d:%d:%d:%d:%d:%d:%d:%d", &hw->x[0] , &hw ->x[1] , &hw ->x[2] , &hw ->x[3] , &hw ->x[4] , &hw ->x[5] , &hw ->x[6] , &hw ->x[7] ) ;
      else
	messcrash ("Missing argument --s, expectirn -wws 1:0:2;...    giving the Dynkin lables of the highest weight ") ;

      if (k != sa->rank)
	messcrash ("HW: rank = %d but you provided %d Dynkin weights\n", sa->rank, k) ;
      break ;
    case -3:
      *hw = sa->evenHw ;
      if (sa->rank1)
	{
	  WW eew, *ew ;
	  *hw = sa->evenHw1 ;
	  memcpy (&eew, &sa->evenHw1, sizeof (WW)) ; ;
	  eew.mult = create ? 1 : 0 ;
	  eew.hw = TRUE ;
	  eew.k = locateWeight (sa, &eew, TRUE) ;
	  ew = arrayp (wws, eew.k, WW) ;
	  *ew = eew ;
	}	  
      if (sa->rank2)
	{
	  WW eew, *ew ;
	  memcpy (&eew, &sa->evenHw2, sizeof (WW)) ; ;
	  eew.mult = create ? 1 : 0 ;
	  eew.hw = TRUE ;
	  eew.k = locateWeight (sa, &eew, TRUE) ;
	  ew = arrayp (wws, eew.k, WW) ;
	  *ew = eew ;
	}
      if (sa->rank3)
	{
	  WW eew, *ew ;
	  memcpy (&eew, &sa->evenHw3, sizeof (WW)) ; ;
	  eew.mult = create ? 1 : 0 ;
	  eew.hw = TRUE ;
	  eew.k = locateWeight (sa, &eew, TRUE) ;
	  ew = arrayp (wws, eew.k, WW) ;
	  *ew = eew ;
	}
      break ;
    default:
      if (sa->negativeOddRoots && type <= arrayMax (sa->negativeOddRoots))
	{
	  WW *wo = arrp (sa->negativeOddRoots, type, WW) ;
	  for (r = 0 ; r < rank ; r++)
	    hw->x[r] = sa->hw.x[r] - wo->x[r] ;
	}
      else
	messcrash ("getHighestWeight received type=%d out of range of odd roots", type) ;
      break ;
    }
  hw->mult = create ? 1 : 0 ;
  hw->hw = TRUE ;
  hw->k = locateWeight (sa, hw, TRUE) ;

  if (type == -2)
    {
      sa->hw = *hw ;
      if (show)
	{
	  printf ("####### Highest weight ") ;
	  for (r = 0 ; r < rank ; r++)
	    printf (" %d", hw->x[r]) ;
	  printf ("\n") ;
	}
    }
  return *hw ;
} /* getHighestWeight */

/******************************************************************************************************/
/* construct the KacCrystal */
static void getKacCrystal (SA *sa, BOOL show)
{
  Array wws ;
  AC_HANDLE h = ac_new_handle () ;
  Array negativeOddRoots = sa->negativeOddRoots ;
  int r, rank = sa->rank ;
  int k, kMax = sa->nOdd ;
  int ii, iiMax = 1 ;
  WW *w, hw ;
  int method = sa->method ;
  int myOddPart ;
  BitSet oddB ;
  
  /* reinitialize on the trivial h.w */
  sa->dict = dictHandleCreate (32, sa->h) ;
  getHighestWeight (sa, 0, 0, show) ;
  wws = sa->wws ;
  hw = arr (wws, 1, WW) ;

  for (k = 0 ; k < kMax ; k++)
    iiMax *= 2 ;   /* 2^negativeOddRoots = size of the KacCrystal */
  oddB = bitSetCreate (kMax, h) ;
  for (ii = 0 ; ii < iiMax ; ii++) /* all points in the Kac Crystal, h.w. already included */
    {
      int layer = 0 ;
      int x = ii ;
      BOOL ok = TRUE ;
      vTXT txt = vtxtHandleCreate (0) ;
      WW w0 ;
      int n21 = 0 ;
      
      memset (&w0, 0, sizeof (w0)) ;
      w0.hw = 1 ;
      myOddPart = 0 ;
      oddB = bitSetReCreate (oddB, kMax) ;
      for (k = 0 ; k < kMax ; k++)
	{
	  int yes  = x % 2 ; /* odd root k is used in ii */
	  x /= 2 ;
	  if (yes == 1)
	    {
	      WW *wodd = arrayp (negativeOddRoots, k + 1, WW) ;
	      wodd->l2 = wwScalarProduct (sa, wodd, wodd) ;

	      if (method % 10 == 2 && arr(sa->atypic,k+1,int) > 0)
		n21 ++ ;
	      if (sa->isBlack || wodd->l2 == 0)
		{
		  bitSet (oddB, k) ;
		  vtxtPrintf (txt, " %d", k+1) ;
		  layer++ ;
		  if (layer == 1)
		    myOddPart = k + 1 ;
		  for (r = 0 ; r < rank ; r++)
		    w0.x[r] += wodd->x[r] ;
		}
	      else
		ok = FALSE ;
	    }
	}
      if (0)
	printf ("....getKacCrystal k=%d ok=%d iiMax=%d kMax=%d\n", ii, ok, iiMax, kMax) ;

      if (ok)
	{
	  int k2 = locateWeight (sa, &w0, TRUE) ;
	  w = arrayp (wws, k2, WW) ; /* new node */
	  w->k = locateWeight (sa, &w0, TRUE)  ;

	  w->hw = TRUE  ;
	  w->layer = layer ;
	  w->odd = layer %2 ;
	  if (layer == n21) /* fully atypic */
	    w->n21 = layer ;
	  if (! n21)
	    w->mult++ ;
	  if (layer && w->n21 == w->layer)  /* kill fully atypic and his friends */
	    w->mult = 0 ;
	  if (layer == 1)
	    w->myOddPart = myOddPart ;
	  w->oddB = bitSetCopy (oddB, sa->h) ;
	  for (r = 0 ; ok && r < rank ; r++)
	    w->x[r] = w0.x[r] ;
	  if (w->mult && show)
	    {
	      printf (".....Kac crystal: ") ;
	      for (r = 0 ; r < rank ; r++)
		printf (" %d", w->x[r]) ;
	      printf ("  :: ii=%d  yes:%s :: ", ii, vtxtPtr (txt)) ; 
	      printf ("\tmult=%d n21=%d k=%d l=%d %s %s\n", w->mult, w->n21, w->k, w->layer, w->odd ? "Odd" : "", w->hw ? "*" : "" ) ;
	    }
	}
      ac_free (txt) ;
    }
  arraySort (sa->wws, wwLayerOrder) ;
  sa->Kac = sa->wws ;
  sa->wws = 0 ;
  ac_free (h) ;
  if (show)
    wwsShow (sa, "Kac crystal dim =", 1+arrayMax(sa->Kac), sa->Kac, &hw) ;
  if (0)
    printf ("# Constructed %d Kac weigthst\n", arrayMax (sa->Kac) - 1) ;
} /* getKacCrystal */

/******************************************************************************************************/

static Array getShiftedCrystal (SA *sa, int atypic, WW *hwAp, BOOL show)
{
  AC_HANDLE h = ac_new_handle () ;
  int rank = sa->rank ;
  int ii, kk, iMax ;
  Array wws ;
  Array Kac = sa->Kac ;
  int kMax = arrayMax (Kac) ;
  WW hwA ;
  BitSet oddB = bitSetCreate (kMax, h) ;
  
  hwA = getHighestWeight (sa, atypic, TRUE, 0) ;
  *hwAp = hwA ;
  /* reinitialize wws and but not the dictionary */
  wws = sa->wws = arrayHandleCreate (256, WW, sa->h) ;
  
  for (kk = 1 ; kk < kMax ; kk++)
    {
      WW *wo = arrp (Kac, kk, WW) ;
      WW w, *ww ;
      int r, k2 ;
	    
      if (! wo->mult) continue ;
      /* set the corrds of the next point of the tensor product */
      memset (&w, 0, sizeof (w)) ;
      w = hwA ;
      for (r = 0 ; r < rank ; r++)
	w.x[r] += wo->x[r] ;
      
      /* locate it to construct the multiplicity */
      k2 = locateWeight (sa, &w, TRUE)  ;
      ww = arrayp (wws, k2, WW) ; /* new node */
      ww->k = k2 ;
      ww->layer = wo->layer ;
      ww->odd = wo->odd ;
      ww->oddB = wo->oddB ;
      ww->myOddPart = wo->myOddPart ;
      ww->mult += hwA.mult * wo->mult ;
      for (r = 0 ; r < rank ; r++)
	ww->x[r] = w.x[r] ;
    }
  /* remove the non-extended  negative even weights */
  iMax = arrayMax (wws) ;
  if (0) wwsShow (sa, "Shifted Kac Crystal", atypic, wws, &hwA) ;
  for (ii = 2 ; ii < iMax ; ii++)
    {
      int r ;
      WW *ww = arrp (wws, ii, WW) ;
      if (! ww->mult) continue ;
      for (r = 0 ; r < rank ; r++)
	if (! sa->odd[r] && /* ! sa->extended[r] && */ ww->x[r] < 0)
	  ww->mult = 0 ;
      if (ww->mult > 0 && ww->layer == 1)
	bitSet (oddB, ww->myOddPart - 1) ;
    }
  /* trim the nodes like beta_1 beta_4 of su(2/3)
   * such that the h.w. dynkin weight beta_4 - beta_3 = 0
   * because they cannot be principal weights since 
   *  alpha-bar beta_4 = beta_3 and alpha_bar lambda = 0
   * test on the small su(2/3) reps
   */
  if (iMax)
    {
      int r, jj ;
      for (r = 0 ; r < rank ; r++)
	if (! sa->odd[r] && sa->hw.x[r] == 0)
	  {
	    for (ii = 2 ; ii < iMax ; ii++)
	      {
		WW *w1 = arrp (wws, ii, WW) ;

		if (w1->mult)
		  {
		    for (jj = ii+1 ; jj < iMax ; jj++)
		      {
			WW *w2 = arrp (wws, jj, WW) ;
			
			if (ii != jj && w2->mult && w1->layer == w2->layer)
			  {
			    int ok = TRUE ;
			    int i ;
			    for (i = 0 ; ok && i < rank ; i++)
			      if (w2->x[i] + array (sa->Cartan, rank*r + i, int) != w1->x[i])
				ok = FALSE ;
			    if (ok)
			      w2->mult = 0 ;
			  }
		      }
		  }
	      }
	  }
    }
  
  /* trim the branches not connected to the root */
  switch (sa->method)
    {
    case 2:
    case 12:
      for (ii = 1 ; ii < iMax ; ii++)
	{
	  int i, n1 ;
	  WW *ww = arrp (wws, ii, WW) ;
	  BOOL ok = ii == 1 ? TRUE : FALSE ;
	  
	  if (! ww->mult) continue ;
	  n1 = bitSetCount (ww->oddB) ;
	  for (i = 1 ; !ok && i < iMax ; i++)
	    {
	      WW *w2 = arrp (wws, i, WW) ;
	      if (w2->layer == ww->layer - 1)
		{
		  int k, nk = 0 ;
		  for (k = 0 ; k < kMax ; k++)
		    if (w2->mult > 0 && bit (ww->oddB, k) && bit (w2->oddB, k))
		      nk++ ;
		  if (nk == n1 - 1)
		    ok = TRUE ;
		}
	    }
	  if (! ok)
	    ww->mult = 0 ;
    }
      break ;
    }
  /* symmetrize the extended  even weights */
  iMax = arrayMax (wws) ;
  if (0)
    {
      int r ;
      for (r = 0 ; r < rank ; r++)
	if (sa->extended[r])
	  demazureEvenOdd (sa, wws, r, TRUE, 0, 0, 0) ; 
    }

  if (show)
    wwsShow (sa, "Shifted Kac Crystal before extended Demazure", atypic, wws, &hwA) ;

  if (1)
    for (ii = 1 ; ii < iMax ; ii++)
      {
	int r, rank = sa->rank ;
	WW *ww = arrp (wws, ii, WW) ;
	if (! ww->mult) continue ;
	for (r = 0 ; r < rank ; r++)
	  if (sa->extended[r] && ww->x[r] > 0)
	    {
	      int r2, k, dk = ww->x[r] ;
	      WW w2 ;
	      int k2, mult2 = 0, dm = 0 ;

	      memset (&w2, 0, sizeof (w2)) ;
	      for (r2 = 0 ; r2 < rank ; r2++)
		w2.x[r2] = ww->x[r2] ;
	      w2.x[r] -= 2 * dk ;
	      k2 = locateWeight (sa, &w2, TRUE) ;
	      mult2 = w2.mult ;
	      dm = ww->mult - mult2 ;
	      if (dm > 0)
		for (k = 1 ; k <= dk ; k++)
		  {
		    WW *ww2 ;

		    memset (&w2, 0, sizeof (w2)) ;
		    for (r2 = 0 ; r2 < rank ; r2++)
		      w2.x[r2] = ww->x[r2] ;

		    w2.odd = ww->odd ;
		    w2.x[r] -= 2 * k ;
		    k2 = locateWeight (sa, &w2, TRUE) ;
		    w2.k = k2 ;
		    ww2  = arrayp (wws, k2, WW) ;
		    if (!ww2->k)
		      *ww2 = w2 ;
		    if (ww2->x[r] < 0)
		      ww2->mult += dm ;
		  }
	    }
      }

  if (show)
    wwsShow (sa, "Shifted Kac Crystal", atypic, wws, &hwA) ;
  ac_free (h) ;
  return wws ;
} /* getShiftedCrystal */

/******************************************************************************************************/

static void intersectCrystals (SA *sa, Array wws, Array xxs)
{
  WW *ww, *xx ;
  int j ;
  
  wwsShow (sa, "INTERSECT wws", 0, wws, 0) ;
  wwsShow (sa, "INTERSECT xxs", 1, xxs, 0) ;
  for (j = 0 ; j < arrayMax (wws) ; j++)
    { /* we use a single dict */
      ww = arrp (wws, j, WW) ;
      xx = arrayp (xxs, j, WW) ;
      if (xx->mult < ww->mult)
	ww->mult = xx->mult ;
    }
  wwsShow (sa, "INTERSECT SUBTRACTED wws", 0, wws, 0) ;
  return ;
} /* intersectCrystals */

/******************************************************************************************************/

static void getHwCrystal (SA *sa, int *dimp, int *sdimp,  BOOL show)
{
  WW hwA ;
  Array wws = getShiftedCrystal (sa, -2, &hwA, TRUE) ;
  Array old = sa->wws ;
  int method = sa->method ;
  
  if (method == 1 && sa->hasAtypic)
    {
      int ii ;
      
      for (ii = 0 ; 1 && ii < arrayMax (sa->atypic) ; ii++)
	if (array (sa->atypic, ii, int))
	  {
	    Array xxs = getShiftedCrystal (sa, ii, &hwA, show) ;
	    intersectCrystals (sa, old, xxs) ;
	    ac_free (xxs) ;
	    
	    if (show)
	      wwsShow (sa, "Get Subtracted Atypical Tensor Product", ii, old, &hwA) ;
	  }
      wws = old ;
    }
  sa->wws = wws ;
  if (dimp)
    {
      int dimE = 0, dimO = 0 ;
      WW *ww ;
      int ii ;
      for (ii = 0 ; ii < arrayMax (wws) ; ii++)
	{
	  ww = arrp (wws, ii, WW) ;
	  if (ww->mult > 0)
	    {
	      if (ww->odd)
		dimO += ww->mult ;
	      else
		dimE += ww->mult ; 
	    }
	}
      *dimp = dimE + dimO ;
      *sdimp = dimE - dimO ;
    }
  arraySort (wws, wwLayerOrder) ;
  if (1)
    {
      Array old = wws ;
      WW *ww, *ww1 ;
      int ii ;

      wws = sa->wws = arrayHandleCreate (64, WW, sa->h) ;
      /* clean up the dict and redeclare the locators */
      sa->dict = dictHandleCreate (32, sa->h) ;
      for (ii = 0 ; ii < arrayMax (old) ; ii++)
	{
	  ww = arrp (old, ii, WW) ;
	  if (ww->mult > 0)
	    {
	      int k = locateWeight (sa, ww, TRUE) ;
	      ww1 = arrayp (wws, k, WW) ;
	      *ww1 = *ww ;
	      ww1->k = k ;
	    }
	}
    }
  if (show)
    wwsShow (sa, "Final Atypical Tensor Product", 999, wws, 0) ;

  return ;
} /* getHwCrystal */

/******************************************************************************************************/

static BOOL demazureEvenOdd (SA *sa, Array wws, int r1, BOOL even, int *dimEvenp, int *dimOddp, BOOL show) 
{
  BOOL new = FALSE ;
  BOOL odd = ! even ;
  int ii, dimEven, dimOdd, r, rank = sa->rank, ra = -1 ;
  
  if (sa->pass++ == 0) new = TRUE ;
    
  if (sa->hasOdd)
    for (r = 0 ; r < sa->rank ; r++)
      {
	if (sa->odd2[r])  /* distinguished even root */
	  ra = r ;
      }

  for (ii = 0 ; ii < arrayMax (wws) ; ii++) 
    {
      WW *w0 = arrayp (wws, ii, WW) ;
      int n1 = w0->mult ;
      int n2 = 0 ;
      int dn = 0 ;
      int jMax = w0->x[r1] ; /* default even values */

      if (ii == 0)
	{ w0->mult = 0 ; continue ; }
      if (! even && jMax != 0)        /* jMax = b = 0 is atypic 1 */
	jMax = 1 ;
      if (jMax  <= 0)
	continue ;
      if (even)
	n1 = w0->mult ;
      else
	n1 = w0->mult - w0->n21 ;
      if (n1 <= 0)
	continue ;
      if (jMax > 0)  /* create or check existence of the new weights along the sl(2) submodule */
	{
	  int j, r, k2 ;
	  WW ww2, *w2 ;

	  memset (&ww2, 0, sizeof(WW)) ;
          for (j = 1 ; j <= jMax ; j++) 
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		ww2.x[r] = w0->x[r] - array(sa->Cartan, rank * r + r1, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      k2 = locateWeight (sa, &ww2, TRUE)  ;
	      
	      /* create w2(k2)  if needed */
	      w2 = arrayp (wws, k2, WW) ;
	      w0 = arrayp (wws, ii, WW) ;  /* needed because the wws array may be relocated in RAM upon extension */ 
	      if (! w2->k)
		{
		  *w2 = *w0 ;
		  w2->k = k2 ;
		  w2->layer = (even ? w0->layer : w0->layer + 1) ;
		  w2->odd = even ? w0->odd : ! w0->odd ;
		  w2->mult = 0 ;
		  w2->n21 = 0 ;
		  w2->hw = FALSE ;
		  /* w2->l2 = wwScalarProduct (sa, w2, w2) ; done in locate */
		}
	      for (r = 0 ; r < rank ; r++)
		w2->x[r] = ww2.x[r] ;
	      n2 = w2->mult ;
	    }
	}

      dn = n1 - n2 ; /* defect multiplicity */
      if (jMax > 0 && dn > 0)   /* populate the new multiplet */
	{
	  int j, r, k2 ;
	  WW ww2, *w2 ;
	  
	  new = TRUE ;
	  memset (&ww2, 0, sizeof(WW)) ;
	  for (j = jMax ; j >= 1 ; j--)
	    {
	      /* position to the new weight */
	      for (r = 0 ; r < rank ; r++)
		ww2.x[r] = w0->x[r] - array(sa->Cartan, rank * r + r1, int) * j ;
	      
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      k2 = locateWeight (sa, &ww2, TRUE) ;
	      w2 = arrayp (wws, k2, WW) ;
	      w0 = arrayp (wws, ii, WW) ;  /* needed because the wws array may be relocated in RAM upon extension */ 
	      /* increase multiplicity of the new multiplet */
	      if (even)
		{
		  w2->mult += dn ;
		  if (sa->method == 12 && w2->mult > w2->crystal)
		    w2->mult = w2->crystal ;
		}
	      if (odd)
		{
		  w0->n21 += dn ;
		  w2->n21 += dn ;
		  ra = ra + 0 ;
		}
	    }
	}
      
    }

  for (dimEven = dimOdd = 0, ii = 1 ; ii < arrayMax (wws) ; ii++)
    {
      WW *ww = arrp (wws, ii, WW) ;
      if (ww->odd)
	dimOdd += ww->mult ;
      else 
	dimEven += ww->mult ;
    }

  arraySort (wws, wwCreationOrder) ;
  
  if (dimEvenp) *dimEvenp = dimEven ;
  if (dimOddp) *dimOddp = dimOdd ; 
  
  return new ;
} /* demazureEvenOdd */

/******************************************************************************************************/

static Array getSU21Crystal (SA *sa, WW *w0, int ra, int rb, AC_HANDLE h)
{
  int jMax = w0->x[ra] + 1 ;
  int rank = sa->rank ;
  Array aa21 = arrayHandleCreate (4 * jMax + 1, WW, h) ;
  int a = w0->x[ra] ;
  int b = w0->x[rb] ;
  int s = - array (sa->Cartan, sa->rank * rb + ra, int) ;
  int jj = 0 ;
  int j, r ;
  WW *w21 ;
      
  if (a >= 0)
    for (j = 0 ; j <= a ; j++)
      {
	w21 = arrayp (aa21, jj++, WW) ;
	/* position to the new weight */
	for (r = 0 ; r < rank ; r++)
	  w21->x[r] = w0->x[r] - array(sa->Cartan, rank * r + ra, int) * j ;
	w21->k = locateWeight (sa, w21, TRUE)  ;
	w21->odd = FALSE ;
	w21->layer = 0 ;
      }
  if (b != 0)
    for (j = 0 ; j <= a + 1 ; j++)
      {
	w21 = arrayp (aa21, jj++, WW) ;
	/* position to the new weight */
	for (r = 0 ; r < rank ; r++)
	  w21->x[r] = w0->x[r] 
	    - array(sa->Cartan, rank * r + ra, int) * j 
	    - array(sa->Cartan, rank * r + rb, int) 
	    ;
	w21->k = locateWeight (sa, w21, TRUE)  ;
	w21->odd = TRUE ;
	w21->layer = 1 ;
      }
  if (s*b + a + 1 != 0)
    for (j = 0 ; j <= a - 1 ; j++)
      {
	w21 = arrayp (aa21, jj++, WW) ;
	/* position to the new weight */
	for (r = 0 ; r < rank ; r++)
	  w21->x[r] = w0->x[r] 
	    - array(sa->Cartan, rank * r + ra, int) * (j+1)
	    - array(sa->Cartan, rank * r + rb, int) 
	    ;
	w21->k = locateWeight (sa, w21, TRUE)  ;
	w21->odd = TRUE ;
	w21->layer = 1 ;
      }
  if (b*(s*b + a + 1) != 0)
    for (j = 0 ; j <= a ; j++)
      {
	w21 = arrayp (aa21, jj++, WW) ;
	/* position to the new weight */
	for (r = 0 ; r < rank ; r++)
	  w21->x[r] = w0->x[r] 
	    - array(sa->Cartan, rank * r + ra, int) * (j+1) 
	    - array(sa->Cartan, rank * r + rb, int) * 2
	    ;
	w21->k = locateWeight (sa, w21, TRUE)  ;
	w21->odd = FALSE ;
	w21->layer = 2 ;
      }
  
  return aa21 ;
} /* getSU21Crystal */

/******************************************************************************************************/

static BOOL demazureSU21 (SA *sa, Array wws, int rb, int *dimEvenp, int *dimOddp, BOOL show) 
{
  AC_HANDLE h = ac_new_handle () ;
  BOOL new = FALSE ;
  int ii, r, dimEven, dimOdd, ra = -1, rank = sa->rank ;
  
  for (r = 0 ; r < sa->rank ; r++)
    {
      if (sa->odd2[r])  /* distinguished even root */
	ra = r ;
    }
  if (ra == -1)
    messcrash ("demazureSU21, odd2 root not defined") ;
  
  if (sa->pass++ == 0) new = TRUE ;

  
  for (ii = 0 ; ii < arrayMax (wws) ; ii++) 
    {
      WW *w0 = arrayp (wws, ii, WW) ;
      Array aa21 ;
      int n1 = w0->mult ;
      int dn = 0 ;
      int n21 = w0->n21 ;
      int jMax ;

      if (ii == 0)
	{ w0->mult = 0 ; continue ; }
      if (n1 <= 0)
	continue ;
      if (n1 == n21)
	continue ;
      if (n1 < n21)
	messcrash ("mult < n21") ;

      dn = n1 - n21 ;  /* number of su21 multiplets to create */
        if (w0->x[ra] < 0)
	messcrash ("SU21 incomplete negative root w0.x[ra=%d]=%d", ra, w0->x[ra]) ;
      aa21 = getSU21Crystal (sa, w0, ra, rb, h) ;
      jMax = arrayMax (aa21) ;

        if (jMax > 0 && dn > 0)   /* populate the new multiplet */
	{
	    int j ;
	  WW *w21 = arrayp (aa21, jMax - 1, WW) ; /* lowest weight of the Crystal */	  
	  WW *w1 = arrayp (wws, w21->k, WW) ;	  
	  int dn = w0->mult - w0->n21 - w1->mult + w1->n21 ;
	  int dn21 = w0->mult - w0->n21 ;

	  if (jMax == 1)
	    { dn = 0 ; dn21 = w0->mult - w0->n21 ; }
	    
	  new = TRUE ;
	  for (j = 0 ; j < jMax ; j++)
	    {
	      w21 = arrayp (aa21, j, WW) ;	  
	      w1 = arrayp (wws, w21->k, WW) ;	  
	      w0 = arrayp (wws, ii, WW) ; /* may be relocates */
	      
	      if (! w1->k)
		{
		  w1->k = w21->k ;
		  for (r = 0 ; r < rank ; r++)
		    w1->x[r] = w21->x[r] ; /* coords */
		  w1->l2 = wwScalarProduct (sa, w1, w1) ;
		}
	      /* increase multiplicity of the new multiplet */
	      if (0 && j) w1->mult += dn ;
	      w1->n21 += dn21 ;
	      if (w1->mult < w1->n21) w1->mult = w1->n21 ;
	      w1->layer = w0->layer + w21->layer ;
	      if (w1->layer < 0)
		messcrash ("Negative layer") ;
	      w1->odd = w21->odd ^ w0->odd ; 
	      if (j == 0)
		w1->hw = w0->hw ;
	    }
	  printf ("++++ added %d su(2/1) states ", jMax) ;
	  wwsShow (sa, "+++ giving", -99, wws, &sa->hw) ;
	}
    }
  
  for (dimEven = dimOdd = 0, ii = 1 ; ii < arrayMax (wws) ; ii++)
    {
      WW *ww = arrp (wws, ii, WW) ;
        if (ww->odd)
	dimOdd += ww->mult ;
      else 
	dimEven += ww->mult ;
    }

  if (1 && new)
    {
      arraySort (wws, wwLayerOrder) ;

      printf ("................Demazure, pass %d, r=%d dimEven = %d dimOdd = %d\n", sa->pass, ra, dimEven, dimOdd) ;
      for (ii = 1 ; ii < arrayMax (wws) ; ii++)
	{
	  int j ;
	  WW *w = arrp (wws, ii, WW) ;
	  for (j = 0 ; j < rank ; j++)
	    printf (" %d", w->x[j]) ;
	  printf ("\tmult=%d n21=%d k=%d l=%d %s\n", w->mult, w->n21, w->k, w->layer, w->hw ? "*" : "" ) ;
	}
      printf ("\n") ;
  
      arraySort (wws, wwCreationOrder) ;
    }
  
  if (dimEvenp) *dimEvenp = dimEven ;
  if (dimOddp) *dimOddp = dimOdd ; 
  
  return new ;
} /* demazureSU21 */

/******************************************************************************************************/

static BOOL demazureOddDoublets (SA *sa, Array wws, int rb, int *dimEvenp, int *dimOddp, BOOL show) 
{
  BOOL new = FALSE ;
  int ii, r, dimEven, dimOdd, rank = sa->rank ;

  int rExtended = -1 ;
  for (r = 0 ; r < sa->rank ; r++)
    if (sa->extended [r])
      rExtended = r ;
  
  
  if (sa->pass++ == 0) new = TRUE ;

  for (ii = 0 ; ii < arrayMax (wws) ; ii++) 
    {
      WW *w0 = arrayp (wws, ii, WW) ;
      int dn, n1 = w0->mult ;
      int n21 = w0->n21 ;

      if (ii == 0)
	{ w0->mult = 0 ; continue ; }
      if (n1 <= 0)
	continue ;
      if (n1 == n21)
	continue ;
      if (w0->x[rExtended] <= 0)
	continue ;
      if (n1 < n21)
	messcrash ("mult < n21") ;

      dn = n1 - n21 ;  /* number of odd doublets to create */

      if (dn > 0)   /* populate the new multiplet */
	{
	  int j ;
	  WW *w1, ww2 ;

	  memset (&ww2, 0, sizeof (WW)) ;
	  new = TRUE ;
	  for (j = 1 ; j < 2 ; j++)
	    {
	      for (r = 0 ; r < rank ; r++)
		ww2.x[r] = w0->x[r] + sa->oddHw.x[r] ;
	      /* locate the k2 index of the WW weight structure of the corresponding member of the multiplet */
	      ww2.k = locateWeight (sa, &ww2, TRUE) ;
	      w1 = arrayp (wws, ww2.k, WW) ;	  
	      w0 = arrayp (wws, ii, WW) ; /* may have been be relocated */
	      if (! w1->k)
		{
		  *w1 = ww2 ;
		  w1->layer = w0->layer + 1 ;
		  w1->odd = ! w0->odd ; 
		}
	      /* increase multiplicity of the new multiplet */
	      w0->n21 += dn ;
	      w1->n21 += dn ;
	      if (w1->mult < w1->n21) w1->mult = w1->n21 ;
	    }
	  printf ("++++ added %d su(2/1) states ", dn) ;
	  wwsShow (sa, "+++ giving", -99, wws, &sa->hw) ;
	}
    }
  
  for (dimEven = dimOdd = 0, ii = 1 ; ii < arrayMax (wws) ; ii++)
    {
      WW *ww = arrp (wws, ii, WW) ;
        if (ww->odd)
	dimOdd += ww->mult ;
      else 
	dimEven += ww->mult ;
    }

  if (1 && new)
    {
      arraySort (wws, wwLayerOrder) ;

      printf ("................Demazure, pass %d, r=%d dimEven = %d dimOdd = %d\n", sa->pass, rb, dimEven, dimOdd) ;
      for (ii = 1 ; ii < arrayMax (wws) ; ii++)
	{
	  int j ;
	  WW *w = arrp (wws, ii, WW) ;
	  for (j = 0 ; j < rank ; j++)
	    printf (" %d", w->x[j]) ;
	  printf ("\tmult=%d n21=%d k=%d l=%d %s\n", w->mult, w->n21, w->k, w->layer, w->hw ? "*" : "" ) ;
	}
      printf ("\n") ;
  
      arraySort (wws, wwCreationOrder) ;
    }
  
  if (dimEvenp) *dimEvenp = dimEven ;
  if (dimOddp) *dimOddp = dimOdd ; 
  
  return new ;
} /* demazureOddDoublets */

/******************************************************************************************************/

static int demazure (SA *sa, int *dimEvenp, int *dimOddp, int method, BOOL show)
{
  BOOL ok = TRUE ;
  BOOL debug = FALSE ;
  int r, dimEven = 0, dimOdd = 0 ;
  
  int rExtended = -1 ;
  for (r = 0 ; r < sa->rank ; r++)
    if (sa->extended [r])
      rExtended = r ;
  
  if (method == 30)
    debug = TRUE ;
  if (debug)
    wwsShow (sa, "Before Demazure", 0, sa->wws, 0) ;

  while (ok)
    {
      BOOL ok2 = FALSE ;
      ok = FALSE ;
      for (r = 0 ; r < sa->rank ; r++)
	{
	  ok2 = FALSE ;
	  switch (method)
	    {
	    case 0:  /* Lie algebra, or even adjoint, apply all even roots */
	      if (! sa->odd[r] || sa->extended[r])
		ok2 = demazureEvenOdd (sa, sa->wws, r, TRUE, &dimEven, &dimOdd, show) ;
	      break ;
	    case 1:  /* negative odd roots */
	      if (! sa->odd[r] && ! sa->extended[r])
		ok2 = demazureEvenOdd (sa, sa->wws, r, TRUE, &dimEven, &dimOdd, show) ;
	      break ;



	    case 10:
	      
	    case 2:
	         if (! sa->odd[r])  /* even adjoint: extended root is accepted */
		   ok2 = demazureEvenOdd (sa, sa->wws, r, TRUE, &dimEven, &dimOdd, show) ;
		 break ;
	    case 12:
	    case 20:  /* generic representation: exclude the su(2/1) distinguished subalgebra  */
	      if (! sa->odd1[r] && ! sa->odd2[r] && ! sa->extended[r])
		ok2 = demazureEvenOdd (sa, sa->wws, r, TRUE, &dimEven, &dimOdd, show) ;
	      else
		if (sa->odd1[r]) /* construct the distinguished SU(2/1) multioplets */
		  ok2 = demazureSU21 (sa, sa->wws, r, &dimEven, &dimOdd, show) ;
	      break ;
	    case 30:  /* generic representation: exclude the su(2/1) distinguished subalgebra  */
	      if (! sa->extended[r] && ! sa->odd[r])
		ok2 = demazureEvenOdd (sa, sa->wws, r, TRUE, &dimEven, &dimOdd, show) ;
	      break ;
	    }
	  if (ok2 && debug)
	    {
	      printf ("... weyl %d\t", r) ;
	      wwsShow (sa, "Inside Demazure", 0, sa->wws, 0) ;
	    }
	  ok |= ok2 ;
	}
      if (sa->hasOdd)
	{
	  switch (sa->method)
	    {
	    case 12:
	    case 10:
	    case 30:
	      for (r = 0 ; r < sa->rank ; r++)
		{
		  ok2 = FALSE ;
		  if (sa->odd[r]) /* construct the odd doubblets */
		    ok2 = demazureOddDoublets (sa, sa->wws, r, &dimEven, &dimOdd, show) ;
		  if (ok2 && debug)
		    {
		      printf ("... weyl %d\t", r) ;
		      wwsShow (sa, "Inside Demazure odd doublets", 0, sa->wws, 0) ;
		    }
		  ok |= ok2 ;
		}
	      break ;
	    case 50: /* obsolete */
	      for (r = 0 ; r < sa->rank ; r++)
		{
		  ok2 = FALSE ;
		  if (sa->odd1[r]) /* construct the distinguished SU(2/1) multioplets */
		    ok2 = demazureSU21 (sa, sa->wws, r, &dimEven, &dimOdd, show) ;
		  if (ok2 && debug)
		    {
		      printf ("... weyl %d\t", r) ;
		      wwsShow (sa, "Inside Demazure su(2/1) multiplets", 0, sa->wws, 0) ;
		    }
		  ok |= ok2 ;
		}
	      break ;
	    }
	}
    }
  switch (method)
    {
    case 30:
      sa->method = 31 ;
      if (0)
	demazureEvenOdd (sa, sa->wws, rExtended, TRUE, &dimEven, &dimOdd, show) ;
      sa->method = 30 ;
      break ;
    }
  if (dimEvenp) *dimEvenp = dimEven ;  
  if (dimOddp) *dimOddp = dimOdd ; 

  if (show)
    wwsShow (sa, "Demazure", 0, sa->wws, 0) ;
  
  return dimEven ;
} /* demazure */

/******************************************************************************************************/

static void getNegativeOddRoots (SA *sa, BOOL show)
{
  int dimEven = 0 ;
  int dimOdd = 0 ;
  int method = sa->method ;

  sa->method = 0 ;			     
  sa->dict = dictHandleCreate (32, sa->h) ;
  
  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, -1, 1, 0) ;
  /* construct the first layer using Demazure */
  demazure (sa, &dimEven, &dimOdd, 1, show) ;
  sa->nOdd = dimOdd ; 
  sa->negativeOddRoots = sa->wws ;
  sa->wws = 0 ;
  sa->method = method ;
  
  if (show)
    wwsShow (sa, "Negative Odd roots", 1, sa->negativeOddRoots, arrp (sa->negativeOddRoots, 1, WW)) ;
  printf ("## Constructed %d odd roots\n", dimOdd) ;

  return ;
} /* getNegativeOddRoots */

/******************************************************************************************************/


static void getEvenAdjoint (SA *sa, BOOL show)
{
  int dimE = 0 ;
  int dimOdd = 0 ;
  int method = sa->method ;

  sa->method = 1 ;			     

  sa->dict = dictHandleCreate (32, sa->h) ;

  /* use as h.w. the lowering simple root */
  getHighestWeight (sa, -3, 1, 0) ;
  /* construct the even adjoint using Demazure */
  sa->nEven = demazure (sa, &dimE, &dimOdd, 0, show) ;

  sa->method = method ;
  sa->evenRoots = sa->wws ;
  sa->wws = 0 ;

  if (show)
    wwsShow (sa, "Even Adjoint ", 1, sa->evenRoots, 0) ;
  printf ("# Constructed %d adjoint even roots\n", dimE) ;

  return ;
} /* getEvenAdjoint */

/******************************************************************************************************/

static void getRho (SA *sa, BOOL show)
{
  WW *ww ;
  int ii, r ;
  int rank = sa->rank ;

  memset (&sa->rho0, 0, sizeof (WW)) ;
  /* rho0: sum of the positive even roots */
  if (sa->evenRoots)
    for (ii = 1 ; ii < arrayMax (sa->evenRoots) ; ii++)
      {
	int x, i, j ;
	BOOL isPositive = TRUE ;
	ww = arrp (sa->evenRoots, ii, WW) ;
	if (show)
	  {
	    printf (".... even root %d :: ", ii) ;
	    for (i = 0 ; i < rank ; i++)
	      printf ("%d ", ww->x[i]) ;
	    printf ("  :: ") ;
	  }
	for (i = 0 ; i < rank ; i++)
	  {
	    x = 0 ;
	    for (j = 0 ; j < rank ; j++)
	      x +=  arr (sa->upperMetric, rank * i + j, int) * arr (sa->scale, j, int) * ww->x[j] ;
	    if (x < 0)
	      isPositive = FALSE ;
	    if (show)
	      printf ("%d ", x) ;
	  }
	if (show)
	  printf ("\n") ;
	if (isPositive)
	  for (r = 0 ; r < rank ; r++)
	    sa->rho0.x[r] += ww->x[r] ; /* use + since we deal with the negative odd roots */
      }
  /* rho1: sum of the null positive odd roots */
  memset (&sa->rho1, 0, sizeof (WW)) ;
  if (sa->negativeOddRoots)
    for (ii = 0 ; ii < arrayMax (sa->negativeOddRoots) ; ii++)
      {
	ww = arrp (sa->negativeOddRoots, ii, WW) ;
	for (r = 0 ; r < rank ; r++)
	  sa->rho1.x[r] += ww->x[r] ; /* use + since we deal with the negative odd roots */
      }
  for (r = 0 ; r < rank ; r++)
    sa->rho.x[r] = sa->rho0.x[r]  + sa->rho1.x[r] ;
  
  if (show)
    {
      int i ;
      printf ("#================================= Sum of the even positive roots : 2 rho_0 = ") ;
      for (i = 0 ; i < rank ; i++)
	printf (" %d", sa->rho0.x[i]) ;
      
      if (sa->negativeOddRoots)
	{
	  printf ("\n#================================= Sum of the odd negative roots : 2 rho_1 = ") ;
	  for (i = 0 ; i < rank ; i++)
	    printf (" %d", sa->rho1.x[i]) ;
	}
      printf ("\n") ;
    }
} /* getRho */

/******************************************************************************************************/

static void getAtypic (SA *sa, BOOL show)
{
  WW hwT ;
  int rank = sa->rank ;
  int ii, r ;
  Array oddRoots = sa->negativeOddRoots ;
  
  sa->atypic = arrayHandleCreate (arrayMax (sa->negativeOddRoots), int, sa->h) ;  

  /* use as h.w. the declared h.w. */
  getHighestWeight (sa, -2, 1, TRUE) ;
  hwT = sa->hw ;  /* the highest weight L */ 
  for (r = 0 ; r < rank ; r++)
    hwT.x[r] = 2 * sa->hw.x[r] + sa->rho.x[r] ;     /* hw.w translated 2(L + rho) */
      
  for (ii = 1 ; ii < arrayMax (oddRoots) ; ii++)
    {
      /* x = < L + rho | beta_i > */
      int i, j, x = 0 ;
      WW *ww = arrp (oddRoots, ii, WW) ;
      
      for (i = 0 ; i < rank ; i++)
	for (j = 0 ; j < rank ; j++)
	  x += arr (sa->scale, i, int) * ww->x[i] * arr (sa->upperMetric, rank*i + j, int) *  arr (sa->scale, j, int) * hwT.x[j] ;
      if (x == 0)
	{
	  array (sa->atypic, ii, int) = 1 ;
	  sa->hasAtypic = TRUE ;
	}
      x = 0 ;
      if (ii == 1 && ! sa->isBlack)
	{
	  for (i = 0 ; i < rank ; i++)
	    for (j = 0 ; j < rank ; j++)
	      x += arr (sa->scale, i, int) * ww->x[i] * arr (sa->upperMetric, rank*i + j, int) *  arr (sa->scale, j, int) * sa->rho.x[j] ;
	  if (x != 0)
	    messcrash ("The trivial representaion is not atypic 1") ;
	}
    }
  
  if (1 || show)
    {
      int i ;
      
      if (sa->hasAtypic)
	{
	  printf ("\n#=================================   Atypical") ;
	  for (i = 0 ; i < arrayMax (sa->atypic) ; i++)
	    if (array (sa->atypic, i, int))
	      printf (" %d", i) ;
	  printf ("  ::  2(hw+rho)=") ;
	}
      else
	printf ("\n#=================================  Typical ::   2(hw+rho)=") ;
      for (i = 0 ; i < rank ; i++)
	printf ("%3d", hwT.x[i]) ;
      printf ("\n") ;
    }

  return ;
} /* getAtypic */

/******************************************************************************************************/
/*************************************************************************************/
/***************************** Public interface **************************************/
/*************************************************************************************/

static void usage (const char commandBuf [], int argc, const char **argv)
{
  int i ;

  fprintf (stderr,
	   "// Usage: irrep -m 2 -n 1 -show -type A -w -... \n"
	   "//      try: -h --help  \n"
	   "// Example:  irrep ...\n"
	   "// -o fileName : output file name, equivalent to redirecting stdout\n"
  
	   "// -gzo : the output file is gziped\n"
	   "// --method [1|2|12|20|30] : default 1 \n"
	   "//     method, used while debugging to test different algorithms\n"
	   "//     -m 10 : even/odd Demazure\n"     
	   "// -t --type A|B|C|D|E|F|G : (super)algebra type\n"
	   "// -m <int> -n <int> : ranks\n"
	   "//    -t A -m <int> :  su(m) Lie algebra\n"
	   "//    -t A -m <int> -n <int>:  su(m/n) Lie-Kac superalgebra\n"
	   "//    -t B -m <int> :  so(2m+1)) Lie algebra\n"
	   "//    -t B -m <int> -n <int> :  osp(2m+1/2n) Lie-Kac superalgebra\n"
	   "//    -t C -m <int> :  sp(2m) Lie algebra\n"
	   "//    -t C -m 1 -n <int>:  osp(2/2n) Lie-Kac superalgebra\n"
	   "//    -t D -m <int> :  so(2m) Lie algebra\n"
	   "//    -t D -m <int> -n <int>:  osp(2m/2n) Lie-Kac superalgebra\n"
	   "//    -t D -m 2 -n 1 -alpha <int>:  D(2/1,alpha) Lie-Kac superalgebra\n"
	   "//    -t E -m 6|7|8 :  E6,E7,E8 Lie algebra\n"
	   "//    -t F -m 4 :  F4 Lie algebra\n"
	   "//    -t F -n 4 :  F4 Lie-Kac superalgebra\n"
	   "//    -t G -m 2 :  G2, Lie algebra\n"
	   "//    -t G -n 3 :  G3 Lie-Kac superalgebra\n"
	   "// -w --weights 1:0:2...  integer Dynkin weights\n"
	   "//    The number of weights must match the rank of the algebra\n"
	   "// --show\n"
	   "// --table\n"
	   ) ;

  
  if (argc > 1)
    {
      fprintf (stderr,
	       "//\n// You said: %s\n", commandBuf) ;
      fprintf (stderr,
	       "// ########## ERROR: I do not understand the argument%s ", argc > 2 ? "s" : "") ;
      for (i = 1 ; i < argc ; i++)
	fprintf (stderr, "%s ", argv[i]) ;
      fprintf (stderr, "\n") ;
    }
  exit (1) ;
} /* usage */

/*************************************************************************************/
/******************************************************************************************************/


int main  (int argc, const char **argv)
{
   AC_HANDLE h = handleCreate () ;
   SA sa ;
   BOOL show = FALSE ;
   BOOL hasOdd ;
   int ix, dimEven = 0, dimOdd = 0 ;
   char *cp, commandBuf [3200] ;
  
   memset (commandBuf, 0, sizeof (commandBuf)) ;
   freeinit () ;


   memset (&sa, 0, sizeof (sa)) ;
   sa.h = h ;

   for ( ix = 0, cp = commandBuf ;  ix < argc && cp + strlen (argv[ix]) + 1 < commandBuf + 31000 ; cp += strlen (cp), ix++)
    sprintf(cp, "%s ", argv[ix]) ;

   if (argc < 2)
     usage (commandBuf, argc, argv) ;
   if (getCmdLineBool (&argc, argv, "-help") ||
       getCmdLineBool (&argc, argv, "-h")
       )
     usage (commandBuf, 1, argv) ;
   
   sa.method = 10 ;   
   sa.type = "toto" ;
   sa.DynkinWeights = "0" ;
   sa.m = sa.n = -99999 ;
   getCmdLineInt (&argc, argv, "-m", &sa.m) ; 
   getCmdLineInt (&argc, argv, "-n", &sa.n) ;
   sa.alpha = 1 ;
   getCmdLineInt (&argc, argv, "--alpha", &sa.alpha) ;
   show = getCmdLineBool (&argc, argv, "--show") ;
   sa.table = getCmdLineBool (&argc, argv, "--table") ;
   getCmdLineInt (&argc, argv, "--method", &(sa.method)) ;
   if (!getCmdLineOption (&argc, argv, "--type", &sa.type))
     getCmdLineOption (&argc, argv, "-t", &sa.type) ;
   if (! getCmdLineOption (&argc, argv, "--weights", &sa.DynkinWeights))
     getCmdLineOption (&argc, argv, "-w", &sa.DynkinWeights) ;
   
   if (argc != 1)
     {
       fprintf (stderr, "unknown argument, sorry\n") ;
       usage (commandBuf, argc, argv) ;
     }
   
   switch (sa.method)
     {
     case 1:
     case 2 :
     case 12:
     case 10:
     case 20:
     case 30:
       break ;
     default:
       messcrash ("argument --method should be 1,2, 10, 12, 20, 30 [default 1]") ;
       break ;
     }

   if (sa.m == -99999)
     messcrash ("Please specify the rank of the algebra  -m -n") ;
   if (sa.m != -99999 && sa.m < 0) messcrash ("rank argument -m m of type %s must be positive or null", sa.type) ;
   if (sa.n != -99999 && sa.n < 0) messcrash ("rank argument -n n of type %s must be positive or null", sa.type) ;
   
   sa.dict = dictHandleCreate (32, sa.h) ;
   
   getCartan (&sa, show) ;
   getMetric (&sa, show) ;
   
   hasOdd = sa.hasOdd ;
   sa.hasOdd = FALSE ;
   if (hasOdd) /* do this first then destroy the dict */
     getNegativeOddRoots (&sa, show) ;
   /* apply Demazure of the even group */
   getEvenAdjoint (&sa, show) ;
   getRho (&sa, show) ;
   sa.hasOdd = hasOdd ;
   
   if (sa.hasOdd)
     {
       if (sa.method < 10) /* atypic */
	 getAtypic (&sa, show) ;

       if (sa.method < 10) /* contruct the kacCrystal */
	 getKacCrystal (&sa, show) ; 
     }
   
   /* construct the h.w of the top layer even module */
   sa.dict = dictHandleCreate (32, sa.h) ;
   getHighestWeight (&sa, -2, TRUE, TRUE) ;
   printf ("*************************** %s m=%d n=%d %s \n "
	   , sa.type, sa.m, sa.n, sa.DynkinWeights) ;

   if (!sa.hasOdd)
     demazure (&sa, &dimEven, &dimOdd, 0, FALSE) ;
   else
     {
       switch (sa.method)
	 {
	 case 1:
	 case 2:
	   getHwCrystal (&sa, &dimEven, &dimOdd, show) ;
	   /* complete the hw Crystal to a full module */
	   demazure (&sa, &dimEven, &dimOdd, sa.method, FALSE) ;
	   break ;
	 case 12: /* use crystal and su21 */
	   getHwCrystal (&sa, &dimEven, &dimOdd, show) ;
	   /* complete the hw Crystal to a full module */
	   sa.method = 2 ;
	   demazure (&sa, &dimEven, &dimOdd, 2, FALSE) ;
	   if (arrayMax (sa.wws))
	   {
	     int i ;     /* transfer mult to crystal, then use it as a limit in the su21 constrution */
	     for (i = 1 ; i <= arrayMax (sa.wws) ; i++)
	       {
		 WW *w = arrp (sa.wws, i, WW) ;
		 w->crystal = w->mult ;
		 if (i > 1) w->mult = 0 ;
	       }
	   }
	   sa.method = 12 ;
	   demazure (&sa, &dimEven, &dimOdd, 12, TRUE) ;
	   break ;
	 case 10: 
	 case 20: 
	 case 30: /* use odd pairs */
	   demazure (&sa, &dimEven, &dimOdd, sa.method, TRUE) ;
	   break ;
	 }
     }
	   
   printf  ("################################################## Final Representation dimEven/dimOdd %d / %d\n",  dimEven, dimOdd) ;
   if (1) arraySort (sa.wws, wwLayerOrder) ;
   if (show) wwsShow (&sa, "Final representation", 1, sa.wws, &sa.hw) ;
   if (sa.hasOdd) getAtypic (&sa, show) ;
   messfree (h) ;
   printf ("A bientot\n") ;
   
   return 0 ;   
} /* main */

/************************************************************************/
/************************************************************************/
/************************************************************************/
/*
Variation:
  i again compute the shifted crystal and set w->mult = 0 is w->mult>x->mult
  this is favorable for F(4) adjoint because it kill the 1:0:0:0 dim 7 at level 0
  but this could be computed by hand, in the exact way
 
   BUG: the h.w. is not recovered by the shifted crystal which is absurd 
run -type D -m 2 -n 1 -alpha 3  -w 0:0:2 -show
*/

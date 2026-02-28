#ifndef FASTITOA__H
#define FASTITOA__H

/* this utility in tight loops is equivalent but enormously faster than
 *  sprintf (buf, "%u", v
 */

// 200-byte table → two digits at a time (00..99)
// Extremely fast for the common case (<100), still decent above
static const char digit_pairs[200] = {
    '0','0','0','1','0','2','0','3','0','4','0','5','0','6','0','7','0','8','0','9',
    '1','0','1','1','1','2','1','3','1','4','1','5','1','6','1','7','1','8','1','9',
    '2','0','2','1','2','2','2','3','2','4','2','5','2','6','2','7','2','8','2','9',
    '3','0','3','1','3','2','3','3','3','4','3','5','3','6','3','7','3','8','3','9',
    '4','0','4','1','4','2','4','3','4','4','4','5','4','6','4','7','4','8','4','9',
    '5','0','5','1','5','2','5','3','5','4','5','5','5','6','5','7','5','8','5','9',
    '6','0','6','1','6','2','6','3','6','4','6','5','6','6','6','7','6','8','6','9',
    '7','0','7','1','7','2','7','3','7','4','7','5','7','6','7','7','7','8','7','9',
    '8','0','8','1','8','2','8','3','8','4','8','5','8','6','8','7','8','8','8','9',
    '9','0','9','1','9','2','9','3','9','4','9','5','9','6','9','7','9','8','9','9'
};

// Fast path for most values (< 100), fallback for larger
// Buffer must have ≥12 bytes

static inline int fast_itoa (char *buf, unsigned int v)
{
  if (v < 100u)
    {
      // Most common case — 0..99
      const char *p = digit_pairs + v * 2 ;
      buf[0] = p[0] ;
      if (v < 10u) return 1 ; /* single digit */
      buf[1] = p[1] ;

      return 2 ;
    }

  /* larger numbers, double-digit steps */
  char *p = buf + 11;
  p-- ;
  
  while (v >= 100u)
    {
      unsigned idx = (v % 100u) << 1 ;
      p -= 2 ;
      *(unsigned short *)p = *(const unsigned short *)(digit_pairs + idx);
      v /= 100u ;
    }

  /*  Remainder */
  if (v > 0)
    {
      const char *pair = digit_pairs + v * 2 ;
      *--p = pair[1] ;
      if (v >= 10) *--p = pair[0] ;
    }
  
  int len = buf + 12 - p ;
  __builtin_memmove(buf, p, len) ;
  return len ;
} /* fast_itoa */

#endif

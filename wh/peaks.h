/*  File: peaks.h
 */

#ifndef DEF_PEAKS_H
#define DEF_PEAKS_H

void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels,
			Array aa ,    // array of unsigned int
			const char *target,
			int posMin, int step, // x = i * step + posMin
			int minCover,   // default 3 median (aa) no zero
			int scale // divide aa values by scale
			) ;

typedef struct cstStruct { int magic ; int dxMax, nu, bestShift ; double u1, v1, u2, v2, bestUv, *uu, *uv ; } *CST ;
CST cisTransCreate (int dxMax, AC_HANDLE h) ;
void cisTransCumulate (CST cst, int step, BigArray aaf, BigArray aar) ;
void cisTransNormalize (CST cst) ;

#endif /* DEF_PEAKS_H */

/**************************** End of File ******************************/

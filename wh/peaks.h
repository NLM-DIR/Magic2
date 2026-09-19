/*  File: peaks.h
 */

#ifndef DEF_PEAKS_H
#define DEF_PEAKS_H

void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels,
			Array aa ,    // array of unsigned int
			const char *target,
			int posMin, int step, // x = i * step + posMin
			int minCover   // default 3 median (aa) no zero
			) ;

#endif /* DEF_PEAKS_H */

/**************************** End of File ******************************/

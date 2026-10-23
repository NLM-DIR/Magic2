/*  File: peaks.h
 */

#ifndef DEF_PEAKS_H
#define DEF_PEAKS_H

void peaksCreateExport (ACEOUT ao, ACEOUT aoLevels,
			Array aa ,    // array of unsigned int
			const char *target,
			int posMin,   // x = i * step + posMin
			int step,     // step used in findPeaks
			int minCover,   // default 10 medianNoZero
			int scale       // divide aa values by scale
			) ;


#endif /* DEF_PEAKS_H */

/**************************** End of File ******************************/

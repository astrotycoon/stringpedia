/*
 * st_mwdc_header.h
 *
 *  Created on: 1 Sep 2010
 *      Author: ben
 */

#ifndef ST_MWDC_HEADER_H_
#define ST_MWDC_HEADER_H_

#include "sp_mwdc.h"

int getLineLength(FILE *file);

void runTest(char *text, char *pattern, char wildCardChar, int n, int m, char *startText, char *endText,int repeats,FILE *dataFile,short printResults,
		void (*functionPointer)(char *,char *,char,int,int,int *,struct SP_MWDC_MATCHING_POSITIONS*));

#define NAIVE 1
#define NAIVE_START "Naive Search"
#define NAIVE_FP &matchNaively

#define NAIVE_CONVOLUTIONS 2
#define NAIVE_CONVOLUTIONS_START "Naive Search Using Convolutions"
#define NAIVE_CONVOLUTIONS_FP &matchNaivelyWithConvolutions

#define FFT 3
#define FFT_START "FFT-based search (O(n log n) flavour) using FFTW"
#define FFT_FP &matchWithFFT

#define FFT_MEASURE 4
#define FFT_MEASURE_START "FFT-based search (O(n log n) flavour) using FFTW (FFTW_MEASURE)"
#define FFT_MEASURE_FP &matchWithFFT_Measure

#define FFT_PAD 5
#define FFT_PAD_START "FFT-based search (O(n log n) flavour) using FFTW (PADDING)"
#define FFT_PAD_FP &matchWithFFT_Pad

#define FFT_MEASURE_PAD 6
#define FFT_MEASURE_PAD_START "FFT-based search (O(n log n) flavour) using FFTW (FFTW_MEASURE,PADDING)"
#define FFT_MEASURE_PAD_FP &matchWithFFT_Measure_Pad

#define FFT_WISDOM 7
#define FFT_WISDOM_START "FFT-based search (O(n log n) flavour) using FFTW (SYSTEM WISDOM)"
#define FFT_WISDOM_FP &matchWithFFT_Wisdom

#define FFT_R2R 8
#define FFT_R2R_START "FFT-based search (O(n log n) flavour) using FFTW (R2R)"
#define FFT_R2R_FP &matchWithFFT_R2R

#define FFT_R2R_WISDOM 9
#define FFT_R2R_WISDOM_START "FFT-based search (O(n log n) flavour) using FFTW (R2R,WISDOM)"
#define FFT_R2R_WISDOM_FP &matchWithFFT_R2R_Wisdom

#define FFT_R2R_NWIT_WISDOM 10
#define FFT_R2R_NWIT_WISDOM_START "FFT-based search (O(n log n) flavour) using FFTW (R2R,WISDOM,NO WILDS IN TEXT)"
#define FFT_R2R_NWIT_WISDOM_FP &matchWithFFT_R2R_Wisdom_NWIT

#define FFT_NWIT 11
#define FFT_NWIT_START "FFT-based search (O(n log n) flavour) using FFTW (NO WILDS IN TEXT)"
#define FFT_NWIT_FP &matchWithFFT_NWIT

#define FFT_NWIT_WISDOM 12
#define FFT_NWIT_WISDOM_START "FFT-based search (O(n log n) flavour) using FFTW (NO WILDS IN TEXT,SYSTEM WISDOM)"
#define FFT_NWIT_WISDOM_FP &matchWithFFT_NWIT_Wisdom

#define FFT_FAST 13
#define FFT_FAST_START "FFT-based search (O(n log m) flavour) using FFTW"
#define FFT_FAST_FP &matchWithFasterFFT

#define FFT_FAST_MEASURE 14
#define FFT_FAST_MEASURE_START "FFT-based search (O(n log m) flavour) using FFTW (FFTW_MEASURE)"
#define FFT_FAST_MEASURE_FP &matchWithFasterFFT_Measure

#define FFT_FAST_PAD 15
#define FFT_FAST_PAD_START "FFT-based search (O(n log m) flavour) using FFTW (PADDING)"
#define FFT_FAST_PAD_FP &matchWithFasterFFT_Pad

#define FFT_FAST_MEASURE_PAD 16
#define FFT_FAST_MEASURE_PAD_START "FFT-based search (O(n log m) flavour) using FFTW (FFTW_MEASURE, PADDING)"
#define FFT_FAST_MEASURE_PAD_FP &matchWithFasterFFT_Measure_Pad

#define FFT_FAST_WISDOM 17
#define FFT_FAST_WISDOM_START "FFT-based search (O(n log m) flavour) using FFTW (SYSTEM WISDOM)"
#define FFT_FAST_WISDOM_FP &matchWithFasterFFT_Wisdom

#define FFT_FAST_MINSIZE 18
#define FFT_FAST_MINSIZE_START "FFT-based search (O(n log m) flavour) using FFTW (FORCED MIN SIZE)"
#define FFT_FAST_MINSIZE_FP &matchWithFasterFFT_MinSize

#define FFT_FAST_NWIT_WISDOM 19
#define FFT_FAST_NWIT_WISDOM_START "FFT-based search (O(n log m) flavour) using FFTW (NO WILDS IN TEXT, WISDOM)"
#define FFT_FAST_NWIT_WISDOM_FP &matchWithFasterFFT_NWIT_Wisdom

#define FFT_FAST_NWIT_MINSIZE 20
#define FFT_FAST_NWIT_MINSIZE_START "FFT-based search (O(n log m) flavour) using FFTW (NO WILDS IN TEXT, FORCED MIN SIZE)"
#define FFT_FAST_NWIT_MINSIZE_FP &matchWithFasterFFT_NWIT_MinSize

#define FFT_FAST_R2R 21
#define FFT_FAST_R2R_START "FFT-based search (O(n log m) flavour) using FFTW (R2R)"
#define FFT_FAST_R2R_FP &matchWithFasterFFT_R2R

#define FFT_FAST_R2R_WISDOM 22
#define FFT_FAST_R2R_WISDOM_START "FFT-based search (O(n log m) flavour) using FFTW (R2R,WISDOM)"
#define FFT_FAST_R2R_WISDOM_FP &matchWithFasterFFT_R2R_Wisdom

#define FFT_FAST_R2R_MINSIZE 23
#define FFT_FAST_R2R_MINSIZE_START "FFT-based search (O(n log m) flavour) using FFTW (R2R,FORCED MIN SIZE)"
#define FFT_FAST_R2R_MINSIZE_FP &matchWithFasterFFT_R2R_minsize

#define FFT_FAST_R2R_NWIT_MINSIZE 24
#define FFT_FAST_R2R_NWIT_MINSIZE_START "FFT-based search (O(n log m) flavour) using FFTW (R2R,FORCED MIN SIZE,NO WILDS IN TEXT)"
#define FFT_FAST_R2R_NWIT_MINSIZE_FP &matchWithFasterFFT_R2R_minsize_NWIT

#define FLINT 25
#define FLINT_START "Search using the FLINT library"
#define FLINT_FP &matchWithFlint

#define FFT_RANDOM 26
#define FFT_RANDOM_START "Randomized, FFT-based search (O(n log m))"
#define FFT_RANDOM_FP &matchWithRandomizedFFT

#define FFT_RANDOM_CHECK 27
#define FFT_RANDOM_CHECK_START "Randomized, FFT-based search (O(n log m)) (NAIVE VERIFICATION)"
#define FFT_RANDOM_CHECK_FP &matchWithRandomizedFFT_CheckNaively

#define FFT_RANDOM_R2R 29
#define FFT_RANDOM_R2R_START "Randomized, FFT-based search (O(n log m)),(R2R)"
#define FFT_RANDOM_R2R_FP &matchWithRandomizedFFT_R2R

#define FFT_RANDOM_R2R_CHECK 30
#define FFT_RANDOM_R2R_CHECK_START "Randomized, FFT-based search (O(n log m)) (R2R,NAIVE VERIFICATION)"
#define FFT_RANDOM_R2R_CHECK_FP &matchWithRandomizedFFT_R2R_CheckNaively

#define KALAI 31
#define KALAI_START "Kalai Algorithm (O(n log m))"
#define KALAI_FP &matchWithKalai

#define KALAI_CHECK 32
#define KALAI_CHECK_START "Kalai Algorithm (O(n log m)) (NAIVE VERIFICATION)"
#define KALAI_CHECK_FP &matchWithKalai_CheckNaively

#define KALAI_R2R 34
#define KALAI_R2R_START "Kalai Algorithm (O(n log m)) (R2R)"
#define KALAI_R2R_FP &matchWithKalai_R2R

#define KALAI_R2R_CHECK 35
#define KALAI_R2R_CHECK_START "Kalai Algorithm (O(n log m)) (R2R,NAIVE VERIFICATION)"
#define KALAI_R2R_CHECK_FP &matchWithKalai_R2R_CheckNaively

#endif /* ST_MWDC_HEADER_H_ */

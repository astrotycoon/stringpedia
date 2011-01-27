#!/bin/bash

if [ $# -ne 1 ]
then
	echo "It is recommended that you run this before using the rest of the stringpedia libraries; you will acheive better performance. You will also see the warning: 'Could not load wisdom!' until this has been run"
	echo "Usage ./wisdomgenerator timelimit. Timelimit is in hours, specify 0 to let this run until completion. Warning! Early termination causes FFTW to loose all accumulated wisdom!"
else
	fftw-wisdom -v -t $1 rof2 rob2 rof4 rob4 rof8 rob8 rof16 rob16 rof32 rob32 rof64 rob64 rof128 rob128 rof256 rob256 rof512 rob512 rof1024 rob1024 rof2048 rob2048 rof4096 rob4096 rof8192 rob8192 rof16384 rob16384 rof32768 rob32768 rof65536 rob65536 rob131072 rof131072 rof262144 rob262144 rof524288 rof524288 rob524288 rof1048576 rob1048576 rof2097152 rob2097152 ko2f ko2b ko4f ko4b ko8f ko8b ko16f ko16b ko32f ko32b ko64f ko64b ko128f ko128b ko256f ko256b ko512f ko512b ko1024f ko1024b ko2048f ko2048b ko4096f ko4096b ko8192f ko8192b ko16384f ko16384b ko32768f ko32768b ko65536f ko65536b ko131072b ko131072f ko262144f ko262144b ko524288f ko524288f ko524288b ko1048576f ko1048576b ko2097152f
fi

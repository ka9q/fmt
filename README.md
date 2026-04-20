**Frequency Measurement Test software**  
*Phil Karn, KA9Q*  
*19 April 2026*

This is software I threw together to analyze the IQ recordings I made during the April 17, 2026
ARRL Frequency Measuring Test. I must have done something right because I came in second place,
so I figure I should share it.

I used ka9q-radio (see https://github.com/ka9q/ka9q-radio) to make simultaneous IQ recordings of every frequency during the entire test period.
I used a 16 kHz sample rate, but any would do as long as the relevant frequencies are included. At the moment, the file must have two channels
(as implied by an IQ recording). It shouldn't be terribly hard to support mono files, ie, those from a SSB receiver, by simply using a real-to-complex
FFT instead of a complex-to-complex FFT.

The analysis involves running overlapping windowed FFTs across the analysis interval, estimating the peak frequency by quadratic interpolation. The FFTs are
zero padded (by a factor of 4 by default) to help the interpolation.

Then the statistical culling begins. First, windows with energy too far below the average are discarded. The default is -15 dB. Next, 
windows with estimated frequency too far from the median are discarded. The default is +/- 1 Hz.
The surviving windows are then sorted by frequency and the tails are trimmed off. The default is 10%.
Finally, a weighted mean of the trimmed windows is calculated along with an RMS error estimate.

Sample use:  
```
fmt -v -w 3.0 -i 1.5 -S '2026-04-17 02:49:06.587' -d 60.078 -l 90. -h 110 7065k2026-04-17T02:38:34.8Z.wav
```

command line arguments:  
```
-v turn on verbose output  
-w [seconds] FFT window duration  
-i [seconds] Amount to slide each FFT  
-S [date]    UTC of keydown time (requires unixstarttime attribute)  
  or  
-s [seconds] keydown time offset from file start  
-d [seconds] duration of keydown  
-l [Hz]      lower limit of search (may be negative)  
-h [Hz]      upper limit of search (may be negative)  
-t [0.0 to .99] Trim fraction (default 0.10)  
-O [oversample] FFT zero-padding ratio, default 4 (pad time domain data to 4x length)  
-o [Hz]      threshold for discarding frequency outliers from median (1 hz default)  
```

Two external file attributes are extracted from the .wav file:

*unixstarttime* (UTC start of file)  
*frequency* (radio frequency corresponding to 0 Hz in the IQ data)

If the *unixstarttime* attribute is present, you may specify the starting time of the measurement period with -S as a UTC time.
Otherwise use -s to give an offset in decimal seconds relative to the start of the file.

If the *frequency* attribute is present, the final estimate will be given as a radio frequency, otherwise it will be relative to zero frequency in the input file.

The -v option enables dumping of the individual FFT results: the time within the analysis
interval, the frequency estimate, the energy relative to the average, and several flags related to the statistical processing:

*outlier* means the estimated frequency was too far from the median estimate (default 1 Hz, change with -o)__
*weak* means the energy in the estimate was below a threshold relative to the average energy of all windows (default 15 dB, change with -m)  
*trimmed* means the (otherwise good) frequency estimate was in the bottom or top 10% (can be changed with -t)


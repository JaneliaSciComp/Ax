ACOUSTIC SEGMENTER (AX)

Given one or more time series of the same source, find tones that are
significantly above background noise.  Multi-taper spectral analysis is
used to identify time-frequency pixels containing signal and machine vision
techniques are used to cluster these into contours.  Output is bounding
boxes consisting of start and stop times, and low and high frequencies.
Designed for audio recordings of behaving animals.


SYSTEM REQUIREMENTS

A recent version of Matlab.  Tested with 2012a and 2013a.  Alternatively,
binary executables are available upon request.

Parallel Computing Toolbox and a computer with a lot cores is highly
recommended.


FILES

ax1.m extracts pixels containing signal from raw time series data.
This is the time limiting step.

ax2.m subsequently combines these pixels into contours.

cluster*.sh automates the analysis of multiple data sets by batching jobs.

compile.sh compiles ax1 and ax2 into binary executables for deployment.


See documentation embedded in scripts for more detail.

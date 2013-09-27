ACOUSTIC SEGMENTER (AX)

Given one or more time series of the same source, find tones that are
significantly above background noise.  Multi-taper harmonic analysis is
used to identify time-frequency pixels containing signal and machine vision
techniques are used to cluster these into contours.  Output are bounding
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


ALGORITHM

Overlapping segments in time are Fourier transformed (Guass 1805) using
multiple discrete prolate spheroidal sequences (Slepian 1961) as windowing
functions.  An F test (Fisher 1920) is used to infer whether the independent
estimates of intensity at each time-frequency point is significantly above
noise (Thomson 1982).  This procedure is performed for multiple segment
lengths on each microphone channel to capture data at different temporal
and spectral scales.  The data are combined in a single sonogram whose
pixel size corresponds to the time resolution of the shortest segment and
frequency resolution of the longest.  This image is then convolved with a
square box to fill in small gaps before the locations of contiguous pixels
exceeding a minimum area are characterized.


See documentation embedded in scripts for more detail.

ACOUSTIC SEGMENTER (AX)

Given one or more time series of the same source, find tones that are
significantly above background noise.  Multi-taper harmonic analysis is
used to identify time-frequency pixels containing signal and machine vision
techniques are used to cluster these into contours.  Output are bounding
boxes consisting of start and stop times, and low and high frequencies.
Designed for audio recordings of behaving animals.


SYSTEM REQUIREMENTS

A recent version of Matlab (at least 2012a), plus the image processing,
signal processing, and statistics toolboxes.  The groundtruthing script
also requires the mapping toolbox.

A computer with a lot cores is highly recommended, or better yet, a
small cluster.  The former requires the parallel computing toolbox,
and the latter additionally requires the distrubted computing
server or the matlab compiler.

Alternatively, ax1, the computational bottleneck, has been ported to both
Julia and Python.  Both are open source (i.e. free) and can access a cluster
without compilation or use of the unix command line.


FILES

ax1.m extracts pixels containing signal from raw time series data.
This is the time limiting step.

ax2.m subsequently combines these pixels into contours.

cluster*.sh automates the analysis of multiple data sets by batching jobs.

compile.sh compiles ax1 and ax2 into binary executables for deployment.


ALGORITHM

Overlapping segments in time are fast Fourier transformed (Guass 1805)
using multiple discrete prolate spheroidal sequences (Slepian 1961) as
windowing functions.  An F test (Fisher 1920) is used to infer whether
each time-frequency point is significantly above noise based on these
independent estimates of intensity (Thomson 1982).  This procedure is
performed for multiple segment lengths on each microphone channel to capture
data at different temporal and spectral scales.  The data are combined in a
single sonogram whose pixel size corresponds to the time resolution of the
shortest segment and frequency resolution of the longest.  This image is
then convolved with a square box to fill in small gaps before the locations
of contiguous pixels exceeding a minimum area are characterized.


See documentation embedded in scripts for more detail.

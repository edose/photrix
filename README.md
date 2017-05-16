# photrix
Photometric trickery for PLANNING and DATA REDUCTION supporting the **New Mexico Mira Project**'s programs in variable-star photometry, largely to support lightcurves administered by the [AAVSO](https://aavso.org%20%22), the American/Awesome/Allgemeine Association of Variable Star Observers.

The photrix code base has two functions:

**1. Planning:** this is in pretty good shape. After the SAS/AAVSO meeting mid-June 2017, I'll duplicate my poster in this repo. It tells all. Essentially, the code generates a roster of LPV and eclipser variable-star targets that:

 - are **available** tonight, for long enough above the horizon, and
   with no moon too close, *and*
 - have sufficient **priority**. This is a dynamic priority, which is lowered for a target if it has been recently observed according to AAVSO's WebObs database.

Then the user drags-and-drops star names from the daily roster to a planning spreadsheet and juggles the names at will. When the .txt summary sheet is sweetened to taste, the ACP (Astronomer's Control Panel) machine-readable plan files written at the same time can be uploaded directly to the telescope PC.

**2. Data reduction:** Still using R for this (repo 'Photometry-R'). I expect to have this ported to python, tested, and in daily use by July-August 2017. It will be kept in this repo, commingled more or less with the planning software.

## The code, such as it is
For 2 years the data-reduction functions were coded in and executed within R / RStudio, as published in repo 'Photometry-R'.

This newer repo is based on python 3.5 (as of May 2017). To be clear right up front, it does depend on some packages: 

 - ephem: for base astronomical calculations that I have no patience to reinvent (not compatible with all python versions)
 - pandas: makes python almost as good as R with no packages at all
 - pytest: skin-pricklingly-shiver-inducingly-awesome unit testing. Python's best feature.
 - requests: polls AAVSO's databases when needed
 - beautifulsoup4: parse raw data that 'requests' downloads
 - matplotlib: primitive graphics
 - statsmodels: mixed-model regression. Like R's lmer, but without all that bothersome clarity and accurate documentation

You might well ask: why port to python at all, given all my whinging & wailing. Yes, it's true, sometimes I throw manuals out the office window (impressive, as they are PDFs). You may read of py-ecstasy and py-agony in Python-Love-Hate.md (this repo).

## Status of this repository

**Planning sections** of this repository are functioning now. Not to say they couldn't be improved. They certainly will be.

**Data-reduction sections** of this repository are a work in progress, and frankly the R code is working so well already that my motivations are limited to future improvements that will require the object orientation and unit testing that python affords and that R doesn't. I expect to have a working data-reduction module process.py completed by July-August 2017, but no promises.

No other sections are expected. (The above two are plenty to deal with.)






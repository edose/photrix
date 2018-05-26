# photrix
Photometric trickery for SCHEDULING and DATA REDUCTION supporting the **New Mexico Mira Project**'s programs in variable-star photometry, largely to support lightcurves administered by the [AAVSO](https://aavso.org%20%22), the American/Awesome/Allgemeine Association of Variable Star Observers. The observing we're talking about here is: CCD-based, with an automated scope, mount, focuser, and camera, high up on a New Mexico mesa and operated from my small, warm, highly caffeinated office 60 miles away just off the Rio Grande in Albuquerque. 

A daily cycle is like this: 
 1. use the **photrix Scheduling modules** to make txt plans describing all the excrutiating details of the next night's observing,
 2. upload these txt plans to the scope PC during the day, 
 3. launch the first plan at dark, 
 4. get some sleep,
 5. download the resulting image files early the next morning,
 6. use the **photrix Data-Reduction modules** to turn the raw image files into shiny new lightcurve data points,
 7. upload the data points to AAVSO and review them in lightcurve context, and
 8. review the data to make adjustments to my target database for better future sessions.

> NOTE: The following is a May 2018 rewrite, made to support my presentation to the June 2018 SAS (Society for Astronomical Sciences) meeting in Ontario, California. This version describes the steps, in order, of my daily scheduling and data-reduction workflow, each step followed by a description of the software function that supports it.

# Workflows supported:
**1. Scheduling:** From my easily editable database of variable-star targets (one "FOV" or field-of-view description text file per main target), these photrix functions help me choose the best targets available for the upcoming observing night. The modules make a .csv (text) file which opens directly in Microsoft Excel, and I just drag and drop the target names I want (in order) to another Excel spreadsheet. Then more photrix modules make a summary file and ready-to-use ACP (Astronomer's Control Panel) plan files. At my fastest, I needed ca. 4 hours to make a night's plans by hand (badly). Now it takes me 30 minutes, and about 80% of my nights are free of errors or inefficiencies. Scroll a *little* ways down from here to read about **Workflow Part the First: SCHEDULING**

**2. Data reduction:** This new python code is ported from R code in old repo 'Photometry-R'). It was ported to python, tested, late July 2017 and has been in quite regular use since then. Scroll a *long* ways down from here to read about **Workflow Part the Second: DATA REDUCTION**


## This "photrix" code (and how I use it):
All code in photrix is written in python 3.5 or 3.6.  Yes, I used to write in R, but R's (1) lack of any testing facilities and (2) joke-class object orientation became two Giant Fails as my code grew. So python it is. You may read of my py-ecstasy and py-agony in Python-Love-Hate.md (this repo).

And this code does depend on some external python packages: 
 - **ephem**: for base astronomical calculations that I have no patience to reinvent (not compatible with all python versions)
 - **pandas**: this package makes python almost as good as R is with no packages at all
 - **pytest**: skin-pricklingly-shiver-inducingly-awesome unit testing. Python's best feature.
 - **requests**: polls AAVSO's databases when needed
 - **beautifulsoup4**: parse raw data that 'requests' downloads
 - **matplotlib**: primitive graphics
 - **statsmodels**: mixed-model regression. Like R's "lmer" package, but without all that bothersome clarity and accurate documentation

The status of this photrix code, both scheduling and data-reduction parts, are simple: working, fully unit tested, and in use for every observing night (> 200 nights/year). Nice.

FWIW, I run everything here within the Python Console window of the PyCharm Pro IDE. It's *way* easier than messing with python directly, it helps manage python versions and packages, and all my coding and testing time is within PyCharm Pro anyway, so I already know it well enough.

Some very useful photrix conventions that you might not have seen before: 
 - "an" or "AN" refers to an "Astronight", which is the local data in format "yyyymmdd" of the day when the nighttime observing session begins. So AN 20180102 is the observing session beginning January 2, 2018 local date, paying no attention to the Universal Time (UTC), and even if all your observations happen after midnight, that is, actually in January 3, 2018. This date convention is the same wise choice made in ACP's design, and removes practically all ambiguities in exactly which night one is referring to. The convention also removes any traditional clumsiness like "the night of January 2/3" etc.
 - "middark" time of an Astronight is the midpoint between the Astronight's sunset and sunrise times. (Surprisingly, it is not quite the same as time of the sun's lowest point below the horizon.) It has the virtue of being easy to calculate; also it avoids confusion with "midnight" as meaning the time 0000 or 12:00 AM and all the daylight savings smell that goes with that.
 - 


## FOV files: the text database of targets##
***Everything*** in photrix depends on these FOV files. So let's have a look now, shall we? Look at the FOV files used for testing, in this repo "*photrix/test/$FOVs_for_testing".*

Each FOV file contains everything "permanent" that photrix needs to know about one main variable star target, for both scheduling and data-reduction). And it's a lot, including:

 - the RA and Dec to point to when observing this main target;
 - the AAVSO chart number, needed to construct the uploaded final report;
 - the target type, e.g., Mira or Eclipser etc;
 - the variable star Period in days;
 - a Julian Date of brightness maximum and minimum (and for eclipsers the secondary minimum if known);
 - V magnitude of the brightness max and min (and for eclipsers that of the 2'min if known);
 - V-I color at the brightness max and min (and for eclipsers that of the 2'min if known);
 - Observing style, 
 - Target's priority, on scale of 0-10 (0=never observe, 4=typical, 8=pet target);
 - Gap score days, which defines how the urgency of observing depends on the number of days since last credibly observed (scale of 0-2, and multiplied by Target's priority to give overall priority); and
 - the STARS section which defines the sky position of targets including the main target, and which defines the sky position and standard magnitudes (and errors) of check and comp stars to be used in data reduction.

## Workflow Part the First: **SCHEDULING**
The python module planning ("pl") is used.

Do these steps (in this order) before each observing night. 

**Step 1.  Make a roster of available targets for the night.**

Start by running the function pl.roster (that is: the roster() function within the pl = planning.py module). A typical invocation would read something like:

    pl.make_an_roster(an_date_string='20180519', 
                      output_directory='C:/Astro/ACP/AN20180519', 
                      user_update_tolerance_days=0.1, 
                      exp_time_factor=0.78)

...where the parameters are:

 - *an_date_string* is the Astronight date string we're scheduling for;    
 - *output_directory* is where you want the roster .csv file to be written (Windows);
 - *user_update_tolerance_days* is the maximum age (in days) of the recent-observations cache that you will tolerate (0.1 days works very well); and
 - *exp_time_factor* is a general fudge factor applied to all this run's computed exposure times. The convention is that 1 is appropriate to give an estimated maximum star ADU of 16,000 at a FWHM of 4 pixels. This is usually a bit severe, and I use typically 0.6 (clear moonless night) to 1.0 (moonlit night with not really enough available targets to fill the night anyway)

...and the output is a .csv file to be opened in Microsoft Excel (see example file in this repo).

This *pl.make_an_roster()* function does the following:

 - makes a list of all FOV files within the FOVs directory;
 - for each FOV that is available this Astronight:
	 - call to AAVSO's WebObs web facility to download the latest 200 observations, and determine how long it's been since the last credible V-filter observations (at least 20 observations in one night for an eclipser),
	 - calculate the priority score for this night, remove from this Astronight's list if this FOV's priority is below minimum threshold;
	 - calculate its times (UTC) of availability, event times (e.g., eclipser minima), optimum exposure times for each filter specified in the FOV file, and total observation time (slewing + guiding startup + exposures + downloads, etc);
 - then compiles all the data into 4 tables, one line per star available this Astronight: standard stars, eclipser-like stars with high priority, eclipser-like stars with lower priority, and Mira/LPV stars; and
 - writes the 4 tables into the output .csv file. Tables are sorted by optimum time during the night, early to late.


**Step 2.  Choose target variable stars for this night.**

This is the fun part. You're going to drag target names from (1) the .csv file made in Step 1 above and opened in Excel to (2) a Microsoft Excel file "planning.xlsx" (from a previous night) as a template. Have them both open in separate instances of Excel.

First, update the AN date string at top of planning.xlsx and delete all the previous targets. Then pick the first high-priority target from the .csv file and drag the name from there into planning.xlsx. You arrange the target names in the order you want, and you can put in #PLAN cells to divide the cells into ACP plan files. Get planning.xlsx the way you think you want it, and then...


**Step 3.  Make the summary file and ACP plan files.**

This is the easiest part. You just run the function pl.plan (that is: the plan() function within the pl = planning.py module). A typical invocation would read something like:

    pl.make_an_plan(plan_excel_path='C:/Astro/ACP/AN20180519/planning.xlsx',
                    exp_time_factor=0.78)

...where the parameters are:

 - *plan_excel_path* is the full path to the planning.xlsx file; and
 - *exp_time_factor* is a general fudge factor applied to all this run's computed exposure times (may differ from the value in Step 1)

...and the output files are the following (see example file in this repo):

 - one summary file, e.g., *Summary_20180519.txt*; and
 - one or more ACP plan files, e.g., *plan_20180519_A.txt*.

Arrange the target-name cells in the order you want, top to bottom, then left to right, in the usual fashion. You'll need to include some standard targets during the night, for later data reduction--drag from the list of standards available that night at the roster file's top. You also have control over when ACP does autofocusing is done, via cells with AUTOFOCUS or AFINTERVAL in them.

The summary file is for your use (example in this repo). Each line of the report will contain a target name, estimated time (UTC) that the target's observing will start, and the degrees above the horizon during its observation. If the target is too close to the moon on this night (separate user-specified minimum), a warning message line is written just below the observing line. Similarly, if the target is very faint on this night (for Mira/LPV targets; user-specified max V-magnitude), a warning message line is written. 

The ACP plan files are for the observatory computer's use. Don't worry about these--they will absolutely be right (as of ACP 8.1) if the summary file is just the way you want it.

It won't be. You'll need to do Step 4.

**Step 4.  Iteratively improve the summary file.**

Basically, you'll do through a few cycles of tweaking the planning.xlsx and re-running *pl.make_an_plan()*, to get the night's sequence of targets just the way you want it. You'll keep in mind:

 - choosing high-priority targets, 
 - keeping the scope busy with them, 
 - avoiding warning messages in the summary file whenever possible, and
 - observing each target as high in the sky as possible.

This shuffling of targets within the planning.xlsx can get a bit tedious, but that is because it really is a tedious job. For most nights it takes 5-20 minutes to get it the way you want it. But when you get it there it really is nearly perfect, and you can do it night after night, even when dark is coming and you're rushed to get the scope started.

Compare this with optimizing by hand, taking 3-4 hours--and still such manual plans end up very far from perfect--bad exposure times, low targets, forgotten standards or autofocusing, typos, etc.

**Step 5.  Start the night's plans.**

Usually, by doing the following in about this order:

 - upload the plan_....txt file(s) to the observatory computer, 
 - print the summary file for your own use, 
 - prep the scope (camera cooled, focuser ok, scope unparked, ACP started and connected),
 - start the scope,
 - watch as it launches the first plan at dark (I use "WAITUNTIL -9" high in planning.xlsx).

If your summary looks good and the weather deities smile on your scope, it will run all night. You can watch a movie then sleep, and still the scope will be used more efficiently then you could ever arrange by hand. Most mornings, you will be rewarded with 100-300 FITS files waiting for data reduction. 

**A word about ACP Scheduler as an alternative.**

Many if not most variable-star observers use ACP Scheduler, a dispatch scheduler that tells ACP to do the best available target in its database, the one that is best at that moment. You can put on constraints including times, minimum sky elevation, etc, but you don't know in advance what it's going to decide to do--not in ordering, and not in absolute time--which is not a problem if your aim is observing the maximum number of high-priority targets in one night.

So why did I go to all the trouble to develop photrix's scheduling workflow and code?

Sigh--maybe I shouldn't have. Especially in light of the several months' full time development effort it ended up requiring. But develop I did, and here's why:

 - For Miras and other LPV targets, exposure times had to be computed in each filter, for each night. This is because these stars change brightness over an extremely wide range, and rather fast. If you    observe them weekly, successive observations may well differ by a whole magnitude.
 - For eclipsers and other time-series targets, I wanted real and detailed control, in particular to know exactly which eclipser would run each night, and between what start and stop times.

Probably I could have used ACP Scheduler as the scope-driving engine. Doing what I needed--to get the required nightly control over exposure times if nothing else--would have required that I write software to generate nightly RTML files for Scheduler's use. And that in turn would require that I develop the FOV system or something like it, and develop the exposure-time prediction system as well. The computation to generate nightly RTML files for ACP Scheduler already required most of what photrix scheduling does. So I decided to go independent. I'll probably never know whether that was the right choice.

And so I don't try to justify photrix as better than ACP Scheduler--for most observers, it probably isn't better. Most observers thinking about this journey should check out ACP Scheduler at DC-3 Dreams, or contact its hyper-talented developer, my friend Bob Denny. I can only justify that if one were going to build a comprehensive *advance* scheduler for variable star observing, one compatible with a robust and unit-tested data-reduction facility, that it would probably end up looking a lot like photrix.

## Workflow Part the Second: **DATA REDUCTION**
The python module process ("pr") is used.

Do these steps (in this order) after each observing night. There are a lot of steps, but fear not: the steps are the same for every observing night. We won't cover obtaining calibration frames (dark, bias, and flat images) for your rig, or how to get them ready to use in calibration each night's frames in MaxIm DL.

**Step 1.  Download all the night's images.**
And make a backup of the images if you'd like.

**Step 2.  Run pr.start().**
This python function arranges all the subdirectories (of the downloaded AN image directory) that you'll need for this workflow. It then moves the images into the "/Uncalibrated" subdirectory, makes a raw copy in the "/Ur" subdirectory. ("Ur" means "absolutely the original"; musicians refer to a composer's original manuscript, set for printing with the minimum number of changes, as "Urtext".) Other subdirectories made for later use are: "/Exclude" (where you'll drag and drop images to exclude them from processing), "/Photometry" (where your reduced data, including the final AAVSO-format upload file will go), and "FOV" (where photrix will store an exact, current-version copy of all the FOV files used in this data reduction session).

    pr.start(an_top_directory=AN_TOP_DIRECTORY, 
             an_rel_directory)

...where the parameters are:

 - *an_top_directory* is the directory holding all AN image data (default  = 'J:/Astro/Images/C14' which you will surely need to change)
 - *an_rel_directory* is the subdirectory under *an_top_directory*, typically like '20180519'

This also renames the files from ACP naming (e.g., "ST Tri-S001-R001-C001-I.fts" to simple and sequential names, by FOV name, in order of the images' date and time regardless of filter (e.g., "ST Tri-0003-I.fts").

**Step 3.  In MaxIm DL: Calibrate your images.**
This calibrates (removes flat and dark artifacts from) the images in "/Uncalibrated", putting the calibration images in "/Calibrated". It all happens in one operation.

First, make sure MaxIm's "Set Calibration" is ready to go with current calibration files, preferably masters. Then in MaxIm DL, under File > Batch Save & Convert, select all files in "/Uncalibrated" (make sure you're in the right night's directory), select destination as the "/Calibrated" subdirectory, make sure the "Perform Calibration" box is checked, and press "OK".

**Step 4.  In MaxIm DL: Look at your images. Exclude the rotten ones.**
Please don't skip this visual inspection step--it's easily set up in MaxIm DL and will save you many, many times its effort in later steps.

Make sure there are no files open in MaxIm. Then from Windows, select all the files in "/Calibrated" and drag them into MaxIm (may take a minute or two). The in MaxIm, click View > Animate, and in the small pop-up box select a time between images--I like 1 second--then press the forward button. Stop when you see a bad image and write down its filename. Keep going to the last image. When you're done, just drag the unwanted images into "/Exclude". They're safe but won't interfere with the later steps.

**Step 5.  Run pr.assess().**
This python function rigorously assesses the files in "/Calibrated" and writes to the screen a summary of the suspect ones. 

    pr.assess(an_top_directory=AN_TOP_DIRECTORY, 
             an_rel_directory)

...where the parameters are the same as for pr.start().

The tests run on each image filename found in "/Calibrated":
 - File name doesn't actually point to a subdirectory
 - Is readable as a FITS file (by astropy)
 - Is calibrated with flats and darks (currently must be calibrated by MaxIm 5 or 6)
 - Has plate solution info (either by PinPoint (MaxIm or ACP), or by TheSkyX's "Image Link")
 - FITS "OBJECT" header must match the file name
 - FWHM of stars must be reasonable
 - Computed focal length must be within 3% of mean FL of all FITS from this Astronight
 - File extension must be .fts (as written by pr.start() above)
 If all tests pass, you will see "ALL nnn FITS FILES APPEAR OK." where nnn is the number of valid FITS file found.

**Step 6.  Get the OK message from pr.assess() by whatever means necessary.**
Iteratively repair, plate-solve, or Exclude FITS files and re-run pr.assess() until you do get the FITS FILE APPEAR OK message.

The basis of all photrix is to pass only clean data to the next step.

**Step 7.  Run pr.make_df_master().**
This python function reads the FITS files and constructs the master Data Frame that will be used in all downstream calculations. The FITS files are not much used again.

    make_df_master(an_top_directory=AN_TOP_DIRECTORY, 
                   an_rel_directory,
                   instrument_name, 
                   ask_user)

...where the parameters are:

 - *an_top_directory* and *an_rel_directory* are the same as for pr.start()
 - *instrument_name* is the name of the instrument used to get the images. The name must refer to a pre-prepared instrument definition file placed in the "/instrument" directory. The file is short; its format is strict JSON. See for example file 'Borea.json' in this repository
 - *ask_user* is False if you want to skip the "Proceed?" prompt in the middle of *pr.make_df_master()*'s execution, otherwise just omit this parameter.

The dataframe is stored as a CSV-format file in the '/Photometry' subdirectory of this Astronight's directory for archival as well as retrieval for later steps. All FOV files used by this steps are written for archival to the '/FOV' subdirectory of this Astronight's directory.

**Step 8.  Run pr.SkyModel() for each filter.**
Here, we model the sky for this Astronight, that is, the zero-point and extinction (or you can use the default extinctions in the Site definition file, e.g., DSW.json in this repository), and include the effects of filter transforms (pre-measured and in the Instrument definition file). The function .SkyModel() is actually an object constructor; the SkyModel object gets a variable name, and it holds all the sky model information you'll need for all the images from this Astronight.

The fit is made only on comp stars (from target images) and standard stars (from standard-FOV images), and even then only on those stars with listed catalog errors of less than some threshold, typically 10 millimagnitudes.

The fit is a mixed-model regression, which allows normal linear regression on several parameters, and has the added advantage of rendering, as a by-product, any number of group parameters. The group parameter of interest here is the averaged deviation of all the comp stars in this

Typical example of running .SkyModel():
    
    V = pr.SkyModel(an_top_directory=AN_TOP_DIRECTORY, 
                    an_rel_directory=20180519, 
                    filter='V',
                    instrument_name='Borea', 
                    site_name='DSW',
                    fit_vignette=True,
                    fit_extinction=True)

...where the parameters are:

 - *an_top_directory* and *an_rel_directory* are the same as for pr.start()
 - *filter* is the name of the desired filter, as given in the FITS files
 - *instrument_name* is the name of the instrument used to get the images (as for *pr.make_df__master()*)
 - *site_name* is the site name (also as for *pr.make_df__master()*)
 - fit_vignette is True only if you want to fit a vignetting parameter (defined as a parabolic magnitude correction centered on the images; typically not needed when calibration flats are properly applied)
 - fit_extinction is True only if you want to fit the extinction to the Airmasses of the images in this filter. In my experience, this is only worthwhile if you have at least 150-200 observations in this filter with the required low catalog errors. Otherwise set this to False, and photrix will look up an extinction value from the Site definition file.

Additional parameters, rarely needed:

 - *max_cat_mag_error* (default=0.01 mag), with which you can set the maximum catalog error to be included in the sky model for this filter (leaving this at default is highly recommended)
 - *max_inst_mag_sigma* (default=0.03 mag), with which you can set the maximum sigma for a given star in a given image to allow its inclusion in the sky model for this filter (leaving this at default is highly recommended) saturation_adu defines the CCD peak ADU above which a star will be excluded from the model
 - *fit_sky_bias* (default=True) is True if a rounding-error bias term is to be included in the model. Very small effect.
 - *fit_xy* (default=False) is True if linear optical gradient correction terms are to be included in the model. Not affected by sky background gradients--is computed only from comp and standard star magnitudes.
 - *fit_transform* (default=False) is only set True when actually measuring Transform values (leaving this at default is nearly mandatory otherwise)
 - *fit_log_adu* (default=True) is True when a log(peak ADU) term is to be included in the model (generally leave this at default)
 - *do_plots* (default=True) is True when you want photrix to draw Q-Q and diagnostic plots to the computer screen

The **Q-Q plot** is a very standard statistical representation of how well a distribution matches a normal (Gaussian) distribution; in this case, the regression errors will be normal except for outliers. These outliers almost always indicate a comp-star or standard-star profile that was damaged by cosmic rays or other artefacts, and these points should be omitted (by including their data-point serial number in omit.txt, found in "/Photometry". Repeat the pr.SkyModel() fit and eliminate outliers until no outliers exist, and you have a consensus fit. Do this for each filter in which you took data.

The **diagnostic plots** are 12 plots arranged in one window. Each plot displays all the regression's input points in a way designed to expose a typical fault in data. For example, in one plot the data points' regression error is plotted vs its X position on the CCD. If the data appear to be higher on the left than on the right or vice versa, this would indicate that a X-gradient in the data was not properly calibrated away. This should be corrected by repeating all image calibrations and beginning again with pr.assess(), or simply by setting fit_xy to True in pr.SkyModel() and rerunning it.

In my experience, when the SkyModel --> Q-Q plot --> omit outliers --> rerun (until no outliers) cycle completes and no the Q-Q plot shows no obvious outliers, then all 12 diagnostic plots will also show no evidence of other regression problems. 

The assigned variable (e.g., "V" to left of the equal sign in the above example) holds the current SkyModel object holding all information about the regression on that filter's data. This is the input to the next step.

Between your refinement iterations of pr.SkyModel(), make your changes in the "/Photometry/omit.txt" text file. pr.SkyModel() reads omit.txt and prunes the input observation data set as its first step. 

Omit.txt starts as a stub file (help comment lines only), to which you add directives to tell pr.SkyModel() and later steps how you want the existing data set to be used. You may want observations:
 - to be omitted (because outliers, e.g.) by adding a #SERIAL directive to the text file "Photometry/omit.txt".

By adding directives to omit.txt, you may also:
 - omit one observation by star ID and image ID (#OBS directive);
 - omit a star from all images, in one filter or in all filters (#STAR directive);
 - omit all observations between two times (#JD directive); or
 - omit all observations from one image (#IMAGE directive).

An example of omit.txt:

    ;----- This is omit.txt for AN directory 20180525
    ;----- Use this file to omit observations from input to SkyModel (all filters).
    ;----- Example directive lines:
    ;
    ;#OBS    Obj-0000-V, 132 ; to omit star 132 from FITS image Obj-0000-V.fts
    ;#STAR   FOV, 132, V     ; to omit star 132 from all FITS with FOV and filter V
    ;#STAR   FOV, 132        ; to omit star 132 from all FITS with FOV and ALL filters
    ;#IMAGE  Obj-0000-V      ; to omit FITS image Obj-0000-V.fts specifically
    ;#JD     0.72, 1         ; to omit fractional JD from 0.72 through 1
    ;#SERIAL 123,77 54   6   ; to omit observations by Serial number (many per line OK)
    ;
    ;----- Add your directive lines:
    ;
    #SERIAL  1482         ; V outlier
    #IMAGE ST Tri-0001-V  ; had cloud interference

When you're happy with all the models in all filters, proceed to the next step.

**Step 9.  Run pr.PredictionSet().**
Here, we estimate magnitudes of all the unknown targets, which after all is the ultimate goal of this data-reduction workflow.

Typical example of running pr.PredictionSet():

    ps = pr.PredictionSet(an_top_directory=AN_TOP_DIRECTORY,
                          an_rel_directory='20180519', 
                          instrument_name='Borea', 
                          site_name='DSW',
                          skymodel_list=[V,R,I])

where the parameters are:
- *an_top_directory* and *an_rel_directory* are the same as for pr.start()
- *instrument_name* and *site_name* are the same as for pr.SkyModel()
- *skymodel_list* is simply a python list of the just-constructed SkyModel objects, the list built in place, via the syntax shown above

additional parameter, rarely used:
 - *max_inst_mag_sigma* (default=0.05): sets the reference instrumental magnitude error--any target data point having error above this reference will be silently omitted from the predictions (but will be   included in other output Data Frames to help in diagnosing errors). The default value is recommended.

Untransformed target (unknown) magnitudes are predicted from the target star flux as measured from the images, and from the SkyModel for each filter. The current version of photrix applies system transform values in V-I color, so for each reported target magnitude, both a V and I magnitude are required to compute the V-I color and thus the transform * color correction. This correction is straightforward for every target having only one V and one I magnitude as for Miras; for time series the color is spline-interpolated in time from V-I color values estimated from every adjacent pair of V and I raw points. At the end of this process, every target magnitude is transformed.

Before it exits, the PredictionSet constructor writes to the screen a list of targets (including check stars) that could not be predicted and which are omitted from reported results. There can be several causes for this, but the vast majority fail because for the V and/or I image for that target:

 - the V and/or I image is missing for that target, often because the user omitted the image for faults (in the pr.assess() loop)
 - the target was saturated, so that flux could not be estimated
 - target had an estimated instrumental magnitude error greater than the max_inst_mag_sigma threshold, often because the star was too faint (for the exposure time used)

There's currently no good system for determining which of the above causes applies to a given failure. In the current photrix version, the author investigates the failed main targets (which appear as the FOV name followed by the target star name, e.g. "ST Tri_ST_Tri") by opening the df_master Data Frame within PyCharm IDE and looking for the above 3 causes. This investigation is not required but is useful to help prevent failures for that target in future Astronights (for example, by adjusting exposure times, guide stars, etc). One may read and sort the stored DataFrame into a local variable by running in the Python Console, for example:

    df_master=pr.get_df_master(an_rel_directory='20180519')
                .sort_values(['ModelStarID', 'Filter'])

The assigned variable (ps in the PredictionSet() example above) is a PredictionSet object containing all information about the prediction set of all raw, transformed target magnitudes and their nominal prediction uncertainties (errors). This PredictionSet object is used in *all* remaining data-reduction and reporting steps.

If there are no time series observations in this Astronight data set, skip to Step 13!

**Step 10.  Run ps.stare_comps() for each filter.**
Here, we ensure that an absolutely uniform comp star set is applied to every reportable data point, within one time-series ("stare") target, within one filter. This is done to eliminate noise in the relative magnitudes caused merely by differences in the comp stars applied. That is, second-order color/transform effects aside, even substantial catalog errors in the comp-star magnitudes do not matter to the target star's *relative* magnitude estimates so long as the same comp stars are used throughout the series--so it pays richly to eliminate this effect as this step does.

Typical example of running ps.stare_comps():

    ps.stare_comps(fov="ST Tri",
                   star_id="ST Tri",
                   this_filter="V")

In the common case that a comp star profile was damaged in a single image (e.g., a by cosmic ray), the user must choose between eliminating the comp star (even better if there are several), or eliminating the time-series data point (again, little loss if the time-series was densely observed). Almost always, reducing the comp star set is better.

The output from ps.stare_comps() looks like this:
XXX

To set the comp star set for a given time-series in a given filter, simply copy the desired line from the *ps.stare_comps()* output, and paste it into *stare_comps.txt* within the "/Photometry" subdirectory (we'll later re-run PredictionSet() to use this change).

To eliminate a datapoint from this time-series in this filter, make a note of its image filename and omit the point later (by including it in a #SERIAL directive in markup_report.txt).

Repeat this process for each filter in which time-series data were taken (probably V and I).

**Step 11.  Re-run pr.PredictionSet().**
Skip this step if stare_comps.txt was not altered.

This step includes the beneficial effects of reducing comp star sets when needed.

**Step 12.  Run ps.stare_plot().**
Here, we persuade photrix to display a time series plot for one target, to detect faulty point not previously detected, as well as to enjoy a big visual payoff for all the hard work up to this point.

Typical example of running ps.stare_plot():

    ps.stare_plot(star_id="ST Tri")

**Step 13.  Run ps.markup_report().**
Here, we make a report, consisting of a table, one line to each target observation that is eligible to be reported to the AAVSO. The idea is to mark up this report, jotting down whether each line is: 
 - to be omitted (because it's an outlier, e.g.) by adding a #SERIAL directive to the text file "Photometry/report_map.txt";
 - to be combined with other lines by weighted averaging (when multiple exposures of the same field in the same filter have been taken) by adding a #COMBINE directive; or
 - to be used just as it appears (the default--if you do nothing, the line will be included as is in the AAVSO upload report).

By adding directives to report_map.txt, you may also:
 - omit all observations between two times (#JD directive); or
 - omit all observations for one target star (#TARGET directive).

An example of omit.txt 

    ;----- This is report_map.txt for AN folder 20180525
    ;----- Use this file to omit and/or combine target observations from AAVSO report.
    ;----- Example directive lines:
    ;
    ;#TARGET  GD Cyg ; to omit this target star altogether from AAVSO report.
    ;#JD  0.233 0.311 ; to omit this JD (fractional) range from AAVSO report.
    ;#SERIAL  34 44,129  32  1202 ; to omit these 5 Serial numbers from AAVSO report.
    ;#COMBINE  80,128 ; to combine (average) these 2 Serial numbers within AAVSO report.
    ;----- Add your directive lines:
    ;
    #TARGET  ST Tri       ;  try again tomorrow night
    #SERIAL  908 909 910  ;  bad comp stars
    #COMBINE  657 669     ; duplicate obs

Other than hand-editing the final AAVSO report (*not* recommended except for adding annotations), this is the user's last chance to affect what is uploaded directly into the AAVSO database. So the user should absolutely proofread this *report_map.txt* against the *markup_report.txt* to be sure that the next and final data-reduction computation step will do exactly what the user wants.
 
 **Step 13.  Run ps.aavso_report(), and upload the report to AAVSO.**
When the user has verified carefully that report_map.txt does in fact represent his intentions, simply run:

    ps.aavso_report()

The resulting file will be named e.g., AAVSOreport-20180525.txt. Proof this file to be sure it didn't get corrupted somehow (very unlikely), then upload it to the AAVSO upload page.

**Step 14 (recommended).  Review the just-uploaded data in context using the AAVSO Light Curve Generator (LCG).**
You won't often find a problem in the just-uploaded data, but it does happen. Especially cosmic-ray affected problems are often found only in context. And sometimes you'll find duplicate observations etc, and you'll want to remove them using AAVSO WebObs, but you can't remove them if you don't know about them.

I use the version 1 LCG found at https://www.aavso.org/lcg . I just can't stand the "Enhanced LCG" (never was a webpage less aptly named), version 2--its layout is an unholy mess, and I only hold my nose and use it when I must, e.g., to identify a data point by its observer. Version 1 gives better graphs in 10% of the time, and it's much, much clearer for our present review pruposes.

Yes, go through all the targets in your markup report. The most recent 500 days are fine for Miras and LPVs, the most recent 1 day is fine for time series. Don't forget to check all the filters you want to see. If you've promptly reduced and uploaded your data, yours will be the most recent data points, at the plot's far right. 

Congratulations! You're done, and the data-pruning, transforming, extinction, averaging, and proofing steps above have ensured that your data is among the best that's being uploaded to AAVSO.

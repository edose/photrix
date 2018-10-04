import os
from math import floor, sqrt, log10, log
from datetime import datetime, timezone
import shutil

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt

from .image import FITS, Image, R_DISC, Aperture
from .user import Instrument, Site
from .util import MixedModelFit, weighted_mean, jd_from_datetime_utc, RaDec
from .fov import Fov, FOV_DIRECTORY

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

THIS_SOFTWARE_VERSION = '2.0.0'  # as of 20170716
AN_TOP_DIRECTORY = 'J:/Astro/Images/C14'
DF_MASTER_FILENAME = 'df_master.csv'
FITS_REGEX_PATTERN = '^(.+)\.(f[A-Za-z]{2,3})$'
MIN_FWHM = 1.0
AAVSO_REPORT_DELIMITER = ','

START_PROCESSING_HERE___________ = ''

#######
#
#  photrix.process workflow:
#    import photrix.planning as pl, photrix.process as pr
#    pr.start(an_rel_directory='20170509')
#    [do MaxIm DL calibration]
#    pr.assess(an_rel_directory='20170509')
#    pr.make_df_master(an_rel_directory='20170509', instrument_name='Borea', ask_user=False)
#    V = pr.SkyModel(an_rel_directory='20170509', instrument_name='Borea', filter='V')
#    R = pr.SkyModel(an_rel_directory='20170509', instrument_name='Borea', filter='R')
#    I = pr.SkyModel(an_rel_directory='20170509', instrument_name='Borea', filter='I')
#       ... [whilte editing omit.txt until all V, R, and I models are right]
#    ps = pr.PredictionSet(an_rel_directory='20170816', skymodel_list=[V,R,I])
#    IF STARES exist in this AN:
#        ps.stare_comps(fov='ST Tri', star_id='ST Tri', this_filter='V')
#        ps.stare_plot(star_id='ST Tri')
#    ps.markup_report()
#    ps.aavso_report()
#    --> [upload to AAVSO]
#    --> [LCG review; edit FOV files if needed]
#    --> [check guiding exp times for possible FOV #CENTER (RaDec) adjustment]
#
#######
# Regex pattern to match Ur FITS filenames (from ACP):
#     r'^(.{3,}?)-S\d{3}-R\d{3}-C\d{3}-([a-zA-Z]{1,2}?)(_dupe-\d{1,4})?.f\w{1,2}'
# Regex pattern to match photrix FITS filenames
#     r'(.{3,}?)-\d{4}-([a-zA-Z]{1,2}?).f\w{1,2}'
#######


def start(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """First step in photrix processing pipeline.
    Always starts with raw directory exactly as downloaded from telescope PC.
    Typical usage: pr.start(an_rel_directory='20180804')
    :param an_top_directory:
    :param an_rel_directory:
    :return: None
    """
    # Copy original files to \Ur as backup (immutable):
    target_subdir = os.path.join(an_top_directory, an_rel_directory, 'Ur')
    if not os.path.exists(target_subdir) and not os.path.isdir(target_subdir):
        fits_subdir = os.path.join(an_top_directory, an_rel_directory)
        shutil.copytree(fits_subdir, target_subdir)
    else:
        print('>>>>> Could not create \'' + target_subdir + '\' as backup. (Already exists?)')

    # Move FITS to \Uncalibrated (not a full backup: some of these files may get Excluded):
    target_subdir = os.path.join(an_top_directory, an_rel_directory, 'Uncalibrated')
    if not os.path.exists(target_subdir) and not os.path.isdir(target_subdir):
        os.mkdir(target_subdir)
    else:
        print('>>>>> Could not create \'' + target_subdir + '\'. (Already exists?)')
    fits_subdir = os.path.join(an_top_directory, an_rel_directory)
    n_moved = 0
    for fits_entry in os.scandir(fits_subdir):
        if fits_entry.is_file():
            shutil.move(os.path.join(fits_entry.path),
                        os.path.join(target_subdir, fits_entry.name))
            n_moved += 1

    # Make remaining needed subdirectories:
    needed_subdirectories = ['Calibrated', 'Exclude', 'Photometry', 'FOV']
    for subdir_name in needed_subdirectories:
        new_subdir_path = os.path.join(an_top_directory, an_rel_directory, subdir_name)
        if not os.path.exists(new_subdir_path) and not os.path.isdir(new_subdir_path):
            os.mkdir(new_subdir_path)
        else:
            print('>>>>> Could not create \'' + new_subdir_path + '\'. (Already exists?)')

    # Arrange any calibration (incl. autoflat) files, write advisory if present:
    autoflat_path = os.path.join(an_top_directory, an_rel_directory, 'AutoFlat')
    if os.path.exists(autoflat_path) and os.path.isdir(autoflat_path):
        # Make calibration subdirectory if it doesn't already exist:
        calibration_path = os.path.join(an_top_directory, an_rel_directory, 'Calibration')
        if not os.path.exists(calibration_path) and not os.path.isdir(calibration_path):
            os.mkdir(calibration_path)
        else:
            print('>>>>> Could not create \'' + calibration_path + '\'. (Already exists?)')
        # Move all autoflat FITS files to /Calibration.
        if os.path.exists(calibration_path) and os.path.isdir(calibration_path):
            for autoflat_entry in os.scandir(autoflat_path):
                if autoflat_entry.is_file():
                    shutil.move(autoflat_entry.path, calibration_path)
        os.rmdir(autoflat_path)

    # Rename FITS files to photrix convention
    #     (within this rel_directory's /Uncalibrated subdirectory):
    _rename_to_photrix(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=an_rel_directory)

    print('.start() has moved', str(n_moved), 'FITS files to /Uncalibrated & has renamed them.')
    print('\n >>>>> Next:')
    print('    1. Calibrate with MaxIm now (File > Batch Save and Convert,')
    print('           from /Uncalibrated to /Calibrated.')
    print('    2. Visually inspect all FITS, e.g., with MaxIm')
    print('    3. Run assess().')


def assess(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Rigorously assess FITS files and directory structure for readiness to construct df_master.
    Collect and print all warnings and summary stats. Makes no changes to data.
    May be run as many times as needed, after start() and before make_df_master().
    Typical usage: pr.assess(an_rel_directory='20180811')
    :param an_top_directory: [string]
    :param an_rel_directory: [string]
    :return: [None]
    """
    # TODO: Add checks for guide exposure time. (?)
    # TODO: Add checks for binning=(1,1), when binning fields become available in FITS objects.
    # Make DataFrame of all files (& dirs) in directory, add some per-file info, and sort:
    filenames, isdir, extensions = [], [], []
    fits_path = os.path.join(an_top_directory, an_rel_directory, 'Calibrated')
    for entry in os.scandir(fits_path):
        filenames.append(entry.name)
        isdir.append(entry.is_dir())
        extensions.append(os.path.splitext(entry.name)[-1].lower())  # file extensions
    df = pd.DataFrame({'Filename': filenames, 'IsDir': isdir, 'Extensions': extensions},
                      index=filenames).sort_values(by=['Filename'])

    # Offer to delete any .src source files found (by-product of TheSkyX plate solutions).
    src_files = df.loc[df['Extensions'] == '.src', 'Filename']
    if len(src_files) >= 1:
        answer = input(' ..... ' + str(len(src_files)) +
                       ' .src files found. Delete them? (y/n, recommend y):')
        if answer.strip().lower()[0] == 'y':
            for filename in src_files:
                fullpath = os.path.join(an_top_directory, an_rel_directory, "Calibrated", filename)
                os.remove(fullpath)
        df = df.loc[~(df['Extensions'] == '.src'), :]  # remove df rows for files just deleted.

    # Directories: should be none; report and remove them from df:
    dirs = df.loc[df['IsDir'], 'Filename']
    if len(dirs) >= 1:
        print('Subdirectories found within /Calibrated (please remove them):')
        for this_dir in dirs:
            print('   ' + this_dir)
        df = df.loc[~df['IsDir'], :]  # remove rows referring to directories.
        del df['IsDir']  # as all rows in df refer to files and not directories.

    # Add empty columns to df:
    df['Valid'] = True
    df['PlateSolved'] = False
    df['Calibrated'] = False
    df['Object'] = ''
    df['ObjectMatchesName'] = False
    df['FovFileReady'] = False
    df['FWHM'] = np.nan
    df['FocalLength'] = np.nan

    # Try to open all filenames as FITS, collect info relevant to errors and warnings:
    fits_subdir = os.path.join(an_rel_directory, 'Calibrated')
    fov_name_cache = []
    for filename in df['Filename']:
        fits = FITS(an_top_directory, fits_subdir, filename)
        df.loc[filename, 'Valid'] = fits.is_valid
        if fits.is_valid:
            df.loc[filename, 'PlateSolved'] = fits.is_plate_solved
            df.loc[filename, 'Calibrated'] = fits.is_calibrated
            df.loc[filename, 'Object'] = fits.object
            df.loc[filename, 'ObjectMatchesName'] = filename.startswith(fits.object + '-')
            fov_proven_ready = False
            if fits.object in fov_name_cache:
                fov_proven_ready = True
            else:
                fov = Fov(fits.object)
                if fov.is_valid:
                    if fov.fov_name == fits.object:
                        fov_proven_ready = True
                        fov_name_cache.append(fits.object)
            df.loc[filename, 'FovFileReady'] = fov_proven_ready
            df.loc[filename, 'FWHM'] = fits.fwhm
            df.loc[filename, 'FocalLength'] = fits.focal_length

    # Non-FITS files: should be none; report and REMOVE THEM from df:
    invalid_fits = df.loc[~ df['Valid'], 'Filename']
    if len(invalid_fits) >= 1:
        print('\nINVALID FITS files:')
        for f in invalid_fits:
            print('    ' + f)
        print('\n')
        df = df.loc[df['Valid'], :]  # keep only rows for valid FITS files.
        del df['Valid']  # as all rows in df now refer to valid FITS files.
    else:
        print('All ' + str(len(df)) + ' files can be read as FITS files.')

    # Now assess all FITS, and report errors & warnings:
    not_calibrated = df.loc[~ df['Calibrated'], 'Filename']
    if len(not_calibrated) >= 1:
        print('\nNOT CALIBRATED:')
        for f in not_calibrated:
            print('    ' + f)
        print('\n')
    else:
        print('All calibrated.')

    not_platesolved = df.loc[~ df['PlateSolved'], 'Filename']
    if len(not_platesolved) >= 1:
        print('\nNO PLATE SOLUTION:')
        for f in not_platesolved:
            print('    ' + f)
        print('\n')
    else:
        print('All platesolved.')

    object_nonmatch = df.loc[~ df['ObjectMatchesName'], 'Filename']
    if len(object_nonmatch) >= 1:
        print('\nFITS Object does not match filename:')
        for f in object_nonmatch:
            # fits_object = FITS(an_top_directory, fits_subdir, f).object
            fits_object = df.loc[f, 'Object']
            print('    ' + f + ' has FITS Object = \'' + fits_object + '\'.')
        print('\n')
    else:
        print('All FITS objects match their filenames.')

    fov_file_not_ready = df.loc[~ df['FovFileReady'], 'Filename']
    if len(fov_file_not_ready) >= 1:
        print('\nFOV files ABSENT:')
        for f in fov_file_not_ready:
            fits_object = df.loc[f, 'Object']
            # fits_object = FITS(an_top_directory, fits_subdir, f).object
            print('    ' + f + ' is missing FOV file \'' + fits_object + '\'.')
        print('\n')
    else:
        print('All FOV files are ready.')

    odd_fwhm_list = []
    for f in df['Filename']:
        fwhm = df.loc[f, 'FWHM']
        if fwhm < 1.5 or fwhm > R_DISC:  # too small, or larger than half the aperture diameter:
            odd_fwhm_list.append((f, fwhm))
    if len(odd_fwhm_list) >= 1:
        print('\nUnusual FWHM (in pixels):')
        for f, fwhm in odd_fwhm_list:
            print('    ' + f + ' has unusual FWHM of ' + '{0:.2f}'.format(fwhm) + ' pixels.')
        print('\n')
    else:
        print('All FWHM values seem OK.')

    odd_fl_list = []
    mean_fl = df['FocalLength'].mean()
    for f in df['Filename']:
        fl = df.loc[f, 'FocalLength']
        if abs((fl - mean_fl)) / mean_fl > 0.03:
            odd_fl_list.append((f, fl))
    if len(odd_fl_list) >= 1:
        print('\nUnusual FocalLength (vs mean of ' + '{0:.1f}'.format(mean_fl) + ' mm:')
        for f, fl in odd_fl_list:
            print('    ' + f + ' has unusual Focal length of ' + str(fl))
        print('\n')
    else:
        print('All Focal Lengths seem OK.')

    # Set all FITS file extensions to '.fts' (MaxIm calibration sets it to '.fit' for some reason):
    _set_fits_extensions(an_top_directory=an_top_directory, fits_subdir=fits_subdir,
                         fits_filenames=df['Filename'])
    print('All FITS extensions are OK (=\'.fts\').')

    # Summarize and write instructions for next steps:
    n_warnings = len(not_calibrated) + len(not_platesolved) + len(object_nonmatch) +\
        len(fov_file_not_ready) + len(odd_fwhm_list) + len(odd_fl_list) + len(invalid_fits)
    if n_warnings == 0:
        print('\n >>>>> ALL ' + str(len(df)) + ' FITS FILES APPEAR OK.')
        print('Now...   1. Visually inspect all FITS files, if not already done.')
        print('         2. Run make_df_master().')
    else:
        print('\n >>>>> ' + str(n_warnings) + ' warnings (see listing above).')
        print('        Correct errors and rerun assess() until no errors remain.')


def make_df_master(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None,
                   instrument_name='Borea', ask_user=True):
    """Make the master DataFrame of all required information for downstream photometric processing.
    Typical usage: pr.make_df_master(an_rel_directory='20180811', ask_user=False)
    :param an_top_directory:
    :param an_rel_directory:
    :param instrument_name: name of Instrument (object) that took data [string]
    :param ask_user: True to ask user before building df_master; False to proceed directly.
    :return: [None] df_master is written as csv file to Photometry/df_master.txt.
    """
    # Build cross-reference DataFrame fits_fov_list:
    fits_fov_list = []
    fov_dict = {}
    an_directory = os.path.join(an_top_directory, an_rel_directory)
    fits_directory = os.path.join(an_directory, 'Calibrated')
    for entry in os.scandir(fits_directory):
        this_fits = FITS(an_top_directory, os.path.join(an_rel_directory, 'Calibrated'), entry.name)
        fov_name = this_fits.object
        fits_fov_list.append((entry.name, fov_name))  # list of 2-tuples
        if fov_name not in fov_dict.keys():
            fov = Fov(fov_name)
            fov_dict[fov_name] = fov  # caching Fov objects
    fits_list, fov_list = zip(*fits_fov_list)  # unzip 2-tuples to 2 parallel lists
    df_fits_fov = pd.DataFrame({'FITSname': fits_list, 'FOVname': fov_list}, index=fits_list)

    # Build display DataFrame one row per FOV file:
    df_ask_user = df_fits_fov.groupby(['FOVname']).count().sort_index().copy()
    df_ask_user.rename(columns={'FITSname': 'N_FITS'}, inplace=True)
    df_ask_user['FOV_file_exists'] = False  # default
    df_ask_user['FOV_file_exists_text'] = '*NA*'  # default
    df_ask_user['N_CheckStars'] = 0  # default
    df_ask_user['CheckMsg'] = ''  # default
    for fov_name in df_ask_user.index:
        this_fov = fov_dict[fov_name]  # these were cached above.
        df_ask_user.loc[fov_name, 'FOV_file_exists'] = this_fov.is_valid
        df_ask_user.loc[fov_name, 'N_CheckStars'] = sum([aa.star_type.lower() == 'check'
                                                         for aa in this_fov.aavso_stars])
    df_ask_user['FOV_file_exists_text'] = ['OK' if ex else 'MISSING'
                                           for ex in df_ask_user['FOV_file_exists']]
    df_ask_user['CheckMsg'] = ['OK' if n == 1 else 'WARNING: Target FOVs must have one Check Star.'
                               for n in df_ask_user['N_CheckStars']]

    # List FOVs and FITS file counts & FOV status, ask user to continue:
    # df_ask_user.drop(['FOV_file_exists', 'N_CheckStars'], axis=1, inplace=True)  # drop columns
    df_ask_user.rename(columns={'FOV_file_exists': 'FOV_exists', 'CheckMsg': 'CheckStar'},
                       inplace=True)
    df_ask_user['FOV'] = df_ask_user.index
    len_fov = max(len(s) for s in df_ask_user.index) + 1
    print('\n' + '  FOV'.ljust(len_fov), 'N FITS', 'FOV exists', 'Check Star')
    for ind in df_ask_user.index:
        print(df_ask_user.loc[ind, 'FOV'].ljust(len_fov),
              str(df_ask_user.loc[ind, 'N_FITS']).rjust(len('N FITS')),
              '    ' + df_ask_user.loc[ind, 'FOV_file_exists_text'].ljust(6),
              df_ask_user.loc[ind, 'CheckStar'])

    all_fovs_exist = all(df_ask_user['FOV_exists'])
    if not all_fovs_exist:
        print(' >>>>> STOPPING: at least one FOV file is missing.')
        return None
    print(' ------------------------------\n    ' + str(len(df_ask_user)) + ' FOVS.')

    if ask_user is True:
        answer = input('    .....Proceed? (y/n): ')
        if answer.strip().lower()[0] != 'y':
            print(' >>>>> STOPPING at user request.')
            return None
    else:
        print('ask_user = False, so make_df_master() continues...')

    print()
    instrument = Instrument(instrument_name)
    df_ur = pd.read_csv(os.path.join(an_top_directory, an_rel_directory,
                                     'Photometry', 'File-renaming.txt'),
                        sep=';', index_col='PhotrixName')
    fov_names = df_fits_fov['FOVname'].drop_duplicates().sort_values()
    df_master_list = []
    n_rows = 0
    for fov_name in fov_names:
        # Prepare data for this FOV:
        print('FOV >', fov_name)
        fov = fov_dict[fov_name]
        star_data = []
        for star in fov.aavso_stars:
                star_data.append((star.star_id, star.ra, star.dec, star.star_type, star.mags))
        star_id_list, ra_list, dec_list, star_type_list, mags_list = zip(*star_data)
        df_star_data_numbered = pd.DataFrame({'Number': range(len(star_data)),
                                              'StarID': star_id_list,
                                              'degRA': ra_list, 'degDec': dec_list,
                                              'StarType': [s.title() for s in star_type_list],
                                              'Mags': mags_list})
        df_star_data_numbered.index = range(len(df_star_data_numbered))
        if len(fov.punches) >= 1:
            punch_id_list, d_north_list, d_east_list = zip(*fov.punches)
            df_punches = pd.DataFrame({'StarID': punch_id_list,
                                       'dNorth': d_north_list, 'dEast': d_east_list})
        else:
            df_punches = pd.DataFrame()

        df_fov_list = []
        fits_names = df_fits_fov.loc[df_fits_fov['FOVname'] == fov_name, 'FITSname']

        for fits_name in fits_names:
            # Construct Image object for this FITS file:
            image = Image.from_fits_path(an_top_directory,
                                         os.path.join(an_rel_directory, 'Calibrated'), fits_name)
            for star_id, ra, dec in zip(df_star_data_numbered['StarID'],
                                        df_star_data_numbered['degRA'],
                                        df_star_data_numbered['degDec']):
                x0, y0 = image.fits.xy_from_radec(RaDec(ra, dec))
                image.add_aperture(star_id, x0, y0)
            image.add_punches(df_punches=df_punches)

            # Build df_apertures:
            ap_list = []
            ap_names = [k for k in image.apertures.keys()]
            for ap_name in ap_names:
                ap_list.append(dict(image.results_from_aperture(ap_name)))
            df_apertures = pd.DataFrame(ap_list, index=ap_names)  # constructor: list of dicts
            df_apertures['StarID'] = df_apertures.index
            df_apertures.rename(columns={'r_disc': 'DiscRadius',
                                         'r_inner': 'SkyRadiusInner',
                                         'r_outer': 'SkyRadiusOuter',
                                         'x_centroid': 'Xcentroid',
                                         'y_centroid': 'Ycentroid',
                                         'annulus_flux': 'SkyADU',
                                         'annulus_flux_sigma': 'SkySigma',
                                         'fwhm': 'FWHM',
                                         'x1024': 'X1024',
                                         'y1024': 'Y1024',
                                         'vignette': 'Vignette',
                                         'sky_bias': 'SkyBias'},
                                inplace=True)
            df_apertures = df_apertures.loc[df_apertures['net_flux'] > 0.0, :]
            df_apertures['InstMag'] = -2.5 * np.log10(df_apertures['net_flux']) + \
                2.5 * log10(image.fits.exposure)
            df_apertures['InstMagSigma'] = (2.5 / log(10)) * \
                                           (df_apertures['net_flux_sigma'] /
                                            df_apertures['net_flux'])  # math verified 20170726.
            df_apertures['ModelStarID'] = image.fits.object + '_' + df_apertures['StarID']
            df_apertures.drop(['n_disc_pixels', 'n_annulus_pixels', 'max_adu',
                               'net_flux', 'net_flux_sigma'],
                              axis=1, inplace=True)  # delete columns

            # For each aperture, add its max ADU from the original ("Ur", uncalibrated) FITS file:
            ur_filename = df_ur.loc[fits_name, 'UrName']
            df_apertures['UrFITSfile'] = ur_filename
            ur_image = Image.from_fits_path(an_top_directory, os.path.join(an_rel_directory, 'Ur'),
                                            ur_filename)
            df_apertures['MaxADU_Ur'] = np.nan
            df_apertures['LogADU'] = np.nan
            for star_id in df_apertures.index:
                ap = Aperture(ur_image, star_id,
                              df_apertures.loc[star_id, 'Xcentroid'],
                              df_apertures.loc[star_id, 'Ycentroid'],
                              df_punches=None)
                df_apertures.loc[star_id, 'MaxADU_Ur'] = ap.max_adu
                if ap.max_adu > 0.0:
                    df_apertures.loc[star_id, 'LogADU'] = log10(ap.max_adu)

            # Add FOV star data to each aperture row:
            df = pd.merge(df_apertures, df_star_data_numbered, how='left', on='StarID')
            df.sort_values(by='Number', inplace=True)  # stars in FOV order (not strictly needed).
            df.index = df['StarID']

            # Add catalog mag, CatMagError, and color index from FOV stars and Image's filter:
            df['CatMag'] = np.nan
            df['CatMagError'] = np.nan
            df['CI'] = np.nan  # old-school V-I color index
            for star_id in df['StarID']:
                mags = df.loc[star_id, 'Mags']
                # We extract from mags (dict) with .get() in case this filter is missing.
                mag_and_error = mags.get(image.fits.filter, (np.nan, np.nan))
                df.loc[star_id, 'CatMag'] = mag_and_error[0]
                df.loc[star_id, 'CatMagError'] = mag_and_error[1]
                # TODO: Make choice of color index passbands more flexible.
                ci_mag_1 = mags.get('V', (np.nan, np.nan))[0]  # old-school V-I color index
                ci_mag_2 = mags.get('I', (np.nan, np.nan))[0]  # old-school V-I color index
                df.loc[star_id, 'CI'] = ci_mag_1 - ci_mag_2    # old-school V-I color index
            df.drop('Mags', axis=1, inplace=True)

            # Add FITS data to all rows:
            df['FITSfile'] = image.fits.filename
            df['Object'] = image.fits.object
            df['JD_start'] = jd_from_datetime_utc(image.fits.utc_start)
            df['UTC_start'] = image.fits.utc_start
            df['Exposure'] = image.fits.exposure
            df['JD_mid'] = jd_from_datetime_utc(image.fits.utc_mid)
            df['Filter'] = image.fits.filter
            df['Airmass'] = image.fits.airmass

            # Add FOV data to all rows:
            df['FOV'] = fov.fov_name
            df['FOV_date'] = fov.fov_date
            df['Chart'] = fov.chart

            # Append this image's dataframe to list, write image line to console:
            df_fov_list.append(df)
            n_rows += len(df)
            if len(df) >= 1:
                print(fits_name, '=', len(df), 'rows',
                      '--> df_master_now_has', n_rows, 'rows.')
            else:
                print(fits_name, '=', len(df), 'rows',
                      '--> NO ROWS from this file ... df_master_now_has', n_rows, 'rows.')

        # Make df_fov by concatenating all and complete df_fov
        df_fov = pd.DataFrame(pd.concat(df_fov_list, ignore_index=True))
        df_fov = _add_ci_values(df_fov, df_star_data_numbered, instrument)
        df_master_list.append(df_fov)

    # Make df_master by concatenating all fov dataframes:
    df_master = pd.DataFrame(pd.concat(df_master_list, ignore_index=True))
    df_master.sort_values(['JD_mid', 'StarType', 'Number'], inplace=True)
    df_master.drop(['Number', 'Object'], axis=1, inplace=True)
    df_master.insert(0, 'Serial', range(1, 1 + len(df_master)))  # inserts in place
    df_master.index = list(df_master['Serial'])

    # For comp stars w/ CatMagError==NA, overwrite w/ largest CatMagError for same FOV and filter:
    # TODO: revisit using largest CatMagError for all CatMagError values. Seems too harsh.
    # TODO: move this loop to CatMagError computation (just above). Avoids .groupby().
    df_groupby = df_master[df_master['StarType'] == 'Comp'].groupby(['FOV', 'Filter'])
    for group_name, df_group in df_groupby:
        serials_nan_error = df_group.loc[np.isnan(df_group['CatMagError']), 'Serial'].tolist()
        serials_zero_error = df_group.loc[df_group['CatMagError'] <= 0.0, 'Serial'].tolist()
        serials_to_update = serials_nan_error + serials_zero_error
        if sum(serials_to_update) >= 1:
            max_catmagerror = df_group['CatMagError'].max()
            df_master.loc[serials_to_update, 'CatMagError'] = max_catmagerror

    # Write df_master to file:
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', DF_MASTER_FILENAME)
    df_master.to_csv(fullpath, sep=';', quotechar='"',
                     quoting=2, index=False)  # quoting=2-->quotes around non-numerics.

    # Finish & exit:
    _archive_fov_files(an_top_directory, an_rel_directory, fov_names)
    _write_omit_txt_stub(an_top_directory, an_rel_directory)
    print()
    # return df_master


class SkyModel:
    def __init__(self, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, filter=None,
                 instrument_name='Borea', site_name='DSW',
                 max_cat_mag_error=0.01, max_inst_mag_sigma=0.03, max_color_vi=+2.5,
                 saturation_adu=None,
                 fit_sky_bias=True, fit_vignette=True, fit_xy=False,
                 fit_transform=False, fit_extinction=True, fit_log_adu=True,
                 do_plots=True):
        """Constructs a sky model using mixed-model regression on df_master.
            Normally used by make_model()
        Typical usage: V = pr.SkyModel(an_rel_directory='20180804', filter='V', fit_extinction=True)
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param filter: name of filter to which this model applies [string, e.g., 'V' or 'R']
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_cat_mag_error: maximum catalog error allowed to stars in model [float]
        :param max_inst_mag_sigma: max instrument magnitude error allowed star observations [float]
        :param max_color_vi: maximum V-I color allowed to stars in model [float]
        :param saturation_adu: ccd ADUs that constitute saturation [float; from Instrument if None] 
        :param fit_sky_bias: True to fit sky bias term [bool]
        :param fit_log_adu: True to fit log(MaxADU_Ur) as CCD nonlinearity measure [bool]
        :param fit_vignette: True to fit vignette term (dist^2 from ccd center) [bool]
        :param fit_xy: True to fit X and Y gradient terms [bool]
        :param fit_transform: True to fit transform terms; else use values from Instument obj [bool]
        :param fit_extinction: True to fit extinction terms; else use values from Site obj [bool]
        Parameter 'fit_star_id' is not included in this version (has never been used in R, and
            would lead to crossed RE vars).
        """
        self.an_top_directory = an_top_directory
        self.an_rel_directory = an_rel_directory
        self.filter = filter
        self.instrument_name = instrument_name
        self.site_name = site_name
        self.max_cat_mag_error = max_cat_mag_error
        self.max_inst_mag_sigma = max_inst_mag_sigma
        self.max_color_vi = max_color_vi
        if saturation_adu is not None:
            self.saturation_adu = saturation_adu
        else:
            instrument = Instrument(self.instrument_name)
            self.saturation_adu = instrument.camera['saturation_adu']
        self.fit_sky_bias = fit_sky_bias
        self.fit_vignette = fit_vignette
        self.fit_xy = fit_xy
        self.fit_transform = fit_transform
        self.fit_extinction = fit_extinction
        self.fit_log_adu = fit_log_adu
        self.dep_var_name = 'InstMag_with_offsets'
        self.df_model = None    # data to/from regression, one row per input pt [pandas DataFrame]
        self.mm_fit = None      # placeholder [MixedModelFit object]
        self.df_star = None     # one row per unique model star [pandas DataFrame]
        self.extinction = None  # scalar result, placeholder
        self.transform = None   # "
        self.vignette = None    # "
        self.x = None           # "
        self.y = None           # "
        self.sky_bias = None    # "
        self.log_adu = None     # "
        self.converged = False  # "
        self.n_obs = None       # "
        self.n_images = None    # "
        self.sigma = None       # "
        self.df_image = None    # one row per image, placeholder

        # Rows from df_master, as curated by user in file omit.txt:
        df, warning_lines = _apply_omit_txt(self.an_top_directory, self.an_rel_directory)

        # Remove rows for several causes of ineligibility to help form a sky model:
        df = df[df['Filter'] == self.filter]
        df = df[df['StarType'] == 'Comp']
        df = df[df['CatMag'].notnull()]
        df = df[df['Airmass'].notnull()]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df = df[df['MaxADU_Ur'].notnull()]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['CI'].notnull()]
        df = df[df['CI'] <= self.max_color_vi]
        df = df[df['CatMagError'].notnull()]
        df = df[df['CatMagError'] <= self.max_cat_mag_error]
        df = df[df['FWHM'] >= MIN_FWHM]
        self.df_model = df

        self._prep_and_do_regression()
        self._build_output()
        if do_plots:
            self.plots()
        self.print_high_cirrus()
        # self.to_json_file()  # GIVE UP on JSON -- it can't handle DataFrames and Series.
        _write_stare_comps_txt_stub(self.an_top_directory, self.an_rel_directory)

    # @classmethod
    # def from_json(cls, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, filter=None):
    #     """
    #     Alternate constructor, reads from JSON file previously written.
    #        Normally used by _predict_fixed_only(), which requires final (immutable) Models.
    #     :param an_top_directory: path to an_rel_folder [str]
    #     :param an_rel_directory: folder for this instrument on this Astronight [string]
    #     :param filter: the filter to which this model applies [string, e.g., 'V' or 'R']
    #     :return: newly constructed Model object [class Model]
    #     """
    #     json_fullpath = os.path.join(an_top_directory, an_rel_directory, "model-" +
    #         filter + ".txt")
    #     with open(json_fullpath, 'r') as f:
    #         d = json.load(f)
    #     # Series and DataFrames (pandas) were stored as dictionaries, so convert them back:
    #
    #
    #
    #     # Distribute dictionary elements to their original object attributes:
    #     return cls(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory,
    #                filter=filter,
    #                instrument_name=d['instrument_name'],
    #                max_inst_mag_sigma=d['max_inst_mag_sigma'],
    #                max_cat_mag_error=d['max_cat_mag_error'],
    #                max_color_vi=+d['max_color_vi'], saturated_adu=d['saturation_adu'],
    #                fit_skyBias=d['fit_sky_bias'],
    #                fit_vignette=d['fit_vignette'], fit_xy=d['fit_xy'],
    #                fit_transform=d['fit_transform'], fit_extinction=d['fit_extinction'])
    #
    # def to_json_file(self):
    #     """
    #     Writes (most of) the current object to a JSON file.
    #        Does NOT include df_master (which is huge, but can be found in the same directory).
    #        Also does not include the statsmodels::MixedModelLm which cannot be serialized.
    #     :return: True if file successfully written, else False.
    #     """
    #     json_fullpath = os.path.join(self.an_top_directory, self.an_rel_directory, 'Photometry',
    #                                  "model-" + self.filter + ".json")
    #     json_dict = vars(self).copy()
    #
    #     # Convert pandas DataFrames to dictionaries without json-illegal int64, etc.
    #     json_dict['df_model'] = convert_pandas_to_json_compatible_dict(json_dict['df_model'])
    #     json_dict['df_star'] = convert_pandas_to_json_compatible_dict(json_dict['df_star'])
    #     # Convert pandas Series to dictionaries without json-illegal int64, etc.
    #     json_dict['mm_fit'].fitted_values = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].fitted_values)
    #     json_dict['mm_fit'].group_values = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].group_values)
    #     json_dict['mm_fit'].residuals = \
    #         convert_pandas_to_json_compatible_dict(json_dict['mm_fit'].residuals)
    #
    #     with open(json_fullpath, 'w') as f:
    #         json.dump(json_dict, f, indent=4)

    def _prep_and_do_regression(self):

        # Initiate dependent-variable offset, which will aggregate all such offset terms:
        dep_var_offset = self.df_model['CatMag'].copy()  # *copy* CatMag, or it will be damaged

        # Build fixed-effect (x) variable list and construct dep-var offset:
        fixed_effect_var_list = []
        if self.fit_transform:
            fixed_effect_var_list.append('CI')
        else:
            instrument = Instrument(self.instrument_name)
            transform_vi = instrument.transform(self.filter, 'V-I')
            dep_var_offset += transform_vi * self.df_model['CI']
            print(' >>>>> Transform (Color Index V-I) not fit: value fixed at',
                  '{0:.3f}'.format(transform_vi))
        if self.fit_extinction:
            fixed_effect_var_list.append('Airmass')
        else:
            site = Site(self.site_name)
            extinction = site.extinction[self.filter]
            dep_var_offset += extinction * self.df_model['Airmass']
            print(' >>>>> Extinction (Airmass) not fit: value fixed at',
                  '{0:.3f}'.format(extinction))
        if self.fit_sky_bias:
            if sum([x != 0 for x in self.df_model['SkyBias']]) > int(len(self.df_model) / 2):
                fixed_effect_var_list.append('SkyBias')
        if self.fit_log_adu:
            fixed_effect_var_list.append('LogADU')
        if self.fit_vignette:
            fixed_effect_var_list.append('Vignette')
        if self.fit_xy:
            fixed_effect_var_list.extend(['X1024', 'Y1024'])

        # Build 'random-effect' variable:
        random_effect_var_name = 'FITSfile'  # cirrus effect is per-image

        # Build dependent (y) variable:
        self.df_model[self.dep_var_name] = self.df_model['InstMag'] - dep_var_offset

        # Execute regression:
        self.mm_fit = MixedModelFit(data=self.df_model, dep_var=self.dep_var_name,
                                    fixed_vars=fixed_effect_var_list,
                                    group_var=random_effect_var_name)
        if self.mm_fit.statsmodels_object.scale != 0.0 and \
                self.mm_fit.statsmodels_object.nobs == len(self.df_model):
            print(self.mm_fit.statsmodels_object.summary())

    def _build_output(self):
        """
        Builds appropriate output attributes for external use. 
        :return: None
        """
        # Add 1/obs regression data as new df_model columns:
        self.df_model['FittedValue'] = self.mm_fit.df_observations['FittedValue']
        self.df_model['Residual'] = self.mm_fit.df_observations['Residual']
        dep_var_offset = self.df_model['InstMag'] - self.df_model[self.dep_var_name]
        self.df_model['FittedInstMag'] = self.df_model['FittedValue'] + dep_var_offset

        # Build df_star (star ID and count only):
        self.df_star = self.df_model[['Serial', 'ModelStarID']].groupby('ModelStarID').count()
        self.df_star['ModelStarID'] = self.df_star.index

        # Build df_image (from Mixed Model random effect), 1 row per FITS file,
        #    index = FITSfile, columns = FITSfile, JD_mid, Value:
        df = self.mm_fit.df_random_effects.copy()
        df = df.rename(columns={'GroupName': 'FITSfile', 'GroupValue': 'Value'})
        df_xref = self.df_model[['FITSfile', 'JD_mid']].drop_duplicates()
        df = pd.merge(df, df_xref, on='FITSfile', how='left', sort=False).sort_values(by='JD_mid')
        self.df_image = df.copy()
        self.df_image.index = self.df_image['FITSfile']

        # Extract and store scalar results:
        if self.fit_transform:
            self.transform = self.mm_fit.df_fixed_effects.loc['CI', 'Value']  # .loc(row, col)
        else:
            instrument = Instrument(self.instrument_name)
            self.transform = instrument.transform(self.filter, 'V-I')

        if self.fit_extinction:
            self.extinction = self.mm_fit.df_fixed_effects.loc['Airmass', 'Value']
        else:
            site = Site(self.site_name)
            self.extinction = site.extinction[self.filter]

        self.vignette = self.mm_fit.df_fixed_effects.loc['Vignette', 'Value'] \
            if self.fit_vignette is True else 0.0
        self.x = self.mm_fit.df_fixed_effects.loc['X1024', 'Value'] if self.fit_xy is True else 0.0
        self.y = self.mm_fit.df_fixed_effects['Y1024', 'Value'] if self.fit_xy is True else 0.0
        self.sky_bias = self.mm_fit.df_fixed_effects.loc['SkyBias', 'Value'] \
            if self.fit_sky_bias is True else 0.0
        self.log_adu = self.mm_fit.df_fixed_effects.loc['LogADU', 'Value'] \
            if self.fit_log_adu is True else 0.0
        self.converged = self.mm_fit.converged
        self.n_obs = len(self.df_model)
        self.n_images = len(self.df_model['FITSfile'].drop_duplicates())
        self.sigma = self.mm_fit.sigma
        print('\n', len(self.df_model), ' observations --> sigma=',
              round((1000.0 * self.sigma), 1), ' mMag')

    def _predict_fixed_only(self, df_predict_input):
        """
        Uses current model to predict best star magnitudes for *observed* InstMag and other inputs.
        FIXED-EFFECTS ONLY: does not include random effect.
        :param df_predict_input: data, all needed input columns for skymodel [pandas DataFrame]
            will NOT include or use random effects.
        :return: a dependent variable prediction for each input row [pandas Series of floats,
            with index = index of predict_input]
        """
        # First, verify that input is a DataFrame and that all needed columns are present.
        #    Names must be same as in model.
        if not isinstance(df_predict_input, pd.DataFrame):
            print('>>>>> SkyModel._predict_fixed_only(): predict_input is not a pandas DataFrame.')
            return None
        required_input_columns = ['Serial', 'FITSfile', 'InstMag', 'CI', 'Airmass']
        if self.fit_sky_bias:
            required_input_columns.append('SkyBias')
        if self.fit_log_adu:
            required_input_columns.append('LogADU')
        if self.fit_vignette:
            required_input_columns.append('Vignette')
        if self.fit_xy:
            required_input_columns.extend(['X1024', 'Y1024'])
        all_present = all([name in df_predict_input.columns for name in required_input_columns])
        if not all_present:
            print('>>>>> SkyModel._predict_fixed_only(): at least one column missing.')
            print('      Current model requires these columns:')
            print('         ' + ', '.join(required_input_columns))
            return None

        # Make parsimonious copy of dataframe; add bogus CatMag column (required by model):
        df_for_mm_predict = (df_predict_input.copy())[required_input_columns]
        bogus_cat_mag = 0.0
        df_for_mm_predict['CatMag'] = bogus_cat_mag  # totally bogus local value, reversed later

        # Execute MixedModelFit.predict(), giving Intercept + bogus CatMag + FEs + REs (pd.Series)
        #    DOES NOT INCLUDE RANDOM EFFECTS (these will be added back as Cirrus Effect terms):
        raw_predictions = self.mm_fit.predict(df_for_mm_predict, include_random_effect=False)

        # Compute dependent-variable offsets for unknown stars:
        dep_var_offsets = pd.Series(len(df_for_mm_predict) * [0.0], index=raw_predictions.index)
        if self.fit_transform is False:
            dep_var_offsets += self.transform * df_for_mm_predict['CI']
        if self.fit_extinction is False:
            dep_var_offsets += self.extinction * df_for_mm_predict['Airmass']

        # Extract best CatMag d'un seul coup, per (eq B - eq A), above:
        predicted_star_mags = \
            df_for_mm_predict['InstMag'] - dep_var_offsets - raw_predictions + bogus_cat_mag
        return predicted_star_mags

    def plots(self):
        # Setup for all Figures.
        obs_is_std = [name.startswith('Std_') for name in self.df_model['FOV']]
        obs_point_colors = ['darkgreen' if obs_is_std[i] is True else 'black'
                            for i, x in enumerate(obs_is_std)]
        image_is_std = [name.startswith('Std_') for name in self.df_image['FITSfile']]
        image_point_colors = ['darkgreen' if image_is_std[i] is True else 'black'
                              for i, x in enumerate(image_is_std)]
        jd_floor = floor(min(self.df_model['JD_mid']))
        obs_jd_fract = self.df_model['JD_mid']-jd_floor
        xlabel_jd = 'JD(mid)-' + str(jd_floor)
        obs_residuals_mmag = self.df_model['Residual'] * 1000.0

        # FIGURE 1 (Q-Q plot):
        from scipy.stats import norm
        df_y = self.df_model.copy()[['Residual', 'Serial']]
        df_y['Residual'] *= 1000.0
        df_y['Colors'] = obs_point_colors  # keep colors attached to correct residuals when sorted
        df_y = df_y.sort_values(by='Residual')
        n = len(df_y)
        t_values = [norm.ppf((k-0.5)/n) for k in range(1, n+1)]

        # Construct Q-Q plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 8))  # (width, height) in "inches"
        ax.grid(True, color='lightgray', zorder=-1000)
        ax.set_title('Q-Q plot of Residuals: ' +
                     self.an_rel_directory + '    ' + self.filter +
                     ' filter', color='darkblue', fontsize=20, weight='bold')
        ax.set_xlabel('t (sigma.residuals = ' + str(round(1000.0 * self.sigma, 1)) + ' mMag)')
        ax.set_ylabel('Residual (mMag)')
        ax.scatter(x=t_values, y=df_y['Residual'], alpha=0.6, color=df_y['Colors'], zorder=+1000)

        # Label potential outliers:
        mean_y = df_y['Residual'].mean()
        std_y = df_y['Residual'].std()
        z_score_y = (df_y['Residual'] - mean_y) / std_y
        df_y['T'] = t_values
        df_to_label = df_y[abs(z_score_y) >= 2.0]
        for x, y, label in zip(df_to_label['T'], df_to_label['Residual'], df_to_label['Serial']):
            ax.annotate(label, xy=(x, y), xytext=(4, -4),
                        textcoords='offset points', ha='left', va='top', rotation=-40)

        # Add reference line:
        x_low = 1.10 * min(df_y['T'])
        x_high = 1.10 * max(df_y['T'])
        y_low = x_low * std_y
        y_high = x_high * std_y
        ax.plot([x_low, x_high], [y_low, y_high], color='gray', zorder=-100, linewidth=1)

        # Add annotation: number of observations:
        fig.text(x=0.5, y=0.87,
                 s=str(len(self.df_model)) + ' observations in model.',
                 verticalalignment='top', horizontalalignment='center',
                 fontsize=12)
        fig.canvas.set_window_title(self.filter + ': Q-Q')
        plt.show()

        # FIGURE 2 (multiplot): Set up plot grid and style parameters:
        fig, axes = plt.subplots(ncols=4, nrows=3, figsize=(16, 10))  # (width, height) in "inches"

        def make_labels(ax, title, xlabel, ylabel, zero_line=True):
            ax.set_title(title, y=0.89)
            ax.set_xlabel(xlabel, labelpad=-27)
            ax.set_ylabel(ylabel, labelpad=-8)
            if zero_line is True:
                ax.axhline(y=0, color='lightgray', linewidth=1, zorder=-100)

        # Cirrus Plot (one point per image):
        ax = axes[0, 0]
        make_labels(ax, 'Image Cirrus Plot', xlabel_jd, 'mMag')
        ax.scatter(x=self.df_image['JD_mid']-jd_floor, y=self.df_image['Value'] * 1000.0,
                   alpha=0.6, color=image_point_colors)

        # Sky background vs JD_mid:
        ax = axes[0, 1]
        make_labels(ax, 'Sky background vs JD_mid', xlabel_jd, 'Sky ADU',
                    zero_line=False)
        ax.scatter(x=obs_jd_fract, y=self.df_model['SkyADU'],
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Instrument Magnitude:
        ax = axes[0, 2]
        make_labels(ax, 'Residuals vs Instrument Mag', 'Instrument Mag', 'mMag')
        ax.scatter(x=self.df_model['InstMag'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Max ADU:
        ax = axes[0, 3]
        xlabel_text = 'Max ADU, uncalibrated [log scale]'
        make_labels(ax, 'Residuals vs Max ADU', xlabel_text, 'mMag')
        ax.set_xlabel(xlabel_text, labelpad=-30)
        ax.set_xscale('log')
        x_scale_min = min(1000.0, 0.9 * min(self.df_model['MaxADU_Ur']))
        x_scale_max = 1.1 * max(self.df_model['MaxADU_Ur'])
        ax.set_xlim(x_scale_min, x_scale_max)
        # ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=20))
        ax.scatter(x=self.df_model['MaxADU_Ur'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Sky background:
        ax = axes[1, 0]
        make_labels(ax, 'Residuals vs Sky background', 'Sky ADUs', 'mMag')
        ax.scatter(x=self.df_model['SkyADU'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs JD:
        ax = axes[1, 1]
        make_labels(ax, 'Residuals vs JD', xlabel_jd, 'mMag')
        ax.scatter(x=obs_jd_fract, y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Color Index:
        ax = axes[1, 2]
        make_labels(ax, 'Residuals vs Color Index', 'Color Index (V-I)', 'mMag')
        ax.scatter(x=self.df_model['CI'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Airmass:
        ax = axes[1, 3]
        make_labels(ax, 'Residuals vs Airmass', 'Airmass', 'mMag')
        ax.scatter(x=self.df_model['Airmass'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Exposure Time:
        ax = axes[2, 0]
        make_labels(ax, 'Residuals vs Exposure Time', 'seconds', 'mMag')
        ax.scatter(x=self.df_model['Exposure'], y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs Vignette:
        ax = axes[2, 1]
        make_labels(ax, 'Residuals vs Vignette', 'pixels from CCD center', 'mMag')
        ax.scatter(x=1024*np.sqrt(self.df_model['Vignette']), y=obs_residuals_mmag,
                   alpha=0.6, color=obs_point_colors)

        # Residuals vs X:
        ax = axes[2, 2]
        make_labels(ax, 'Residuals vs X', 'X pixels from CCD center', 'mMag')
        ax.scatter(x=1024*self.df_model['X1024'], y=self.df_model['Residual'] * 1000.0, alpha=0.6,
                   color=obs_point_colors)

        # Residuals vs Y:
        ax = axes[2, 3]
        make_labels(ax, 'Residuals vs Y', 'Y pixels from CCD center', 'mMag')
        ax.scatter(x=1024*self.df_model['Y1024'], y=self.df_model['Residual'] * 1000.0, alpha=0.6,
                   color=obs_point_colors)

        # Finish the figure, and show the entire plot:
        fig.tight_layout(rect=(0, 0, 1, 0.925))
        fig.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.85, wspace=0.25, hspace=0.25)
        fig.suptitle(self.an_rel_directory +
                     '                   ' + self.filter + ' filter                ' +
                     '{:%Y-%m-%d     %H:%M  utc}'.format(datetime.now(timezone.utc)),
                     color='darkblue', fontsize=20, weight='bold')
        fig.canvas.set_window_title(self.filter + ': 12 plots')
        plt.show()

    def print_high_cirrus(self):
        df_cirrus = self.df_image.copy()
        df_cirrus['AbsValue'] = abs(df_cirrus['Value'])
        lines_to_print = max(6, int(0.05 * len(df_cirrus)))
        df_cirrus = df_cirrus.sort_values(by='AbsValue', ascending=False).head(lines_to_print)
        df_cirrus['mmag Cirrus'] = 1000.0 * df_cirrus['Value']
        print('--------------------------------------------')
        print(pd.DataFrame(df_cirrus['mmag Cirrus']))


class PredictionSet:
    def __init__(self, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None,
                 instrument_name='Borea', site_name='DSW',
                 max_inst_mag_sigma=0.05, skymodel_list=None):
        """        Constructs a prediction set, i.e., a set of best estimates of comp, check, and target star
        magnitudes, ready for marking up (curating) and reporting to the AAVSO (for example).
        Usages:
            ps = pr.PredictionSet(an_rel_directory='20170816', skymodel_list=[V,R,I])
            ps.markup_report()
            ps.aavso_report()
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_inst_mag_sigma: max instrument magnitude error allowed star observations [float]
        :param skymodel_list: SkyModel objects ready to be used [list of SkyModel objs]
        """
        self.an_top_directory = an_top_directory
        self.an_rel_directory = an_rel_directory
        self.instrument_name = instrument_name
        self.instrument = Instrument(self.instrument_name)
        self.site_name = site_name
        self.site = Site(self.site_name)
        self.max_inst_mag_sigma = max_inst_mag_sigma
        # ---- we'll deal with transform preferences in a later version.
        #      for now, we'll stick with V-I for all filters
        # if transform_preferences is None:
        #     transform_preferences = self.get_transform_preferences(self.instrument)
        # else:
        #     self.transform_preferences = transform_preferences
        self.saturation_adu = self.instrument.camera['saturation_adu']
        self.skymodels = dict((skymodel.filter, skymodel) for skymodel in skymodel_list)
        self.df_all_eligible_obs = None
        self.df_all_curated_obs = None
        # df_comps_mags: redefined 5/2017 to include all comps (not only from images w/targets)
        self.df_comp_mags = None
        self.df_cirrus_effect = None
        self.df_transformed = None
        # TODO: make the next two attributes either: both lists or both Serials.
        self.images_with_eligible_comps = None     # pd.Series of strings
        self.images_with_targets_and_comps = None  # list of strings

        # Do workflow steps:
        self.df_all_eligible_obs, warning_lines = _apply_omit_txt(self.an_top_directory,
                                                                  self.an_rel_directory)

        self.df_all_curated_obs, warning_lines = _curate_stare_comps(self.an_top_directory,
                                                                     self.an_rel_directory,
                                                                     self.df_all_eligible_obs)
        self.df_comp_mags = self.compute_comp_mags()

        self.df_cirrus_effect = self._compute_cirrus_effect(self.df_comp_mags)

        df_transformed_without_errors = self._compute_transformed_mags()

        print('\nWait for it...\n')

        self.df_transformed = self._compute_all_errors(df_transformed_without_errors)

        _write_aavso_report_map_stub(self.an_top_directory, self.an_rel_directory)

        self._write_summary_to_console()

    def compute_comp_mags(self):
        """
        Get raw, predicted mag estimates for ALL eligible comp-star observations, in all filters.
        This means for ALL images with or without targets, and with eligible comp observations.
        5/26/2017: Includes standards and other comp-containing images without targets.
        :return: dataframe of 
        """
        df = self.df_all_curated_obs.copy()
        # Select rows, then remove rows that suffer various causes of ineligibility:
        df = df[df['StarType'] == 'Comp']
        df = df[df['MaxADU_Ur'].notnull()]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df = df[df['CatMag'].notnull()]
        df = df[df['CI'].notnull()]
        df = df[df['Airmass'].notnull()]
        df = df[df['FWHM'] >= MIN_FWHM]
        df_comp_obs_with_estimates = df

        # Make lists of images with comps, and of images with comps&targets:
        self.images_with_eligible_comps = df_comp_obs_with_estimates['FITSfile'].drop_duplicates()
        images_with_targets = \
            (self.df_all_curated_obs[self.df_all_curated_obs['StarType'].str.lower() ==
                                     'target'])['FITSfile'].drop_duplicates()  # remove std FOVs etc
        self.images_with_targets_and_comps = [im for im in images_with_targets
                                              if im in self.images_with_eligible_comps.values]

        # Get best magnitude estimates for all eligible comps:
        df_list = []
        for filter_name, skymodel in self.skymodels.items():
            comp_obs_to_include = (df_comp_obs_with_estimates['Filter'] == filter_name) & \
                                  (df_comp_obs_with_estimates['FITSfile']
                                      .isin(self.images_with_eligible_comps))
            df_predict_input = df_comp_obs_with_estimates[comp_obs_to_include]
            df_estimates_this_skymodel = \
                df_predict_input[['Serial', 'ModelStarID', 'FITSfile', 'StarID', 'Chart',
                                  'Xcentroid', 'Ycentroid', 'InstMag', 'InstMagSigma', 'StarType',
                                  'CatMag', 'CatMagError', 'Exposure',
                                  'JD_mid', 'Filter', 'Airmass', 'CI', 'SkyBias', 'LogADU',
                                  'Vignette']]
            # column 'EstimatedMag' does NOT include image/cirrus effect!
            df_estimates_this_skymodel.loc[:, 'EstimatedMag'] = \
                skymodel._predict_fixed_only(df_predict_input)
            df_list.append(df_estimates_this_skymodel)  # list.append() performs in-place
        # Collect and return the dataframe:
        df_comp_mags = pd.concat(df_list, ignore_index=True)
        df_comp_mags['UseInEnsemble'] = True  # default, may be reset later
        df_comp_mags.index = df_comp_mags['Serial']
        return df_comp_mags

    def _compute_cirrus_effect(self, df_comp_mags):
        """
        Generates per-image cirrus effect with more careful inclusion/exclusion criteria than
            were used in the MixedModelFit, esp. in rejecting outlier comp observations.
        As of 5/26/2017: now includes images without targets, rendering resulting df
            a SUPERSET of R::df_cirrus_effect.
        :param: df_comp_mags: comprehensive comp-magnitude data [pandas DataFrame]. As of
            5/2017 includes all eligible comps from images with or without targets.
        :return: data for best per-image cirrus-effects, in magnitudes [pandas DataFrame] 
        """
        df_row_list = []
        for image in self.images_with_eligible_comps:
            df_estimates_comps_this_image = \
                df_comp_mags[df_comp_mags['FITSfile'] == image].copy()
            cirrus_effect_from_comps = \
                df_estimates_comps_this_image['EstimatedMag'] - \
                df_estimates_comps_this_image['CatMag']
            sigma2 = df_estimates_comps_this_image['CatMagError']**2 + \
                df_estimates_comps_this_image['InstMagSigma']**2
            least_allowed_sigma = 0.01
            # raw_weights ~ 1/s^2 where s cannot be less than least_allowed_sigma:
            raw_weights = pd.Series([1.0 / max(s2, least_allowed_sigma**2) for s2 in sigma2],
                                    index=sigma2.index)
            normalized_weights = raw_weights / sum(raw_weights)

            # Compute cirrus effect of mean value & sigma of mean (not sigma of indiv comp values):
            cirrus_effect_this_image, _, cirrus_sigma_this_image = \
                weighted_mean(cirrus_effect_from_comps, normalized_weights)
            if len(df_estimates_comps_this_image) == 1:
                cirrus_sigma_this_image = df_estimates_comps_this_image['CatMagError'].iloc[0]
            comp_ids_used = df_estimates_comps_this_image['StarID']  # starting point = use all
            num_comps_used = len(df_estimates_comps_this_image)      # "
            num_comps_removed = 0                                    # "

            # Reject this image's worst comp stars and recalculate, in the case
            #    (we start with at least 4 comp stars) AND (criterion1 >= 16 OR criterion2 >= 20).
            resid2 = (cirrus_effect_from_comps - cirrus_effect_this_image)**2  # pd.Series
            criterion1 = 0  # how many times worse is the worst comp vs avg of other comps
            criterion2 = 0  # square of worst comp's effective t-value, relative to CatMagError
            if len(df_estimates_comps_this_image) >= 4:
                x = normalized_weights * resid2  # pd.Series; will be >= 0
                c1 = x / ((sum(x) - x) / (len(x) - 1))  # pd.Series
                criterion1 = max(c1)
                c2 = raw_weights * resid2  # a pd.Series; will be >= 0
                criterion2 = max(c2)
                if (criterion1 >= 16) or (criterion2 >= 20):
                    # score > 1 for each row (comp) that may be removed:
                    score = pd.Series([max(c_1/16.0, c_2/20.0) for (c_1, c_2) in zip(c1, c2)],
                                      index=c1.index)
                    max_to_remove = floor(len(score) / 4)
                    to_remove = (score >= 1) & \
                                (score.rank(ascending=True, method='first') >
                                 (len(score) - max_to_remove))

                    # Now, set weights to zero for the worst comps, recalc cirrus effect & sigmas:
                    raw_weights[to_remove] = 0.0
                    normalized_weights = raw_weights / sum(raw_weights)
                    cirrus_effect_this_image, _, cirrus_sigma_this_image = \
                        weighted_mean(cirrus_effect_from_comps, normalized_weights)
                    n_nonzero_weights = sum([w != 0 for w in normalized_weights])
                    if n_nonzero_weights == 1:
                        cirrus_sigma_this_image = \
                            sum([sigma for (nwt, sigma)
                                 in zip(normalized_weights,
                                        df_estimates_comps_this_image['CatMagError'])
                                 if nwt != 0])

                    comp_ids_used = (df_estimates_comps_this_image['StarID'])[~ to_remove]
                    num_comps_used = len(comp_ids_used)
                    num_comps_removed = len(df_estimates_comps_this_image) - num_comps_used

                    # Record these omitted comp stars in the master comp-star dataframe.
                    removed_serials = (df_estimates_comps_this_image['Serial'])[to_remove]
                    df_comp_mags.loc[(df_comp_mags['Serial'].isin(removed_serials)),
                                     'UseInEnsemble'] = False

            # Insert results into this image's row in df_cirrus_effect:
            df_row_this_image = {'Image': image,
                                 'CirrusEffect': cirrus_effect_this_image,
                                 'CirrusSigma': cirrus_sigma_this_image,
                                 'Criterion1': criterion1,
                                 'Criterion2': criterion2,
                                 'NumCompsUsed': num_comps_used,
                                 'CompIDsUsed': ','.join(comp_ids_used),
                                 'NumCompsRemoved': num_comps_removed}
            df_row_list.append(df_row_this_image)
        df_cirrus_effect = pd.DataFrame(df_row_list)
        df_cirrus_effect.index = df_cirrus_effect['Image']
        return df_cirrus_effect

    def _compute_transformed_mags(self):
        """
        Return best mag estimates for all target and check stars, all observations, all filters.
        :return: 
        """
        # Construct df_input_checks_targets:
        df = self.df_all_eligible_obs.copy()
        df = df[[(st in ['Check', 'Target']) for st in df['StarType']]]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df['CI'] = 0.0  # because they are presumed unknown; populated later
        df['CatMagSaved'] = df['CatMag'].copy()
        df['CatMag'] = 0.0  # estimated below by imputation from predicted magnitudes
        df_input_checks_targets = df

        # Make df_estimates_checks_targets (which will account for Airmass(extinction),
        #    but not (yet) for Color Index (transforms) which will be handled below:
        df_filter_list = []
        for this_filter, skymodel in self.skymodels.items():
            rows_to_select = (df_input_checks_targets['Filter'] == this_filter) & \
                             (df_input_checks_targets['FITSfile']\
                              .isin(self.images_with_targets_and_comps))
            df_input_this_skymodel = (df_input_checks_targets.copy())[rows_to_select]
            predict_output = skymodel._predict_fixed_only(df_input_this_skymodel)
            df_estimates_this_filter = df_input_this_skymodel.copy()
            df_estimates_this_filter['PredictedMag'] = predict_output
            df_filter_list.append(df_estimates_this_filter)
        df_estimates_checks_targets = pd.concat(df_filter_list, ignore_index=True)
        df_estimates_checks_targets.index = df_estimates_checks_targets['Serial']  # to ensure
        columns_post_predict = ["Serial", "ModelStarID", "FITSfile", "StarID", "Chart",
                                "Xcentroid", "Ycentroid", "InstMag", "InstMagSigma", "StarType",
                                "CatMag", "CatMagSaved", "CatMagError", "Exposure", "JD_mid",
                                "Filter", "Airmass", "CI", "SkyBias", "Vignette", "LogADU",
                                "PredictedMag"]
        df_estimates_checks_targets = df_estimates_checks_targets[columns_post_predict]
        df_estimates_checks_targets['CatMag'] = df_estimates_checks_targets['CatMagSaved']
        df_estimates_checks_targets = df_estimates_checks_targets.drop(['CatMagSaved'], axis=1)
        df_estimates_checks_targets['UseInEnsemble'] = None

        # CIRRUS CORRECTION: Apply per-image cirrus-effect to checks and targets (for all filters):
        df_predictions_checks_targets = pd.merge(left=df_estimates_checks_targets,
                                                 right=self.df_cirrus_effect,
                                                 how='left', left_on='FITSfile', right_on='Image')
        # df_predictions_checks_targets.sort_values(by=['FOV', 'ModelStarID', 'Serial'],
        #                                           inplace=True)
        df_predictions_checks_targets.index = df_predictions_checks_targets['Serial']

        df_predictions_checks_targets['UntransformedMag'] = \
            df_predictions_checks_targets['PredictedMag'] - \
            df_predictions_checks_targets['CirrusEffect']

        # COLOR CORRECTION: interpolate Color Index values, then apply them (transform):
        transforms = {k: v.transform for (k, v) in self.skymodels.items()}  # a dict of transforms
        df_predictions_checks_targets = _impute_target_ci(df_predictions_checks_targets,
                                                          ci_filters=['V', 'I'],
                                                          transforms=transforms)
        df_transforms = pd.DataFrame([transforms], index=['Transform']).transpose()  # lookup table
        df = pd.merge(df_predictions_checks_targets, df_transforms,
                      how='left', left_on='Filter', right_index=True)
        df['TransformedMag'] = df['UntransformedMag'] - df['Transform'] * df['CI']
        df_transformed_without_errors = df.loc[~np.isnan(df['TransformedMag']), :]
        return df_transformed_without_errors

    def _compute_all_errors(self, df_transformed_without_errors):
        """
        Compute 3 error contributors and total sigma for each target and check star in each image. 
           model_sigma: from mixed-model regression (same for all observations in this filter).
           cirrus_sigma: from variance in ensemble comp stars (same for all obs in this image).
           inst_mag_sigma: from shot noise & background noise in specific obs (unique to each obs).
        :param: df_transformed_without_errors:
        :return: df_transformed, including errors [pandas DataFrame]
        """
        df_transformed = df_transformed_without_errors.copy()

        # Make new empty columns in df_transformed, to be populated later:
        #    (note: column InstMagSigma already exists from photometry)
        df_transformed['ModelSigma'] = np.float64(np.NaN)
        df_transformed['CirrusSigma'] = np.float64(np.NaN)
        df_transformed['TotalSigma'] = np.float64(np.NaN)
        for this_filter, skymodel in self.skymodels.items():
            images_this_filter = (self.df_comp_mags[self.df_comp_mags['Filter'] ==
                                                    this_filter])['FITSfile'].drop_duplicates()
            for image in images_this_filter:
                n = max(1, len(self.df_comp_mags[(self.df_comp_mags['FITSfile'] == image) &
                                                 (self.df_comp_mags['UseInEnsemble'] == True)]))
                model_sigma = skymodel.sigma / sqrt(n)
                cirrus_sigma = \
                    float(self.df_cirrus_effect.loc[self.df_cirrus_effect['Image'] == image,
                                                    'CirrusSigma'])
                # df_targets_checks = (df_transformed[df_transformed['FITSfile'] == image])\
                #     [['Serial', 'InstMagSigma']]  # USE NEXT LINE...(uses .loc[]
                df_targets_checks = df_transformed.loc[df_transformed['FITSfile'] == image,
                                                       ['Serial', 'InstMagSigma']]

                for serial in df_targets_checks['Serial']:
                    inst_mag_sigma = \
                        float(df_targets_checks.loc[df_targets_checks['Serial'] == serial,
                                                    'InstMagSigma'])
                    # Add up total error in quadrature:
                    total_sigma = sqrt(model_sigma**2 + cirrus_sigma**2 + inst_mag_sigma**2)
                    # Write new data into correct cells in df_transformed:
                    this_row = (df_transformed['Serial'] == serial)
                    df_transformed.loc[this_row, 'InstMagSigma'] = inst_mag_sigma
                    df_transformed.loc[this_row, 'ModelSigma'] = model_sigma
                    df_transformed.loc[this_row, 'CirrusSigma'] = cirrus_sigma
                    df_transformed.loc[this_row, 'TotalSigma'] = total_sigma

        # Finish forming df_transformed (return object):
        df_transformed = df_transformed.drop(['PredictedMag', 'Criterion1',
                                              'Criterion2', 'UntransformedMag', 'Transform'],
                                             axis=1)  # remove columns as a bit of cleanup
        df_columns_to_join = self.df_all_curated_obs[['Serial', 'FOV', 'MaxADU_Ur', 'FWHM',
                                                      'SkyADU', 'SkySigma']]  # for markup report
        df_transformed = pd.merge(left=df_transformed, right=df_columns_to_join, on='Serial')
        df_transformed.index = df_transformed['Serial']
        df_transformed = df_transformed.sort_values(by=['ModelStarID', 'JD_mid'])

        # ############# Start temporary code block: mimic R df:
        # r_columns = ['Serial', 'ModelStarID', 'FITSfile', 'StarID', 'Chart',
        #              'Xcentroid', 'Ycentroid',
        #              'InstMag', 'InstMagSigma', 'StarType', 'CatMag', 'CatMagError',
        #              'Exposure', 'JD_mid', 'Filter', 'Airmass', 'CI', 'SkyBias', 'Vignette',
        #              'LogADU', 'CirrusEffect', 'CirrusSigma', 'NumCompsUsed',
        #              'CompIDsUsed', 'NumCompsRemoved', 'JD_num', 'TransformedMag', 'ModelSigma',
        #              'TotalSigma', 'FOV', 'MaxADU_Ur', 'FWHM', 'SkyADU', 'SkySigma']
        # df_r = (df_transformed.copy())[r_columns]
        # df_r.to_csv(r'C:/24hrs/df_r.txt', sep=';', quotechar='"')  # '.txt', else Excel misbehaves
        # ############# End temporary code to mimic R df.

        return df_transformed

    def _write_summary_to_console(self):
        counts_by_star_type = self.df_transformed[['ModelStarID', 'StarType']]\
            .groupby(['StarType']).count()['ModelStarID']
        n_target_stars = len(self.df_transformed.loc[self.df_transformed['StarType'] == 'Target',
                                                     'ModelStarID'].unique())
        print('\nThis PredictionSet yields',
              str(counts_by_star_type['Target']), 'raw Target obs for',
              str(n_target_stars), 'targets and',
              str(counts_by_star_type['Check']), 'Check observations.\n')
        print('Now you are ready to:\n',
              '    1. (Optional) run '
              'df_master=pr.get_df_master(an_rel_directory=\'' + self.an_rel_directory +
              '\').sort_values([\'ModelStarID\', \'Filter\'])\n',
              '    2. (If stares) run ps.stare_comps(fov=\'XX Xxx\', star_id=\'\', '
              'this_filter=\'\')\n',
              '    3. run ps.markup_report(), then combine/reject target obs in'
              ' report_map.txt\n',
              '    4. run ps.aavso_report() and submit it to AAVSO.\n')

    def stare_comps(self, fov, star_id, this_filter):
        lines = get_stare_comps(self.df_transformed, fov, star_id, this_filter)
        print('\n'.join(lines))

    def stare_plot(self, star_id):
        # TODO: add this_filter parm and facility, default this-filter=None for all filters.

        # Setup data:
        df = (self.df_transformed.copy())
        df = df.loc[df['StarID'] == star_id, :]
        if len(df) == 0:
            print('This PredictionSet has no data for star_id \'' + star_id + '\'.')
            return
        floor_x = floor(df['JD_mid'].min())
        x = df['JD_mid'] - floor_x
        y = df['TransformedMag']
        color_dict = {'V': 'green', 'R': 'orange', 'I': 'red', 'B': 'blue'}
        colors = [color_dict[f] for f in df['Filter']]

        # Construct & draw stare plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(12, 8))  # (width, height) in "inches"
        ax.grid(True, color='lightgray', zorder=-1000)
        ax.set_title(star_id + '      ' + self.an_rel_directory,
                     color='darkblue', fontsize=20, weight='bold')
        ax.set_xlabel('JD(mid) - ' + str(floor_x))
        ax.set_ylabel('Best Mag')
        ax.scatter(x=x, y=y, color=colors, alpha=0.8, zorder=+1000)
        plt.gca().invert_yaxis()  # per custom of plotting magnitudes brighter=upward
        fig.canvas.set_window_title(star_id)
        plt.show()

    def markup_report(self):
        """
        Makes markup report from current PredictionSet object.
        :return: text string ready for printing (e.g., by copy/paste into Microsoft Word).
        """
        print('\nWriting \'markup_report.txt\'...', end='', flush=True)
        # First, make check-star dataframe:
        df = self.df_transformed.copy()
        df_check_stars = df.loc[df['StarType'] == 'Check',
                                ['FITSfile', 'StarID', 'TransformedMag', 'CatMag']]
        df_check_stars = df_check_stars.rename(columns={
            'StarID': 'Check', 'TransformedMag': 'CkMag', 'CatMag': 'CkCat'})

        # Make df for markup report:
        df = self.df_transformed.copy()
        df = df[df['StarType'] == 'Target']
        df = df[['Serial', 'FITSfile', 'StarID', 'Filter', 'Exposure', 'TransformedMag',
                 'InstMagSigma', 'ModelSigma', 'CirrusSigma', 'TotalSigma', 'MaxADU_Ur',
                 'FWHM', 'JD_num', 'FOV']]
        df = df.rename(columns={'StarID': 'Target', 'Exposure': 'Exp', 'TransformedMag': 'Mag',
                                'MaxADU_Ur': 'MaxADU'})
        df = pd.merge(df, df_check_stars, how='left', on='FITSfile')
        df.index = df['Serial']
        # df = df.sort_values(by=['Target', 'FITSfile', 'Filter', 'Exp'])
        # df = df.sort_values(by=['Target', 'FOV', 'JD_num'])
        # TODO: verify same order as for aavso_report.
        df = df.sort_values(by=['Target', 'JD_num'])  # same as for aavso_report (verify)

        # A nested helper function:
        def format_column(iterable, decimal_pts=None, min_width=0, left_pad=1):
            if decimal_pts is not None:
                this_list = [('{0:.' + str(decimal_pts) + 'f}').format(x) for x in iterable]
            else:
                this_list = [str(x) for x in iterable]
            n_chars = max(min_width, max([len(x) for x in this_list])) + left_pad
            return pd.Series([x.rjust(n_chars) for x in this_list])

        # Make dataframe df_text of text columns ~ ready to print (column FOV omitted from text):
        df_text = pd.DataFrame()
        df_text['Serial'] = format_column(df['Serial'], min_width=4)
        df_text['Target'] = format_column(df['Target'], min_width=6)
        df_text['FITSfile'] = format_column(df['FITSfile'], min_width=8)
        df_text['FITSfile'] = [f.split('.fts')[0]
                               for f in df_text['FITSfile']]  # remove '.fts'
        df_text['Filt'] = format_column(df['Filter'], min_width=4)
        df_text['Exp'] = format_column(df['Exp'], decimal_pts=1, min_width=5)
        df_text['Exp'] = [s[:-2] + '  ' if s.endswith('.0') else s
                          for s in df_text['Exp']]  # remove any trailing decimal point and zero
        df_text['Mag'] = format_column(df['Mag'], decimal_pts=3)
        df_text['MaxADU'] = format_column(round(df['MaxADU'].astype(int)), min_width=6)
        df_text['FWHM'] = format_column(df['FWHM'], decimal_pts=2, min_width=4)
        df_text['JD_fract'] = format_column(df['JD_num'], decimal_pts=4, min_width=6)
        df_text['Check'] = format_column(df['Check'], min_width=3)
        df_text['CkMag'] = format_column(df['CkMag'], decimal_pts=3, min_width=5)
        df_text['CkCat'] = format_column(df['CkCat'], decimal_pts=3, min_width=5)
        df_text['Inst'] = format_column(round(df['InstMagSigma']*1000.0).astype(int), min_width=3)
        df_text['Model'] = format_column(round(df['ModelSigma']*1000.0).astype(int), min_width=3)
        df_text['Cirr'] = format_column(round(df['CirrusSigma']*1000.0).astype(int), min_width=3)
        df_text['Sigma'] = format_column(round(df['TotalSigma']*1000.0).astype(int), min_width=3)

        # Make dict of just-sufficient column widths:
        column_list = ['Serial', 'Target', 'FITSfile', 'Filt', 'Exp', 'Mag', 'MaxADU',
                       'FWHM', 'JD_fract', 'Check', 'CkMag', 'CkCat',
                       'Inst', 'Model', 'Cirr', 'Sigma']
        dict_widths = {col: max(len(col), max([len(ss) for ss in df_text[col]]))
                       for col in column_list}

        # Make text lines of report:
        lines = ['MARKUP REPORT for ' + self.an_rel_directory +
                 '        generated by photrix ' + THIS_SOFTWARE_VERSION +
                 '  at  ' + '{:%Y-%m-%d %H:%M  UTC}'.format(datetime.now(timezone.utc)) +
                 '        ' + str(len(df)) + ' raw observations.']
        left_spacer = 2 * ' '
        nan_text = ' - '
        header_line = left_spacer + ' '.join([col.rjust(dict_widths[col]) for col in column_list])
        line_length = len(header_line)
        lines.extend(['', line_length * '_', '', header_line, ''])
        for i in range(len(df)):
            this_row = pd.Series([s if s.strip().lower() != 'nan' else nan_text
                                  for s in df_text.iloc[i]], index=column_list)
            line = left_spacer + ' '.join(this_row[col].rjust(dict_widths[col])
                                          for col in column_list)
            lines.append(line)
            # Add blank line between targets:
            if i < len(df) - 1:
                this_target = df_text.iloc[i]['Target']
                next_target = df_text.iloc[i+1]['Target']
                if next_target != this_target:
                    lines.append('')
        lines.extend(['', line_length * '_'])
        lines = [line + '\n' for line in lines]

        output_fullpath = os.path.join(self.an_top_directory, self.an_rel_directory,
                                       'Photometry', 'markup_report.txt')
        with open(output_fullpath, 'w') as this_file:
            this_file.write(''.join(lines))
        print('Done.\n' + 'Written to file \'' + output_fullpath + '\'.', flush=True)

    def aavso_report(self, write_file=True, return_df=False):
        """
        Construct AAVSO report (Enhanced Format) from this PredictionSet object.
        Writes this report as a text file in current (PredictionSet's) directory.
        :param write_file: True to write text file to current dir, False to not write.
        :param return_df: True to return a DataFrame of results, False to return None.
        :return: table of results if requested [DataFrame], else None.
        """

        df_report = self._apply_report_map_txt()
        df_report = df_report.sort_values(by=['TargetName', 'JD'])  # same as for markup_report (verify)

        # Construct text header lines:
        header = ["#TYPE=Extended",
                  "#OBSCODE=DERA",  # DERA = Eric Dose's observer code @ AAVSO
                  "#SOFTWARE=custom python scripts available at"
                  " https://github.com/edose/photrix, tag/version=" +
                  THIS_SOFTWARE_VERSION,
                  "#DELIM=" + AAVSO_REPORT_DELIMITER,
                  "#DATE=JD",
                  "#OBSTYPE=CCD",
                  "#This report of " + str(len(df_report)) + " observations was generated " +
                  '{:%Y-%m-%d %H:%M:%S UTC}'.format(datetime.now(timezone.utc)) +
                  " from raw data in directory " + self.an_rel_directory + ".",
                  "#Eric Dose, New Mexico Mira Project, ABQ, NM",
                  "#",
                  "#NAME,DATE,MAG,MERR,FILT,TRANS,MTYPE,CNAME,CMAG,KNAME,KMAG,AMASS," +
                  "GROUP,CHART,NOTES"
                  ]

        # Format observation report text fields, building df_formatted by columns left to right:
        df_formatted = pd.DataFrame()  # empty
        if len(df_report) == 0:
            obs_lines = ['\n\n\n              >>>>>>>>> NO OBSERVATIONS TO PRINT\n']
        else:
            df_formatted['TargetName'] = df_report['TargetName'].str.strip().str.upper()
            df_formatted['JD'] = df_report['JD'].astype(np.float64).map('{:.5f}'.format)
            df_formatted['Mag'] = df_report['Mag'].map('{:.3f}'.format)
            df_formatted['MagErr'] = df_report['TotalSigma'].map('{:.3f}'.format)
            df_formatted['Filter'] = df_report['Filter'].str.strip().str.upper()
            df_formatted['Transformed'] = 'YES'  # we always transform our reported data
            df_formatted['MType'] = 'STD'  # we use std comp stars, not "differential mode"
            df_formatted['CompName'] = df_report['CompName'].str.strip().str.upper()
            df_formatted['CompMag'] = ['{:.3f}'.format(mag) if not np.isnan(mag) else 'na'
                                       for mag in df_report['CompMag']]
            df_formatted['CheckName'] = [name.strip().upper() if isinstance(name, str) else 'na'
                                         for name in df_report['CheckName']]
            df_formatted['CheckMag'] = ['{:.3f}'.format(mag) if not np.isnan(mag) else 'na'
                                        for mag in df_report['CheckMag']]
            df_formatted['Airmass'] = df_report['Airmass'].map('{:.4f}'.format)
            df_formatted['Group'] = 'na'  # we don't use the observation grouping facility
            df_formatted['Chart'] = df_report['Chart'].str.upper()
            df_formatted['Notes'] = [note.strip() if note != '' else 'na'
                                     for note in df_report['Notes']]

            # Make all observation text lines:
            obs_column_list = ['TargetName', 'JD', 'Mag', 'MagErr', 'Filter', 'Transformed',
                               'MType', 'CompName', 'CompMag', 'CheckName', 'CheckMag',
                               'Airmass', 'Group', 'Chart', 'Notes']
            obs_lines = df_formatted.loc[:, obs_column_list].\
                apply(lambda x: AAVSO_REPORT_DELIMITER.join(x), axis=1).tolist()

        # Write file if requested:
        if write_file:
            lines_to_write = [line + '\n' for line in (header + obs_lines)]
            filename = 'AAVSOreport-' + self.an_rel_directory + '.txt'
            fullpath = os.path.join(self.an_top_directory, self.an_rel_directory,
                                    'Photometry', filename)
            with open(fullpath, 'w') as f:
                f.writelines(lines_to_write)
            print('AAVSO report for AN ' + self.an_rel_directory + ' written to: ' +
                  fullpath + '\n   = ' + str(len(df_formatted)) + ' reportable observations.')

        # Return DataFrame if requested:
        if return_df:
            return df_report
        return None

    def _apply_report_map_txt(self):
        """
        [Called only by other PredictionSet method.]
        :return: df_report: all observation data to construct AAVSO photometry report [DataFrame].
        """

        df_report = self.df_transformed.copy()
        df_report = df_report[['Serial', 'StarType', 'StarID', 'JD_mid',
                               'TransformedMag', 'TotalSigma', 'InstMagSigma',
                               'ModelSigma', 'CirrusSigma', 'Filter']]
        df_report.rename(columns={'StarID': 'TargetName', 'JD_mid': 'JD', 'TransformedMag': 'Mag'},
                         inplace=True)
        df_report['CompName'] = ''
        df_report['CompMag'] = np.float64(np.nan)
        df_report['NComps'] = 0
        # df_report['CheckName'] = ''
        # df_report['CheckMag'] = np.float64(np.nan)
        df_report['Airmass'] = self.df_transformed['Airmass']  # to get column order
        df_report['Chart'] = self.df_transformed['Chart']  # "
        df_report['Notes'] = ['obs#' + str(s) for s in df_report['Serial']]
        df_report = df_report[df_report['StarType'] == 'Target']
        df_report.drop('StarType', axis=1, inplace=True)  # remove column no longer needed

        # Get report_map.txt
        fullpath = os.path.join(self.an_top_directory, self.an_rel_directory, 'Photometry',
                                'report_map.txt')
        if not os.path.exists(fullpath):
            _write_report_map_stub(self.an_top_directory, self.an_rel_directory)
        with open(fullpath) as f:
            lines = f.readlines()
        lines = [line for line in lines if line is not None]  # remove empty list elements
        lines = [line.split(";")[0] for line in lines]  # remove all comments
        lines = [line.strip() for line in lines]  # remove lead/trail blanks
        lines = [line for line in lines if line != '']  # remove empty lines
        lines = [line for line in lines if line.startswith('#')]  # keep only directive lines

        # Apply #TARGET omissions:
        omit_lines = [line for line in lines if line.startswith('#TARGET')]
        for this_line in omit_lines:
            parms, warning_lines = _get_line_parms(this_line, '#TARGET', False, 1, 1)
            if warning_lines is not None:
                print('>>>>> Can\'t parse line:', warning_lines)
            else:
                target_to_omit = parms[0]
                if target_to_omit is not None:
                    rows_to_keep = [t.lower() != target_to_omit.lower()
                                    for t in df_report['TargetName']]
                    df_report = df_report[rows_to_keep]
                    print('Target removed:', target_to_omit, '.')

        # Apply #SERIAL omissions:
        rows_before = len(df_report)
        omit_lines = [line for line in lines if line.startswith('#SERIAL')]
        for this_line in omit_lines:
            parms, warning_lines = _get_line_parms(this_line, '#SERIAL', True, 1, None)
            if warning_lines is not None:
                print('>>>>> Can\'t parse line:', warning_lines)
            else:
                if parms is not None:
                    serials_to_omit = [int(p) for p in parms]
                    rows_to_keep = ~ df_report['Serial'].isin(serials_to_omit)

                    df_report = df_report[rows_to_keep]
        rows_after = len(df_report)
        print(str(rows_before - rows_after), 'serials removed altogether.')

        # Apply #JD directives:
        omit_lines = [line for line in lines if line.startswith('#JD')]
        for this_line in omit_lines:
            raw_parms, warning_lines = _get_line_parms(this_line, '#JD', True, 2, 2)
            if warning_lines is not None:
                print('>>>>> Can\'t parse line:', warning_lines)
            else:
                min_jd_fract, max_jd_fract = np.float64(raw_parms)
                if (min_jd_fract >= 0.0) and (max_jd_fract < 2.0):
                    floor_jd = floor(min(df_report['JD']))
                    rows_to_keep = ((df_report['JD'] < floor_jd + min_jd_fract)|\
                                     (df_report['JD'] > floor_jd + max_jd_fract))
                    rows_before_omit = len(df_report)
                    df_report = df_report[rows_to_keep]
                    rows_after_omit = len(df_report)
                    print('Omitted fractional JD range:', '{:.4f}'.format(min_jd_fract), 'to',
                          '{:.4f}'.format(max_jd_fract), '=',
                          str(rows_before_omit - rows_after_omit), 'observations.')

        # Add check-star names and mags  as new columns (look up from self.df_transformed):
        df_checks = ((self.df_transformed.copy())[self.df_transformed['StarType'] == 'Check'])\
            [['JD_mid', 'StarID', 'TransformedMag']]
        df_checks.rename(columns={'JD_mid': 'JD', 'StarID': 'CheckName',
                                  'TransformedMag': 'CheckMag'}, inplace=True)
        df_report = pd.merge(left=df_report, right=df_checks, how='left', on='JD')
        df_report.index = df_report['Serial']

        # Apply comp-star names and mags:
        df_comps = self.df_comp_mags[['Serial', 'JD_mid', 'StarID', 'Filter']].copy()
        df_comps['ObsMag'] = self.df_comp_mags['EstimatedMag']
        report_jds = df_report['JD'].drop_duplicates().sort_values()
        for this_jd in report_jds:
            df_comps_this_jd = df_comps[df_comps['JD_mid'] == this_jd]
            n_comps_this_jd = len(df_comps_this_jd)
            rows_to_update = (df_report['JD'] == this_jd)
            df_report.loc[rows_to_update, 'NComps'] = n_comps_this_jd
            if n_comps_this_jd <= 0:
                print('>>>>> No comp stars in model for JD = ', this_jd)
                return None
            if n_comps_this_jd > 1:
                # 'Ensemble' photometry case:
                df_report.loc[rows_to_update, 'CompName'] = 'ENSEMBLE'  # 'CompMag' remains NA.
                df_report.loc[rows_to_update, 'Notes'] += ' / ' + str(n_comps_this_jd) + ' comps'
            else:
                # Single comp-star case:
                df_report.loc[rows_to_update, 'CompName'] = (df_comps_this_jd['StarID']).iloc[0]
                df_report.loc[rows_to_update, 'CompMag'] = (df_comps_this_jd['ObsMag']).iloc[0]

        # Nested function for convenience:
        def all_same(list_or_series):
            this_list = list(list_or_series)
            return this_list.count(this_list[0]) == len(this_list)

        # Apply #COMBINE directives last (& verify check and comp stars are in fact associated):
        combine_lines = [line for line in lines if line.startswith('#COMBINE')]
        for this_line in combine_lines:
            raw_parms, warning_lines = _get_line_parms(this_line, '#COMBINE', True, 1, None)
            if warning_lines is not None:
                print('>>>>> Can\'t parse line:', warning_lines)
                continue
            else:
                serials_to_combine = [int(rp) for rp in raw_parms if int(rp) >= 1]
                df_combine = df_report[df_report['Serial'].isin(serials_to_combine)]
                # Verify that user-selected combine obs are in fact eligible to be combined:
                if len(df_combine) <= 1:
                    print('>>>>> Fewer than 2 obs to combine for line: \'' + this_line + '\'')
                    continue  # skip this #COMBINE line.
                if not all_same(df_combine['TargetName']):
                    print('>>>>> Non-uniform target names for line: \'' + this_line +
                          '\'...Combine is skipped.')
                    continue
                if not all_same(df_combine['Filter']):
                    print('>>>>> Non-uniform Filter for line: \'' + this_line +
                          '\'...Combine is skipped.')
                    continue
                if not all_same(df_combine['CompName']):
                    print('>>>>> Non-uniform Comp Names for line: \'' + this_line +
                          '\'...Combine is skipped.')
                    for i in range(len(df_combine)):
                        print(df_combine['Serial'].iloc[i],
                              df_combine['NComps'].iloc[i],
                              df_combine['CompName'].iloc[i])
                    continue
                real_check_names = [cn for cn in df_combine['CheckName'] if cn is not None]
                if len(real_check_names) >= 2:
                    if not all_same(real_check_names):
                        print('>>>>> Non-uniform Check Stars for line: \'' + this_line +
                              '\'...Combine is skipped.')
                        continue
                if not all_same(df_combine['Chart']):
                    print('>>>>> Non-uniform Chart IDs for line: \'' + this_line +
                          '\'...Combine is skipped.')
                    continue
                jd = df_combine['JD'].astype(np.float64)
                if max(jd) - min(jd) > 1.0 / 24.0:  # one hour
                    print('>>>>> Range of JD times is too large for line: \'' + this_line +
                          '\'...Combine is skipped.')
                    continue
                airmass = df_combine['Airmass'].astype(np.float64)
                if max(airmass) - min(airmass) > 0.400:
                    print('>>>>> Range of Airmasses is too large to combine for line: \'' +
                          this_line + '\'...Combine is skipped.')
                    continue

            # This #COMBINE line has passed all the tests, now execute it:
            df_new = df_combine.iloc[0:1].copy()  # a 1-row df; will combine rows into this.
            serial_to_replace = df_new['Serial'].iloc[0]  # scalar
            serials_to_delete = df_combine.loc[df_combine['Serial'] != serial_to_replace, 'Serial']
            # TODO: Consider weighted mean for combinations (w/limits on differences in weights?).
            df_new['JD'] = df_combine['JD'].astype(np.float64).mean()
            df_new['Mag'] = df_combine['Mag'].mean()
            n_combine = len(df_combine)
            # Instrument Mag sigmas are independent.
            inst_mag_sigma = sqrt(df_combine['InstMagSigma'].
                                  clip_lower(0.001).pow(2).mean()) / sqrt(n_combine)
            # Model sigma is uniform, but not independent and so not decreased by multiple images.
            model_sigma = df_combine['ModelSigma'].iloc[0]  # uniform across images to combine.
            # Cirrus_sigma ~ independent if no extreme outlier comp stars (which would correlate).
            cirrus_sigma = sqrt(df_combine['CirrusSigma'].
                                clip_lower(0.001).pow(2).mean()) / sqrt(n_combine)
            df_new['InstMagSigma'] = inst_mag_sigma
            df_new['ModelSigma'] = model_sigma
            df_new['CirrusSigma'] = cirrus_sigma
            df_new['TotalSigma'] = sqrt(model_sigma**2 + cirrus_sigma**2 + inst_mag_sigma**2)
            df_new['CheckMag'] = df_combine['CheckMag'].mean()
            df_new['NComps'] = df_combine['NComps'].min()
            df_new['Airmass'] = df_combine['Airmass'].mean()
            df_new['Notes'] = str(n_combine) + ' obs  >= ' + \
                str(int(df_new['NComps'].iloc[0])) + ' comps'
            df_report.update(df_new)
            df_report.drop(serials_to_delete, inplace=True)  # drop rows by index
            print('Combination of Serials ' +
                  ' '.join(df_combine['Serial'].astype(int).astype(str)) + ': done.')
        return df_report


END_PROCESSING_HERE____________ = ''


class TransformModel:
    def __init__(self, an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None,
                 filter=None, ci_type=None, fovs_to_include="All",
                 instrument_name='Borea', site_name='DSW',
                 max_cat_mag_error=0.03, max_inst_mag_sigma=0.03, max_ci=+2.5,
                 saturation_adu=None, fit_sky_bias=True,
                 fit_extinction=True, fit_log_adu=True):
        """
        Constructs a transform on filter against ci filter description ci_type,
            given a df_master.csv in the given directory.
        :param an_top_directory: e.g., 'J:\Astro\Images\C14' [string]
        :param an_rel_directory: e.g., '20170504'. The dir 'Photometry' is subdir of this. [string]
        :param filter: name of filter to which this model applies [string, e.g., 'V' or 'R']
        :param ci_type: two color-index filter names separated by minus sign,
            e.g., 'V-I' for V-I index [string]
        :param fovs_to_include: defines which eligible rows (by FOV) will be included:
            Choices:
                "All" will use all eligible rows in df_master;
                "Standards" will use all eligible rows from all Standard FOVs in df_master;
                One FOV name as a string will use that one FOV only;
                A list of strings will use those FOV(s) only.
        :param instrument_name: name of Instrument, e.g., 'Borea' [string; name of Instrument obj]
        :param site_name: name of observing site, e.g., 'DSW' [string; name of Site object]
        :param max_cat_mag_error: maximum catalog error allowed to stars in model [float]
        :param max_inst_mag_sigma: max instrument magnitude error allowed star observations [float]
        :param max_ci: maximum color index allowed to stars in model [float]
        :param saturation_adu: ccd ADUs that constitute saturation [float; None if from Instrument]
        :param fit_sky_bias: True to fit sky bias term [bool]
        :param fit_log_adu: True to fit log(MaxADU_Ur) as CCD nonlinearity measure [bool]
        :param fit_extinction: True to fit extinction terms; else use values from Site obj.
            Used only if image_name==False and more than all images are used [bool]
        """
        self.an_top_directory = an_top_directory
        self.an_rel_directory = an_rel_directory
        self.filter = filter
        self.ci_filters = [s.strip().upper() for s in (ci_type.strip().split('-'))]
        if len(self.ci_filters) != 2:
            print(' >>>>> Invalid ci_type \'' + ci_type + '\'')
            self.is_valid = False
            return
        self.fovs_to_include = fovs_to_include
        self.instrument_name = instrument_name
        self.site_name = site_name
        self.max_cat_mag_error = max_cat_mag_error
        self.max_inst_mag_sigma = max_inst_mag_sigma
        self.max_ci = max_ci
        if saturation_adu is not None:
            self.saturation_adu = saturation_adu
        else:
            instrument = Instrument(self.instrument_name)
            self.saturation_adu = instrument.camera['saturation_adu']
        self.fit_sky_bias = fit_sky_bias
        self.fit_extinction = fit_extinction
        self.fit_log_adu = fit_log_adu
        self.dep_var_name = 'InstMag_with_offsets'
        self.df_model = None
        self.fov_list = None
        self.image_list = None
        self.statsmodels_object = None
        self.mm_fit_object = None  # populated only if >=2 images and mixed model used.
        self.fitted_values = None
        self.param_values = None
        self.residuals = None
        self.sigma = None
        self.image_effect = None
        self.is_valid = False  # default until object constructed.

        # Execute steps:
        self._make_and_curate_model_dataframe(fovs_to_include=self.fovs_to_include)
        self._add_ci_column()  # to df_model, add new column CI (color index) from catalog mags.
        self._prep_and_do_regression()
        self._build_output()

    def _make_and_curate_model_dataframe(self, fovs_to_include):
        df, warning_lines = _apply_omit_txt(self.an_top_directory, self.an_rel_directory)
        master_fov_list = df['FOV'].drop_duplicates().tolist()

        # Interpret fovs_to_include, ultimately yielding fov_list:
        if isinstance(fovs_to_include, list):
            fov_list = [fn for fn in master_fov_list if fn in fovs_to_include]
        elif isinstance(fovs_to_include, str):
            if fovs_to_include.lower() == 'all':
                fov_list = master_fov_list.copy()
            elif fovs_to_include.lower() in ['standard', 'standards']:
                fov_list = [fn for fn in master_fov_list
                            if Fov(fn).target_type.lower() == 'standard']
            else:
                fov_list = [fn for fn in master_fov_list if fn in [fovs_to_include]]
        else:
            print("Couldn't interpret input parm 'fovs_to_include'.")
            self.is_valid = False
            return
        if len(fov_list) <= 0:
            print("Input parm 'fovs_to_include' yielded no fovs for folder '",
                  self.an_rel_directory, "'.")
            self.is_valid = False
            return

        # Retain only the appropriate rows of raw dataframe:
        df = df[df['FOV'].isin(fov_list)]
        df = df[df['Filter'] == self.filter]
        df = df[df['StarType'].isin(['Comp', 'Check'])]

        # Retain only the needed columns of dataframe:
        df = df[['Serial', 'StarID', 'FITSfile', 'InstMagSigma', 'FWHM', 'MaxADU_Ur', 'StarType',
                 'JD_mid', 'Filter', 'Airmass', 'FOV', 'CatMag', 'CatMagError', 'InstMag',
                 'SkyBias', 'LogADU']]

        # Curate rows for various quality measures:
        df = df[df['CatMag'].notnull()]
        df = df[df['Airmass'].notnull()]
        df = df[df['InstMagSigma'] <= self.max_inst_mag_sigma]
        df = df[df['MaxADU_Ur'].notnull()]
        df = df[df['MaxADU_Ur'] <= self.saturation_adu]
        df = df[df['CatMagError'].notnull()]
        df = df[df['CatMagError'] <= self.max_cat_mag_error]
        df = df[df['FWHM'] >= MIN_FWHM]
        self.df_model = df

    def _add_ci_column(self):
        df = self.df_model.copy()  # local working copy
        df['CI'] = None  # new blank column
        fov_names = df['FOV'].drop_duplicates()
        for fov_name in fov_names:
            df_fov = (df.copy())[df['FOV'] == fov_name]
            this_fov = Fov(fov_name)
            fov_stars = this_fov.aavso_stars  # a list of AavsoSequenceStar_WithMagError objects
            star_ids = df_fov['StarID']
            for star_id in star_ids:
                # Extract proper AavsoSS_WME object & mags for this star (a dict):
                star_list = [fs for fs in fov_stars if (fs.star_id == star_id and fs.is_valid)]
                if len(star_list) == 1:
                    fov_star = star_list[0]  # object found.
                    # fov_star.mags is tuple (mag, mag_err)
                    absent = (None, None)
                    ci_mag_first = fov_star.mags.get(self.ci_filters[0], absent)[0]
                    ci_mag_second = fov_star.mags.get(self.ci_filters[1], absent)[0]
                else:
                    ci_mag_first, ci_mag_second = (None, None)
                if (ci_mag_first is not None) and (ci_mag_second is not None):
                    this_ci = ci_mag_first - ci_mag_second
                else:
                    this_ci = None
                if this_ci is not None:
                    rows_to_update = (df['FOV'] == fov_name) & (df['StarID'] == star_id)
                    df.loc[rows_to_update, 'CI'] = this_ci
        rows_with_ci = [(ci is not None) for ci in df['CI']]
        self.df_model = df[rows_with_ci]
        self.fov_list = self.df_model['FOV'].drop_duplicates().tolist()
        self.image_list = self.df_model['FITSfile'].drop_duplicates().tolist()

    def _prep_and_do_regression(self):
        if len(self.image_list) <= 0:
            print("No images in image list.")
            self.is_valid = False
            return

        # Build variable list & dep-var offset:
        x_var_list = []  # fixed-effects only
        dep_var_offset = self.df_model['CatMag'].copy()  # *copy* CatMag, or df_model risks damage.
        if len(self.image_list) >= 2:  # if one image, do nothing (fold extinction into zero-point).
            if self.fit_extinction:
                x_var_list.append('Airmass')
            else:
                site = Site(self.site_name)
                extinction = site.extinction[self.filter]
                dep_var_offset += extinction * self.df_model['Airmass']
        if self.fit_sky_bias:
            if sum([x != 0 for x in self.df_model['SkyBias']]) > int(len(self.df_model) / 2):
                x_var_list.append('SkyBias')
        if self.fit_log_adu:
            x_var_list.append('LogADU')
        x_var_list.append('CI')  # to get transform value (the point of this class).
        for var in x_var_list:
            self.df_model[var] = np.float64(self.df_model[var])  # ensure float64 for statsmodels.

        # Build the regression model (ordinary least-squares or mixed-model):
        self.df_model['DepVar'] = self.df_model['InstMag'] - dep_var_offset  # dependent variable.
        if len(self.image_list) == 1:
            # Model with ordinary least squares (OLS, as no group/image effects possible):
            formula = 'DepVar ~ ' + ' + '.join(x_var_list)
            self.statsmodels_object = smf.ols(formula, data=self.df_model).fit()
        else:
            # Model with photrix.util.MixedModelFit (random effect is per-image, "cirrus effect"):
            random_effect_var_name = 'FITSfile'  # cirrus effect is per-image
            self.mm_fit_object = MixedModelFit(data=self.df_model, dep_var='DepVar',
                                               fixed_vars=x_var_list,
                                               group_var=random_effect_var_name)
            self.statsmodels_object = self.mm_fit_object.statsmodels_object

        if self.statsmodels_object.scale != 0.0 and \
                self.statsmodels_object.nobs == len(self.df_model):
            print(self.statsmodels_object.summary())
            self.is_valid = True

    def _build_output(self):
        if len(self.image_list) == 1:
            # Ordinary least squares case...image/cirrus effect folded into zero-point term:
            so = self.statsmodels_object
            self.fitted_values = so.fittedvalues
            self.param_values = so.params
            self.residuals = so.resid
            self.sigma = so.mse_resid ** 0.5
            self.transform_value = so.params['CI']
            self.transform_sigma = so.bse['CI']
            self.image_effect = None
        else:
            # Mixed-model case, data from util.MixedModelFit object:
            so = self.mm_fit_object
            self.fitted_values = so.df_observations['FittedValue']
            self.param_values = so.df_fixed_effects['Value']
            self.residuals = so.df_observations['Residual']
            self.sigma = so.sigma
            self.transform_value = self.param_values['CI']
            self.transform_sigma = so.df_fixed_effects.loc['CI', 'Stdev']
            self.image_effect = so.df_random_effects['GroupValue']


def get_df_master(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Simple utility to read df_master.csv file and return the DataFrame.
    :param an_top_directory: path to an_rel_directory [str]
    :param an_rel_directory: directory for this instrument on this Astronight [string] 
    :return: pandas DataFrame with all comp, check, and target star raw photometric data
         (which dataframe is generated and csv-written by R, as of May 2017).
         The DataFrame index is set to equal column Serial (which was already a kind of index).
    """
    if an_rel_directory is None:
        return None
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', DF_MASTER_FILENAME)
    if not os.path.exists(fullpath):
        return None
    df_master = pd.read_csv(fullpath, sep=';')
    df_master.index = df_master['Serial']
    return df_master


def _add_ci_values(df_fov, df_star_data_numbered, instrument):
    # Killswitch is next line:
    return df_fov

    df_fov['CI_type'] = ''
    df_fov['CI_value'] = np.nan
    df_by_filter = df_fov.groupby('Filter')
    for this_filter, df in df_by_filter.groups:
        transforms = instrument.transforms(this_filter)
        # Decide best CI type for this filter:
        if len(instrument.transforms(this_filter)) == 1:
            # Easy case: only one transform for this filter, so use it:
            ci_type = instrument.transforms(this_filter)[0]
            ci_filters = [s.strip().upper() for s in (ci_type.strip().split('-'))]
            ci_filter_1 = ci_filters[0]
            ci_filter_2 = ci_filters[1]
        else:
            # Hard case: have to decide which transform to use for this FOV and filter:
            for this_transform in instrument.transforms(this_filter):
                this_ci_type, this_ci_value = this_transform  # unpack tuple
                these_ci_filters = [s.strip().upper() for s in (this_ci_type.strip().split('-'))]
                this_ci_filter_1 = these_ci_filters[0]
                this_ci_filter_2 = these_ci_filters[1]
                # Score for this filter: how many target obs satisfy these two criteria:
                # Criterion 1:
                # Criterion 2:

        # Fill in CI_type and CI_value (CI_type decided above):
        # CI type described by ci_filter_1 and ci_filter_2
        star_ids = df['StarID'].unique()
        for star_id in star_ids:
            mags = df_star_data_numbered.loc[df_star_data_numbered['StarID'] == star_id, 'Mags']
            ci_mag_1 = mags[ci_filter_1][0]
            ci_mag_2 = mags[ci_filter_2][0]
            df.loc[df_star_data_numbered['StarID'] == star_id, 'CI_type'] = ci_type
            df.loc[df_star_data_numbered['StarID'] == star_id, 'CI_value'] = ci_mag_1 - ci_mag_2
    return df_fov


def get_stare_comps(df_transformed, fov=None, star_id=None, this_filter=None):
    df = df_transformed.copy()
    df = df[df['FOV'] == fov]
    df = df[df['StarID'] == star_id]
    df = df[df['Filter'] == this_filter]
    joined_comp_stars = ','.join(df['CompIDsUsed'])
    comp_stars_available = pd.Series(joined_comp_stars.split(',')).drop_duplicates()
    df_comps = pd.DataFrame({'CompID': comp_stars_available, 'IsIncluded': False})
    result_lines = ["EDIT file 'pre-predict' with one of the following lines:"]
    if len(df) <= 1:
        return result_lines + \
               ['   >>> One or zero qualifying images in dataframe.']
    for num_to_test in range(1, 1 + len(df_comps)):
        base = list((df_comps[df_comps['IsIncluded']])['CompID'])
        test = list((df_comps[~ df_comps['IsIncluded']])['CompID'])
        max_images_qualifying = 0
        best_test = test[0]  # to give a default so comparisons don't fail
        for this_test in test:
            num_images_qualifying = 0
            test_comps = base.copy()
            test_comps.append(this_test)
            for i_row in range(len(df)):
                this_set = (df.iloc[i_row])['CompIDsUsed'].split(',')
                if this_test in this_set:
                    num_images_qualifying += 1
            if num_images_qualifying > max_images_qualifying:
                # TODO: break ties by choosing lower catmagerror (or initially sorting on them)?
                best_test = this_test
                max_images_qualifying = num_images_qualifying
        df_comps.loc[(df_comps['CompID'] == best_test), 'IsIncluded'] = True
        this_line = '   ' + str(sum(df_comps['IsIncluded'])) + ' comps -> ' + \
                    str(max_images_qualifying) + ' images qualify  -->   #COMPS ' + fov + \
                    ', ' + this_filter + ', ' + \
                    ', '.join(df_comps.loc[df_comps['IsIncluded'], 'CompID'])
        result_lines.append(this_line)
    return result_lines


def _rename_to_photrix(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    # Construct DataFrame of files to rename:
    fits_path = os.path.join(an_top_directory, an_rel_directory, 'Uncalibrated')
    ur_names, object_list, jd_mid_list, filter_list, photrix_names = [], [], [], [], []
    for entry in os.scandir(fits_path):
        if entry.is_file():
            ur_names.append(entry.name)
            this_fits = FITS(an_top_directory, os.path.join(an_rel_directory, 'Uncalibrated'),
                             entry.name)
            object_list.append(this_fits.object)
            jd_mid_list.append(jd_from_datetime_utc(this_fits.utc_mid))
            filter_list.append(this_fits.filter)
    df = pd.DataFrame({'UrName': ur_names, 'Object': object_list,
                       'JD_mid': jd_mid_list, 'Filter': filter_list})
    df = df.sort_values(by=['Object', 'JD_mid'])

    # Construct new photrix names and add them to DataFrame:
    serial_number = 1
    for i in range(len(df)):
        this_object = df['Object'].iloc[i]
        this_filter = df['Filter'].iloc[i]
        if i >= 1:
            if this_object != df['Object'].iloc[i-1]:
                serial_number = 1
            else:
                serial_number += 1
        photrix_name = '-'.join([this_object, '{:04d}'.format(serial_number), this_filter]) + '.fts'
        photrix_names.append(photrix_name)
    df['PhotrixName'] = photrix_names
    df.index = photrix_names

    # Rename all the FITS files:
    for old_name, new_name in zip(df['UrName'], df['PhotrixName']):
        old_path = os.path.join(an_top_directory, an_rel_directory, 'Uncalibrated', old_name)
        new_path = os.path.join(an_top_directory, an_rel_directory, 'Uncalibrated', new_name)
        os.rename(old_path, new_path)

    # Write renaming table to Photometry subdirectory as csv file:
    renaming_fullpath = os.path.join(an_top_directory, an_rel_directory,
                                     'Photometry', 'File-renaming.txt')
    df.to_csv(renaming_fullpath, sep=';')


def _set_fits_extensions(an_top_directory, fits_subdir, fits_filenames=None):
    from photrix.image import FITS_EXTENSIONS
    target_ext = '.fts'
    for filename in fits_filenames:
        f, ext = os.path.splitext(filename)
        if (ext.replace('.', '') in FITS_EXTENSIONS) and (ext != target_ext):
            old_fullpath = os.path.join(an_top_directory, fits_subdir, filename)
            new_fullpath = os.path.join(an_top_directory, fits_subdir, f + target_ext)
            if not os.path.exists(new_fullpath):
                os.rename(old_fullpath, new_fullpath)


def _archive_fov_files(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, fov_names=None):
    from_dir = FOV_DIRECTORY
    to_dir = os.path.join(an_top_directory, an_rel_directory, 'FOV')
    for fov_name in fov_names:
        from_fullpath = os.path.join(from_dir, fov_name + '.txt')
        to_fullpath = os.path.join(to_dir, fov_name + '.txt')
        shutil.copy2(from_fullpath, to_fullpath)


def _apply_omit_txt(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Gets df_master and omit.txt and returns filtered DataFrame.
    :param an_top_directory: [string]
    :param an_rel_directory: [string]
    :return: rows of df_master whose omission is NOT requested by user (in omit.txt) [DataFrame] 
    """
    if an_rel_directory is None:
        return None
    df_master = get_df_master(an_top_directory=an_top_directory, an_rel_directory=an_rel_directory)
    if df_master is None:
        return None

    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    if not os.path.exists(fullpath):
        _write_omit_txt_stub(an_top_directory, an_rel_directory)
        return df_master.copy(), []  # no change or warnings, since omit.txt absent

    with open(fullpath) as f:
        lines = f.readlines()
    lines = [line.split(";")[0] for line in lines]  # remove all comments
    lines = [line.strip() for line in lines]        # remove lead/trail blanks
    lines = [line for line in lines if line != '']  # remove empty lines
    df_filtered = df_master.copy()  # start with a copy, omit lines per user requests.
    warning_lines = []

    for line in lines:
        warning_line = None
        rows_to_omit = len(df_filtered) * [False]  # default to be overwritten
        if line.startswith('#OBS'):
            parms, warning_line = _get_line_parms(line, "#OBS", False, 2, 2)
            if parms is not None:
                fits_file_name = parms[0] + '.fts'
                star_id = parms[1]
                rows_to_omit = (df_filtered['FITSfile'] == fits_file_name) & \
                               (df_filtered['StarID'] == star_id)  # a pandas Series
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#STAR'):
            parms, warning_line = _get_line_parms(line, "#STAR", False, 2, 3)
            if parms is not None:
                fov = parms[0]
                star_id = parms[1]
                rows_to_omit = (df_filtered['FOV'] == fov) & (df_filtered['StarID'] == star_id)
                if len(parms) == 3:
                    # further restrict omission request
                    filter = parms[2]
                    rows_to_omit = rows_to_omit & (df_filtered['Filter'] == filter)
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#SERIAL'):
            parms, warning_line = _get_line_parms(line, "#SERIAL", True, 1, None)
            if parms is not None:
                serials_to_omit = [int(p) for p in parms]
                rows_to_omit = df_filtered['Serial'].isin(serials_to_omit)
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#IMAGE'):
            parms, warning_line = _get_line_parms(line, "#IMAGE", False, 1, 1)
            if parms is not None:
                image_to_omit = parms[0] + '.fts'
                rows_to_omit = df_filtered['FITSfile'] == image_to_omit
            else:
                warning_lines.append(warning_line)
        elif line.startswith('#JD'):
            parms, warning_line = _get_line_parms(line, "#JD", True, 2, 2)
            if parms is not None:
                jd_floor = floor(min(df_filtered['JD_mid']))
                jd_start = float(parms[0]) + jd_floor
                jd_end = float(parms[1]) + jd_floor
                rows_to_omit = (df_filtered['JD_mid'] >= jd_start) & \
                    (df_filtered['JD_mid'] <= jd_end)
            else:
                warning_lines.append(warning_line)
        else:
            warning_line = 'Directive not understood: \'' + line + '\'.'
            warning_lines.append(warning_line)

        if sum(rows_to_omit) >= 1:
            df_filtered = df_filtered[~ rows_to_omit]  # remove rows as user requested.
        else:
            if warning_line is None:
                warning_lines.append('No rows omitted: \'' + line + '\'.')

    for warning_line in warning_lines:
        print(warning_line)
    print(str(len(df_master) - len(df_filtered)) + ' rows removed via omit.txt. ' +
          str(len(df_filtered)) + ' rows remain.')
    return df_filtered, warning_lines


def _write_omit_txt_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Will NOT overwrite existing omit.txt.
    """
    lines = [';----- This is omit.txt for AN directory ' + an_rel_directory,
             ';----- Use this file to omit observations from input to SkyModel (all filters).',
             ';----- Example directive lines:',
             ';',
             ';#OBS    Obj-0000-V, 132 ; to omit star 132 from FITS image Obj-0000-V.fts',
             ';#STAR   FOV, 132, V     ; to omit star 132 from all FITS with FOV '
             'and filter V',
             ';#STAR   FOV, 132        ; to omit star 132 from all FITS with FOV '
             'and ALL filters',
             ';#IMAGE  Obj-0000-V      ; to omit FITS image Obj-0000-V.fts specifically',
             ';#JD     0.72, 1         ; to omit fractional JD from 0.72 through 1',
             ';#SERIAL 123,77 54   6   ; to omit observations by Serial number (many per line OK)',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'omit.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written


def _curate_stare_comps(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None, df_in=None):
    """
    Using user's stare_comps.txt in this an_rel_directory,
       remove unwanted (stare) comp observations from further use.
    :return: data for all observations remaining eligible after this curation [DataFrame],
                and warning messages [list of strings]
    """
    if df_in is None:
        return None, None

    # Read & parse control file stare_comps.txt:
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'stare_comps.txt')
    if not os.path.exists(fullpath):
        _write_stare_comps_txt_stub(an_top_directory, an_rel_directory)
        return None  # no change or warnings as stare_comps.txt absent
    with open(fullpath) as f:
        lines = f.readlines()
    lines = [line.split(";")[0] for line in lines]  # remove all comments
    lines = [line.strip() for line in lines]  # remove lead/trail blanks
    lines = [line for line in lines if line != '']  # remove empty lines

    # Apply directive lines from stare_comps.txt:
    df = df_in.copy()   # starting point
    warning_lines = []  # "
    for line in lines:
        warning_line = None
        if line.startswith("#COMPS"):
            # Do not split on spaces, because FOV name could contain spaces.
            parms, warning_line = _get_line_parms(line, "#COMPS", False, 3, None)
            if parms is not None:
                fov_name = parms[0]
                filter_name = parms[1]
                comp_ids = ' '.join(parms[2:]).split()  # parse whether comma- or space-separated
                rows_to_remove = (df['StarType'] == 'Comp') & \
                                 (df['FOV'] == fov_name) & \
                                 (df['Filter'] == filter_name) & \
                                 (~ df['StarID'].isin(comp_ids))
                df = df[~ rows_to_remove]
            else:
                warning_lines.append(warning_line)
        else:
            warning_line = 'Directive not understood: \'' + line + '\'.'
            warning_lines.append(warning_line)
    return df, warning_lines


def _write_stare_comps_txt_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Will NOT overwrite existing stare_comps.txt file.
    """
    lines = [';----- This is stare_comps.txt for AN directory ' + an_rel_directory,
             ';----- Select comp stars (by FOV, filter, StarID) from input to '
             'rerun of PredictionSet() ',
             ';----- Example directive line:',
             ';',
             ';#COMPS  Obj, V, 132, 133 144    ; to KEEP from FOV \'Obj\': '
             'comp stars \'132\' \'133\' and \'144\' in filter \'V\'',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'stare_comps.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written


def _write_report_map_stub(an_top_directory=AN_TOP_DIRECTORY, an_rel_directory=None):
    """
    Will NOT overwrite existing report.map file.
    """
    lines = [';----- This is report_map.txt for AN directory ' + an_rel_directory,
             ';----- Use this file to omit and/or combine target observations, for AAVSO report',
             ';----- Example directive line:',
             ';',
             ';#TARGET  GD Cyg    ; to omit this target star altogether from AAVSO report',
             ';#JD      0.65 0.67 ; to omit this JD (fractional) range from AAVSO report',
             ';#SERIAL  34 44,129  32 1202 ; to omit these 5 Serial numbers from AAVSO report',
             ';#COMBINE 80 128    ; to combine (average) these 2 Serial numbers in AAVSO report',
             ';',
             ';----- Add your directive lines:',
             ';']
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry', 'report_map.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written


def _get_line_parms(line, directive, sep_by_spaces=False, nparms_min=None, nparms_max=None):
    """
    Take directive line and return parms and possible warning text lines.
    :param line: the line to parse [string]
    :param directive: the directive [string, e.g., '#SERIAL']
    :param sep_by_spaces: False if parms may include spaces [e.g., '#IMAGE' and '#TARGET',
        which may include spaces, e.g., 'ST Tri';
        True if spaces may separate parms [e.g., '#SERIAL' and '#COMBINE' directives]  [bool]
    :param nparms_min:  min number of parms to accept, or None to ignore [int]
    :param nparms_max:  max number of parms to accept, or None to ignore [int]
    :return: 2-tuple of parm list [list of strings] and warning lines [list of strings]
    """
    directive_text = line.split(';')[0].strip()  # remove any comment text
    if directive_text.startswith(directive) is False:
        return None, 'Line does not begin with correct directive: \'' + line + '\'.'
    line_parms = [p.strip() for p in directive_text[(len(directive)):].split(',')]
    if sep_by_spaces is True:
        line_parms = ' '.join(line_parms).split()
    valid_num_line_parms = True  # until falsified
    if nparms_min is not None:  # ignore nparms_min if None.
        valid_num_line_parms = valid_num_line_parms & (len(line_parms) >= nparms_min)
    if nparms_max is not None:  # ignore nparms_max if None.
        valid_num_line_parms = valid_num_line_parms & (len(line_parms) <= nparms_max)
    if valid_num_line_parms:
        return line_parms, None
    else:
        return None, 'Line has wrong number of parameters: \'' + line + '\'.'


def _impute_target_ci(df_predictions_checks_targets, ci_filters, transforms):
    """
    Impute Color Index value for each target and check star, by time-interpolation from known
       (comp) Color Index values. This will REPLACE CI values for Target and Check stars
       (which probably had been set to zero for targets but catalog CI for checks).
       CALLS _extract_ci_points() to get list of ci_filter observations to interpolate.
    :param df_predictions_checks_targets: [pandas DataFrame] 
    :param ci_filters: ['V', 'I'] for the time being (May 2017).
    :param transforms: transforms for all skymodels [dict filter:transform(V-I)]
    :return: updated_df_predictions_checks_targets: updated with CI color index for targets 
    """
    jd_floor = floor(min(df_predictions_checks_targets['JD_mid']))
    df_predictions_checks_targets['JD_num'] = df_predictions_checks_targets['JD_mid'] - jd_floor

    # Only target & check stars:
    target_and_check_rows = df_predictions_checks_targets['StarType'].isin(['Target', 'Check'])
    star_ids_targets_checks = \
        ((df_predictions_checks_targets.copy())[target_and_check_rows])['ModelStarID']
    # Sorted for convenience in testing:
    star_ids_targets_checks = star_ids_targets_checks.drop_duplicates().sort_values()

    # Replace CI color index values for one target star at a time:
    for this_star_id in star_ids_targets_checks:
        rows_this_star_id = df_predictions_checks_targets['ModelStarID'] == this_star_id
        df_star_id = (df_predictions_checks_targets[rows_this_star_id])\
            [['Serial', 'ModelStarID', 'Filter', 'JD_num', 'CI', 'UntransformedMag']]\
            .sort_values(by='JD_num')
        df_ci_points = _extract_ci_points(df_star_id, ci_filters, transforms)
        if len(df_ci_points) <= 0:
            df_predictions_checks_targets.loc[rows_this_star_id, 'CI'] = None
            print(">>>>> ModelStarID=", this_star_id,
                  " no CI points returned by input_target_ci()")
            continue
        df_ci_points = df_ci_points.sort_values(by='JD_num')

        # Interpolate CI, put values into df_star_id:
        #    (Interpolation method depends on number of interpolation points avaialable.)
        if len(df_ci_points) == 1:  # normal case for LPVs
            df_star_id['CI'] = df_ci_points.loc[0, 'CI']  # 1 point --> all CIs set to single value
        elif len(df_ci_points) in [2, 3]:  # 2 or 3 points --> linear fit of CI vs time
            x = df_ci_points['JD_num']
            y = df_ci_points['CI']
            this_linear_fit = np.polyfit(x, y, 1)  # (x,y,deg=1 thus linear)
            this_fit_function = np.poly1d(this_linear_fit)
            df_star_id['CI'] = this_fit_function(df_star_id['JD_num'])  # does the linear interpol
            # Enforce no extrapolation:
            stars_before_jd_range = df_star_id['JD_num'] < min(x)
            stars_after_jd_range = df_star_id['JD_num'] > max(x)
            df_star_id.loc[stars_before_jd_range, 'CI'] = this_fit_function(min(x))
            df_star_id.loc[stars_after_jd_range, 'CI'] = this_fit_function(max(x))
        else:  # here, 4 or more CI points to use in interpolation (prob a stare)
            # TODO: should this be interpolation in V mag rather than (or as alt to) interp in time?
            x = df_ci_points['JD_num']
            y = df_ci_points['CI']
            weights = len(x) * [1.0]
            smoothness = len(x) * 0.03**2  # i.e, N * (sigma_CI)**2
            # Construct spline; ext=3 -> no extrapolation, rather fixed at boundary values:
            spline = UnivariateSpline(x=x, y=y, w=weights, s=smoothness, ext=3)
            df_star_id['CI'] = spline(df_star_id['JD_num'])

        # Insert CI values for this target star:
        indices_to_update = [x
                             if x in df_predictions_checks_targets['Serial'] else None
                             for x in df_star_id['Serial']]
        df_predictions_checks_targets.loc[indices_to_update, 'CI'] = df_star_id['CI']

    return df_predictions_checks_targets


def _extract_ci_points(df_star_id, ci_filters, transforms):
    """
    Derive credible raw Color Index points from continguous observations in the Color Index filters.
    :param df_star_id: rows from df_predictions holding one ModelStarID [pandas DataFrame] 
    :param ci_filters: exactly two filters, in order, defining Color Index [e.g., ['V','I']]
    :param transforms: transforms in the ci_filters base, from which the relevant transform is
               extracted [dict, e.g., {'V': 0.025, 'I': -0.044} ]
    :return: very small dataframe of JD and Color Index for this star.
    """
    max_diff_jd = 60.0 / (24 * 60)  # 60 minutes in days; max time between adjacent obs for color
    rows_to_keep = df_star_id['Filter'].isin(ci_filters)
    df = (df_star_id.copy())[rows_to_keep].sort_values(by='JD_num')

    ci_point_list = []
    if len(df) <= 1:
        return pd.DataFrame()  # can't extract realistic pairs if there's only one point
    # Extract a color (one filter's magnitude minus the other's mag) whenever adjacent observations
    #    of the same ModelStarID are in different filters:
    filter = list(df['Filter'])
    jd = list(df['JD_num'])
    u_mag = list(df['UntransformedMag'])
    for ind in range(0, len(df) - 1):
        if filter[ind+1] != filter[ind]:  # if adj obs differ in filter...
            if (jd[ind+1] - jd[ind]) < max_diff_jd:  # ...and if not too different in time
                new_point_jd = (jd[ind+1] + jd[ind]) / 2.0
                new_point_ci = \
                    _solve_for_real_ci(
                        untransformed_mags={filter[ind]: u_mag[ind], filter[ind+1]: u_mag[ind+1]},
                        ci_filters=ci_filters,
                        transforms=transforms)
                ci_point_list.append({'JD_num': new_point_jd, 'CI': new_point_ci})  # dict for 1 row
    df_ci_point = pd.DataFrame(ci_point_list)
    return df_ci_point


def _solve_for_real_ci(untransformed_mags, ci_filters, transforms):
    """
    Solves for best estimate of real (transformed) Color Index.
    Dicts are passed (unordered); this function safely deduces the sign of the Color Index.
    :param untransformed_mags: exactly two raw magnitudes, by filter
               [dict, e.g., {'V': 12.6, 'I': 8.7} ]
    :param ci_filters: exactly two filters, in order, defining Color Index [e.g., ['V','I']]
    :param transforms: transforms in the ci_filters base, from which the relevant transform is
               extracted [dict, e.g., {'V': 0.025, 'I': -0.044} ]
    :return: best estimate of real, transformed Color Index [float]
    """
    mags_ordered = [untransformed_mags[f] for f in ci_filters]
    transforms_ordered = [transforms[f] for f in ci_filters]
    real_ci = (mags_ordered[0] - mags_ordered[1]) / \
              (1.0 + transforms_ordered[0] - transforms_ordered[1])
    return real_ci


def _write_aavso_report_map_stub(an_top_directory, an_rel_directory):
    lines = [";----- This is report_map.txt for AN folder " + an_rel_directory,
             ";----- Use this file to omit and/or combine target observations from AAVSO report.",
             ";----- Example directive lines:",
             ";",
             ";#TARGET  GD Cyg ; to omit this target star altogether from AAVSO report.",
             ";#JD  0.233 0.311 ; to omit this JD (fractional) range from AAVSO report.",
             ";#SERIAL  34 44,129  32  1202 ; to omit these 5 Serial numbers from AAVSO report.",
             ";#COMBINE  80,128 ; to combine (average) these 2 Serial numbers within AAVSO report.",
             ";----- Add your directive lines:",
             ";"]
    lines = [line + '\n' for line in lines]
    fullpath = os.path.join(an_top_directory, an_rel_directory, 'Photometry',
                            'report_map.txt')
    if not os.path.exists(fullpath):
        with open(fullpath, 'w') as f:
            f.writelines(lines)
        lines_written = len(lines)
    else:
        lines_written = 0
    return lines_written




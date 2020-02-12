import os
import re
from datetime import datetime, timezone, timedelta
from math import floor, ceil, cos, pi, sqrt, log

import numpy as np
import pandas as pd
from scipy.stats import trim_mean
import astropy.io.fits

from photrix.util import RaDec, ra_as_degrees, dec_as_degrees

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

TOP_DIRECTORY = 'C:/Astro/Images/Borea Photrix'
FITS_REGEX_PATTERN = '^(.+)\.(f[A-Za-z]{2,3})$'
FITS_EXTENSIONS = ['fts', 'fit', 'fits']  # allowed filename extensions
ISO_8601_FORMAT = '%Y-%m-%dT%H:%M:%S'
FWHM_PER_SIGMA = 2.0 * sqrt(2.0 * log(2))
SUBIMAGE_MARGIN = 1.5  # subimage pixels around outer annulus, for safety

# R_DISC altered 10 -> 9 Aug 16 2019 for new L-500 mount.
R_DISC = 9  # for aperture photometry, likely to be adaptive (per image) later.
R_INNER = 15  # "
R_OUTER = 20  # "


class Image:
    """
    Holds an astronomical image and apertures for photometric processing.
    Contains a FITS object, but doesn't know of its implementation and doesn't change it.
    """
    def __init__(self, fits_object):
        """
        Main constructor when FITS object is already available.
        :param fits_object: an object of the FITS class (this module).
        """
        self.fits = fits_object
        self.top_directory = fits_object.top_directory
        self.rel_directory = fits_object.rel_directory
        self.plate_solution = fits_object.plate_solution
        self.image = self.fits.image
        self.xsize = self.image.shape[0]
        self.ysize = self.image.shape[1]
        self.apertures = dict()  # initially empty dictionary of Aperture objects
        self.df_punches = pd.DataFrame()

    @classmethod
    def from_fits_path(cls, top_directory=TOP_DIRECTORY, rel_directory=None, filename=None):
        """
        Alternate constructor that starts by fetching the FITS object via its given filename.
        """
        this_fits = FITS(top_directory, rel_directory, filename)
        return Image(this_fits)

    def add_aperture(self, star_id, x0, y0):
        """
        Make one aperture from position (x,y) in image, refine its sub-pixel position by
            computing its bkgd-adjusted flux centroid, and append the aperture to self.apertures
            Will replace if aperture already exists for this starID.
        :param star_id: this aperture's name, e.g., '114_1', 'ST Tri'. Unique to this Image [string]
        :param x0: initial x position of aperture center (will be refined) [float]
        :param y0: initial y position of aperture center (will be refined) [float]
        :return: [None]
        """
        if len(self.df_punches) >= 1:
            df_ap_punches = self.df_punches.loc[self.df_punches['StarID'] == star_id, :]
        else:
            df_ap_punches = None
        self.apertures[star_id] = Aperture(self, star_id, x0, y0, df_ap_punches)
        self._recenter_aperture(star_id)

    def add_punches(self, df_punches):
        """
        Add all punches to this Image's dataframe of punches, then update all affected apertures.
        :param df_punches: new punches, columns=[StarID, dNorth, dEast] [pandas DataFrame]
        :return: [None]
        """
        # Apply punches if any; then for any apertures affected:
        self.df_punches = self.df_punches.append(df_punches)  # simple append, duplicates OK.
        if len(self.df_punches) >= 1:
            ap_names_affected = set(df_punches['StarID'])
        else:
            ap_names_affected = []

        # Replace all affected apertures (incl. updating centroid and results):
        for ap_name in ap_names_affected:
            if ap_name in self.apertures:
                ap_previous = self.apertures[ap_name]
                self.add_aperture(ap_name, ap_previous.xcenter, ap_previous.ycenter)  # replace it.
            else:
                print('>>>>> Warning: df_punch StarID \'' + ap_name +
                      '\' is not a valid aperture name in ' + self.fits.filename + '.')

    def results_from_aperture(self, star_id):
        """
        Return tuple of best positiion, fluxes etc for this aperture.
        :param star_id: which aperture [string]
        :return: Series of results [indexed pandas Series of floats]
        """
        ap = self.apertures[star_id]
        return pd.Series({'r_disc': ap.r_disc,
                          'r_inner': ap.r_inner,
                          'r_outer': ap.r_outer,
                          'n_disc_pixels': ap.n_disc_pixels,
                          'n_annulus_pixels': ap.n_annulus_pixels,
                          'annulus_flux': ap.annulus_flux,
                          'annulus_flux_sigma': ap.annulus_flux_sigma,
                          'net_flux': ap.net_flux,
                          'net_flux_sigma': ap.net_flux_sigma,
                          'x_centroid': ap.x_centroid,
                          'y_centroid': ap.y_centroid,
                          'fwhm': ap.fwhm,
                          'x1024': ap.x1024,
                          'y1024': ap.y1024,
                          'vignette': ap.vignette,
                          'sky_bias': ap.sky_bias,
                          'max_adu': ap.max_adu})

    def _recenter_aperture(self, ap_name, max_cycles=2, pixels_convergence=0.05):
        """
        For one Aperture object, reset center to previously calculated centroid, update entry in
            Image's dict of Aperture objects.
        :param ap_name:
        :param max_cycles: max number of recentering cycles [int]
        :param pixels_convergence: movement by fewer pixels than this stops the refinement. [float]
        :return: [None]
        """
        for i_cycle in range(max_cycles):
            ap = self.apertures[ap_name]
            x_previous, y_previous = ap.xcenter, ap.ycenter
            x_new, y_new = ap.x_centroid, ap.y_centroid
            distance_to_move = sqrt((x_new - x_previous) ** 2 + (y_new - y_previous) ** 2)
            if distance_to_move <= pixels_convergence:
                break
            # Here, the movement distance warrants a new Aperture object:
            self.apertures[ap_name] = ap.yield_recentered()

    def _recenter_all_apertures(self):
        for ap_name in self.apertures.keys():
            self._recenter_aperture(ap_name)


class Aperture:
    """
    Used directly only by class Image. Contains everything about one aperture.
    """
    def __init__(self, image_obj, star_id, x0, y0, df_punches=None):
        """
        :param image_obj: Image to which this Aperture applies [Image class object]
        :param star_id: name of this aperture [string]
        :param x0: initial x center, in pixels [float]
        :param y0: initial y center, in pixels [float]
        :param df_punches: one row for each punch columns=[StarID, dNorth, dEast].
            Only rows with StarID matching star_id will be applied. [pandas DataFrame]
        """
        self.image_obj = image_obj  # reference to Image class object for this aperture.
        self.image = image_obj.image  # reference to (x,y) array holding the image data.
        self.star_id = star_id
        self.xcenter = float(x0)
        self.ycenter = float(y0)
        self.df_punches = None  # default if no punch lines available for this aperture.
        if df_punches is not None:
            if len(df_punches) >= 1:
                self.df_punches = df_punches.loc[df_punches['StarID'] == star_id, :]
        self.r_disc = R_DISC
        self.r_inner = R_INNER
        self.r_outer = R_OUTER

        # Aperture evaluation fields, with default (no-flux) values:
        self.n_disc_pixels, self.n_annulus_pixels = 0, 0
        self.net_flux = 0.0
        self.net_flux_sigma = 0.0
        self.annulus_flux = 0.0
        self.annulus_flux_sigma = 0.0
        self.sn = 0.0
        self.x_centroid = self.xcenter
        self.y_centroid = self.ycenter
        self.fwhm = 0.0
        self.sky_bias = 0.0
        self.max_adu = 0.0

        # Compute needed boundaries of subimage around this aperture:
        image_xsize, image_ysize = self.image.shape
        test_radius = self.r_outer + SUBIMAGE_MARGIN
        xlow = int(floor(self.xcenter - test_radius))
        xhigh = int(ceil(self.xcenter + test_radius))
        ylow = int(floor(self.ycenter - test_radius))
        yhigh = int(ceil(self.ycenter + test_radius))

        # Compute whether needed subimage will fall entirely within image, or not:
        subimage_within_image = (xlow >= 0) & (xhigh <= image_xsize - 1) & \
                                (ylow >= 0) & (yhigh <= image_ysize - 1)

        # Compute values only if subimage entirely contained in current image:
        if subimage_within_image:
            self.subimage = self.image[xlow:xhigh + 1, ylow:yhigh + 1].copy()

            # Construct mask arrays to represent disc and annulus (both same shape as subimage):
            nx = xhigh - xlow + 1  # number of columns in subimage.
            ny = yhigh - ylow + 1  # number of rows in subimage.
            self.ygrid, self.xgrid = np.meshgrid(ylow + np.arange(ny), xlow + np.arange(nx))
            dx = self.xgrid - self.xcenter
            dy = self.ygrid - self.ycenter
            dist2 = dx**2 + dy**2
            self.disc_mask = np.clip(np.sign(self.r_disc**2 - dist2), 0.0, 1.0)
            inside_outer_edge = np.clip(np.sign(self.r_outer**2 - dist2), 0.0, 1.0)
            outside_inner_edge = np.clip(np.sign(dist2 - self.r_inner**2), 0.0, 1.0)
            self.annulus_mask = inside_outer_edge * outside_inner_edge

            # Apply punches:
            if df_punches is not None:
                if len(df_punches) >= 1:
                    self._apply_punches(image_obj.plate_solution)  # only punches for this aperture.

            # Evaluate and store several new fields:
            self.evaluate()
            del self.subimage, self.xgrid, self.ygrid, self.disc_mask, self.annulus_mask

        # Add other fields useful to calling code:
        image_center_x = self.image.shape[0] / 2.0
        image_center_y = self.image.shape[1] / 2.0
        self.x1024 = (self.xcenter - image_center_x) / 1024.0
        self.y1024 = (self.ycenter - image_center_y) / 1024.0
        self.vignette = self.x1024**2 + self.y1024**2  # no sqrt...meant to be parabolic term

    def evaluate(self):
        """
        Compute and several fields in this Aperture object. Put them in Aperture object.
        :return: [None]
        """
        self.n_disc_pixels = np.sum(self.disc_mask)
        self.n_annulus_pixels = np.sum(self.annulus_mask)
        self.annulus_flux = self._eval_sky_005()  # average adus / pixel, sky background
        estimated_background = self.n_disc_pixels * self.annulus_flux
        disc_values = np.ravel(self.subimage[self.disc_mask > 0])  # only values in mask.
        self.max_adu = np.max(disc_values)
        this_net_flux = np.sum(disc_values) - estimated_background
        if this_net_flux > 0:
            self.net_flux = this_net_flux
            gain = 1.57  # TODO: this should come from Instrument object.
            annulus_values = np.ravel(self.subimage[self.annulus_mask > 0])
            self.annulus_flux_sigma = np.std(annulus_values)
            # net_flux_sigma equation after APT paper, PASP 124, 737 (2012), but pi/2 in 3rd term
            # set to 1 as pi/2 seems hard to justify, and as 1 gives S/N closer to VPhot's values.
            self.net_flux_sigma = sqrt((self.net_flux / gain) +
                                       (self.n_disc_pixels * self.annulus_flux_sigma**2) +
                                       1.0 * ((self.n_disc_pixels*self.annulus_flux_sigma)**2 /
                                              self.n_annulus_pixels))
            self.sn = self.net_flux / self.net_flux_sigma
            # Compute centroid (x,y) of net flux:
            net_flux_grid = self.disc_mask * (self.subimage - self.annulus_flux)
            normalizor = np.sum(net_flux_grid)
            if (self.x_centroid is not None) and (self.y_centroid is not None):
                self.xcenter = self.x_centroid  # new subimage center
                self.ycenter = self.y_centroid  # "
            self.x_centroid = np.sum(net_flux_grid * self.xgrid) / normalizor
            self.y_centroid = np.sum(net_flux_grid * self.ygrid) / normalizor

            # Other evaluation results:
            self.fwhm = self._eval_fwhm()
            sky_flux_bias = self.n_disc_pixels * self.annulus_flux_sigma
            self.sky_bias = abs(-2.5 * (sky_flux_bias / self.net_flux) / log(10.0))

    def yield_recentered(self):
        x_new, y_new = self.x_centroid, self.y_centroid
        return Aperture(self.image_obj, self.star_id, x_new, y_new, self.df_punches)

    def _apply_punches(self, plate_solution):
        """
        Apply punches to (remove appropriate pixels from) this Aperture's annulus mask.
        :param plate_solution:
        :return: [None]
        """
        dnorth_dx = 3600.0 * plate_solution['CD2_1']  # in arcseconds northward /pixel
        dnorth_dy = 3600.0 * plate_solution['CD2_2']  # "
        deast_dx = 3600.0 * plate_solution['CD1_1']   # in arcseconds eastward (not RA) /pixel
        deast_dy = 3600.0 * plate_solution['CD1_2']   # "
        ann_mask_new = self.annulus_mask.copy()  # to begin.
        for dnorth, deast in zip(self.df_punches['dNorth'], self.df_punches['dEast']):
            coefficients = np.array([[dnorth_dx, dnorth_dy], [deast_dx, deast_dy]])
            dep_vars = np.array([dnorth, deast])
            solution = np.linalg.solve(coefficients, dep_vars)
            dx_punch, dy_punch = solution[0], solution[1]
            x_punch = self.xcenter + dx_punch
            y_punch = self.ycenter + dy_punch
            x_dist = self.xgrid - x_punch  # x distance from center
            y_dist = self.ygrid - y_punch  # y distance from center
            dist2 = x_dist**2 + y_dist**2
            punch_mask = np.clip(np.sign(dist2 - self.r_disc**2), 0.0, 1.0)
            ann_mask_new = ann_mask_new * punch_mask  # do the punch (pixels within punch set to 0).
        self.annulus_mask = ann_mask_new

    def _eval_sky_005(self):
        """
        Winning sky-background measurement strategy of 2015 tournament of strategies.
        Insensitive to cosmic rays and background stars in or near the annulus.
        :return robust estimate of sky background in adu/pixel [float]
        """
        slice_list = self._make_sky_slices(n_slices=12, method='trimmed_mean')
        sky_adu = trim_mean(slice_list, proportiontocut=0.3)
        return sky_adu

    def _make_sky_slices(self, n_slices=12, method='trimmed_mean'):
        radians_per_slice = (2.0 * pi) / n_slices
        min_pixels_per_slice = 0.5 * (self.n_annulus_pixels / n_slices)
        angle_grid = np.arctan2(self.ygrid-self.ycenter, self.xgrid-self.xcenter)  # -pi to +pi
        slice_list = []
        for i_slice in range(n_slices):
            # Radians delimiting this slice:
            angle_min = i_slice * radians_per_slice - pi
            angle_max = (i_slice + 1) * radians_per_slice - pi
            above_min = np.clip(np.sign(angle_grid - angle_min), 0.0, 1.0)
            below_max = np.clip(np.sign(angle_max - angle_grid), 0.0, 1.0)
            slice_mask = above_min * below_max * self.annulus_mask
            n_slice_pixels = np.sum(slice_mask)
            if n_slice_pixels >= min_pixels_per_slice:
                slice_values = np.ravel(self.subimage[slice_mask > 0])  # only values in mask.
                slice_mean = trim_mean(slice_values, 0.4)
                slice_list.append(slice_mean)
        return slice_list

    def _eval_fwhm(self):
        # TODO: Probably need a better FWHM algorithm.
        """
        Returns estimate of Full Width at Half-Maximum from mean dist2 (=2*sigma^2) of net flux.
        This algorithm may be replaced later:
            overestimates FWHM compared to MaxIm, PinPoint, and sometimes to even visual inspection.
        :return: estimate of FWHM in pixels.
        """
        dx = self.xgrid - self.x_centroid
        dy = self.ygrid - self.y_centroid
        dist2 = dx ** 2 + dy ** 2
        net_flux_xy = (self.disc_mask * (self.subimage - self.annulus_flux))
        mean_dist2 = max(0.0, np.sum(net_flux_xy * dist2) / np.sum(net_flux_xy))
        sigma = sqrt(mean_dist2 / 2.0)
        # this math is verified 20170723, but yields larger FWHM than does MaxIm.
        return FWHM_PER_SIGMA * sigma


class FITS:
    """
    Holds data from a FITS file. Immutable. Used mostly by an Image object (class Image).
    Isolates details of FITS implementation from calling code.
    """
    def __init__(self, top_directory, rel_directory, filename):
        """
        :param top_directory:
        :param rel_directory:
        :param filename:
        """
        # If filename has FITS extension already, use it:
        actual_fits_fullpath = None
        test_fullpath = os.path.join(top_directory, rel_directory, filename)
        if os.path.exists(test_fullpath):
            actual_fits_fullpath = test_fullpath

        # If no file with FITS extension, try to find a matching FITS filename (w/extension):
        if actual_fits_fullpath is None:
            for fits_ext in FITS_EXTENSIONS:
                test_fullpath = os.path.join(top_directory, rel_directory,
                                             filename + '.' + fits_ext)
                if os.path.exists(test_fullpath):
                    actual_fits_fullpath = test_fullpath
                break

        if actual_fits_fullpath is None:
            print("Not a valid file name: '" + filename + "'")
            self.is_valid = False
            return

        self.fullpath = actual_fits_fullpath
        try:
            hdulist = astropy.io.fits.open(self.fullpath)
        except IOError:
            self.is_valid = False
            return

        self.header = hdulist[0].header
        # FITS convention = (vert/Y, horiz/X), pixel (1,1) at bottom left -- NOT USED by photrix.
        # MaxIm/Astrometrica convention = (horiz/X, vert/Y) pixel (0,0 at top left). USE THIS.
        # NB: self.image_fits, self.image_xy, and self.image are different views of the SAME array.
        #     They are meant to be read-only--changing any one of them *will* change the others.
        self.image_fits = hdulist[0].data.astype(np.float64)
        self.image_xy = np.transpose(self.image_fits)  # x and y axes as expected (not like FITS).
        self.image = self.image_xy  # alias
        hdulist.close()

        self.top_directory = top_directory
        self.rel_directory = rel_directory
        self.filename = filename
        self.all_header_keys = self.header.keys()
        self.object = self.header_value('OBJECT')
        self.is_calibrated = self._is_calibrated()
        self.focal_length = self._get_focal_length()
        self.exposure = self.header_value(['EXPTIME', 'EXPOSURE'])  # seconds
        self.temperature = self.header_value(['SET-TEMP', 'CCD-TEMP'])  # deg C
        self.utc_start = self._get_utc_start()
        self.utc_mid = self.utc_start + timedelta(seconds=self.exposure / 2.0)
        self.filter = self.header_value('FILTER')
        self.airmass = self.header_value('AIRMASS')
        self.guide_exposure = self.header_value('TRAKTIME')  # seconds
        self.fwhm = self.header_value('FWHM')  # pixels

        self.plate_solution = self._get_plate_solution()  # a pd.Series
        self.is_plate_solved = not any(self.plate_solution.isnull())
        self.ra = ra_as_degrees(self.header_value(['RA', 'OBJCTRA']))
        self.dec = dec_as_degrees(self.header_value(['DEC', 'OBJCTDEC']))

        self.is_valid = True  # if it got through all that initialization.
        # self.is_valid = all(x is not None
        #                     for x in [self.object, self.exposure, self.filter,
        #                               self.airmass, self.utc_start, self.focal_length])

    def header_value(self, key):
        """
        :param key: FITS header key [string] or list of keys to try [list of strings]
        :return: value of FITS header entry, typically [float] if possible, else [string]
        """
        if isinstance(key, str):
            return self.header.get(key, None)
        for k in key:
            value = self.header.get(k, None)
            if value is not None:
                return value
        return None

    def header_has_key(self, key):
        return key in self.header

    def xy_from_radec(self, radec):
        """
        Computes zero-based pixel x and y for a given RA and Dec sky coordinate.
            May be outside image's actual boundaries.
            Assumes flat image (no distortion, i.e., pure Tan projection).
        :param radec: sky coordinates [RaDec class object]
        :return: x and y pixel position, zero-based, in this FITS image [2-tuple of floats]
        """
        cd11 = self.plate_solution['CD1_1']
        cd12 = self.plate_solution['CD1_2']
        cd21 = self.plate_solution['CD2_1']
        cd22 = self.plate_solution['CD2_2']
        crval1 = self.plate_solution['CRVAL1']
        crval2 = self.plate_solution['CRVAL2']
        crpix1 = self.plate_solution['CRPIX1']  # 1 at edge (FITS convention)
        crpix2 = self.plate_solution['CRPIX2']  # "

        d_ra = radec.ra - crval1
        d_dec = radec.dec - crval2
        deg_ew = d_ra * cos((pi / 180.0) * radec.dec)
        deg_ns = d_dec
        a = cd22 / cd12
        dx = (deg_ns - deg_ew * a) / (cd21 - cd11 * a)
        dy = (deg_ew - cd11 * dx) / cd12
        x = crpix1 + dx
        y = crpix2 + dy
        return x - 1, y - 1  # FITS image origin=(1,1), but our (MaxIm/python) convention=(0,0)

    def radec_from_xy(self, x, y):
        pass
    #     """
    #     Computes RA and Dec for a give x and y pixel count. Assumes flat image (no distortion,
    #         i.e., pure Tan projection).
    #     :param x: pixel position in x [float]
    #     :param y: pixel position in y [float]
    #     :return: RA and Dec [RaDec object]
    #     """
    #     cd11 = self.plate_solution['CD1_1']
    #     cd12 = self.plate_solution['CD1_2']
    #     cd21 = self.plate_solution['CD2_1']
    #     cd22 = self.plate_solution['CD2_2']
    #     crval1 = self.plate_solution['CRVAL1']
    #     crval2 = self.plate_solution['CRVAL2']
    #     crpix1 = self.plate_solution['CRPIX1']
    #     crpix2 = self.plate_solution['CRPIX2']
    #     # Do the calculation (inverse of self.xy_from_radec(self, radec)):
    #     return RaDec(0, 0)

    def _is_calibrated(self):
        calib_fn_list = [self._is_calibrated_by_maxim_5_6()]  # may add more fns when available.
        return any([is_c for is_c in calib_fn_list])

    def _is_calibrated_by_maxim_5_6(self):
        hval = self.header_value('CALSTAT')
        if hval is not None:
            if hval.strip().upper() == 'BDF':  # calib. by MaxIm DL v. 5 or 6
                return True
        return False

    def _get_focal_length(self):
        # If FL available, return it. Else, compute FL from plate solution.
        value = self.header_value('FOCALLEN')
        if value is not None:
            return value  # mm
        x_pixel = self.header_value('XPIXSZ')
        y_pixel = self.header_value('YPIXSZ')
        x_scale = self.header_value('CDELT1')
        y_scale = self.header_value('CDELT2')
        if any([val is None for val in [x_pixel, y_pixel, x_scale, y_scale]]):
            return None
        fl_x = x_pixel / abs(x_scale) * (206265.0 / (3600 * 1800))
        fl_y = y_pixel / abs(y_scale) * (206265.0 / (3600 * 1800))
        return (fl_x + fl_y) / 2.0

    def _get_utc_start(self):
        utc_string = self.header_value('DATE-OBS')
        utc_dt = datetime.strptime(utc_string, ISO_8601_FORMAT).replace(tzinfo=timezone.utc)
        return utc_dt

    def _get_plate_solution(self):
        plate_solution_index = ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                                'CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2']
        plate_solution_values = [np.float64(self.header_value(key))
                                 for key in plate_solution_index]
        return pd.Series(plate_solution_values, index=plate_solution_index)


def all_fits_files(top_directory, rel_directory, validate_fits=False):
    """
    Return list of all FITS files in given directory.
    :param top_directory:
    :param rel_directory:
    :param validate_fits: If True, open FITS files and include only if valid.
        If False, include filename if it appears valid without opening the FITS file.
    :return: List of all FITS files in given directory [list of strings]
    """
    pass
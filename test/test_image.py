import os

import numpy as np
import pandas as pd
from datetime import datetime, timezone, timedelta
import pytest

from photrix import image

__author__ = "Eric Dose :: New Mexico Mira Project, Albuquerque"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_TOP_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test")


def test_class_fits():
    test_rel_directory = '$data_for_test'

    # Test failure on missing file:
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory=test_rel_directory, filename='$no_file.txt')
    assert fits.is_valid is False

    # Test exception handling for non-FITS file format:
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory=test_rel_directory, filename='dummy.txt')
    assert fits.is_valid is False

    # Open FITS file with known extension:
    given_filename = 'CE Aur-0001-V.fts'
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory=test_rel_directory,
                      filename=given_filename)
    assert fits.is_valid
    assert fits.fullpath == os.path.join(TEST_TOP_DIRECTORY, test_rel_directory, given_filename)
    assert fits.header_has_key('NAXIS')
    assert not fits.header_has_key('XXX')
    assert fits.object == 'CE Aur'
    assert fits.is_calibrated
    assert fits.is_plate_solved
    assert fits.focal_length == pytest.approx(2702, abs=1)
    assert fits.exposure == pytest.approx(587, abs=1)
    assert fits.temperature == pytest.approx(-35, abs=0.1)
    target_start_utc = datetime(2017, 4, 24, 4, 0, 31).replace(tzinfo=timezone.utc)
    diff_seconds = (fits.utc_start - target_start_utc).total_seconds()
    assert abs(diff_seconds) < 1
    target_mid_utc = target_start_utc + timedelta(seconds=fits.exposure / 2.0)
    diff_seconds = (fits.utc_mid - target_mid_utc).total_seconds()
    assert abs(diff_seconds) < 1
    assert fits.filter == 'V'
    assert fits.airmass == pytest.approx(1.5263, abs=0.0001)
    assert fits.guide_exposure == pytest.approx(1.4, abs=0.001)
    assert fits.fwhm == pytest.approx(5.01, abs=0.01)
    assert fits.plate_solution['CD1_1'] == pytest.approx(-1.92303985969E-006)
    assert fits.plate_solution['CD2_1'] == pytest.approx(1.90588522664E-004)
    assert fits.plate_solution['CRVAL1'] == pytest.approx(1.03834010522E+002)
    assert fits.ra == pytest.approx(103.83791666666667)
    assert fits.dec == pytest.approx(46.28638888888889)

    # Test .image (== .image_xy):
    assert fits.image.shape == (3072, 2047)  # *image* (x,y), which is *array* (n_rows, n_columns)
    assert fits.image[0, 0] == 275  # upper-left corner
    assert fits.image[0, 2046] == 180  # lower-left corner
    assert fits.image[3071, 2046] == 265  # lower-right corner
    assert fits.image[3071, 0] == 285  # upper-right corner

    # Open FITS file without known FITS extension (which FITS constructor must first find):
    given_filename = 'CE Aur-0001-V'
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory='$data_for_test',
                      filename=given_filename)
    assert fits.is_valid
    assert fits.fullpath == os.path.join(TEST_TOP_DIRECTORY, test_rel_directory,
                                         given_filename +'.fts')
    assert fits.object == 'CE Aur'
    assert fits.airmass == pytest.approx(1.5263, abs=0.0001)

    # Open FITS file with no calibration (at least not by MaxIm 5/6) and no plate solution:
    given_filename = 'AD Dra-S001-R001-C001-I.fts'
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory='$data_for_test', filename=given_filename)
    assert fits.is_valid
    assert fits.is_calibrated is False
    assert fits.is_plate_solved is False


def test_class_image():
    test_rel_directory = '$data_for_test'

    # Open FITS file with known extension:
    given_filename = 'CE Aur-0001-V.fts'
    fits_obj = image.FITS(TEST_TOP_DIRECTORY, rel_directory=test_rel_directory,
                          filename=given_filename)
    im = image.Image(fits_obj)
    assert im.fits.object == 'CE Aur'
    assert im.top_directory == TEST_TOP_DIRECTORY
    assert im.rel_directory == test_rel_directory
    assert im.xsize, im.ysize == fits_obj.image_xy.shape  # .shape is in nrows, ncols
    # Image dimensions are x,y == *image* cols,rows, the reverse of numpy storage.
    # Images are zero based; [0,0] -> upper-left, [n, 0] is on top edge of *image* (not of storage).

    assert im.image.shape == (3072, 2047)
    assert im.image[0, 0] == 275        # upper-left corner
    assert im.image[0, 2046] == 180     # lower-left corner
    assert im.image[3071, 2046] == 265  # lower-right corner
    assert im.image[3071, 0] == 285     # upper-right corner

    # Aperture very simple case: near image center, no punches or interfering signals:
    im.add_aperture('dummy_1', 1523, 1011)  # star near image center, no punches.
    assert len(im.apertures) == 1
    this_ap = im.apertures['dummy_1']
    assert this_ap.x_centroid == pytest.approx(1524.784, abs=0.005)
    results = im.results_from_aperture('dummy_1')
    assert results['x_centroid'] == this_ap.x_centroid
    assert results['fwhm'] == pytest.approx(6.42, abs=0.02)
    assert set(results.index) == set(['r_disc', 'r_inner', 'r_outer', 'n_disc_pixels',
                                      'n_annulus_pixels', 'net_flux', 'net_flux_sigma',
                                      'annulus_flux', 'annulus_flux_sigma',
                                      'x_centroid', 'y_centroid', 'fwhm',
                                      'x1024', 'y1024', 'vignette'])

    # Aperture case: near image center, two punches:
    im.add_aperture('dummy_2', 1535, 979)
    df_punches = pd.DataFrame({'StarID': 'dummy_2',
                               'dNorth': [-11.1, +9.6],
                               'dEast': [0.0, +3.4]})
    im.add_punches(df_punches)
    assert len(im.apertures) == 2
    this_ap = im.apertures['dummy_2']
    assert [this_ap.x_centroid, this_ap.y_centroid] == pytest.approx([1534.456, 978.697], abs=0.005)
    results = im.results_from_aperture('dummy_2')
    assert results['x_centroid'] == this_ap.x_centroid
    assert results['fwhm'] == pytest.approx(6.15, abs=0.02)

    # Aperture case: far from image center, one punch:
    im.add_aperture('dummy_3', 510, 483)
    df_punches = pd.DataFrame({'StarID': ['dummy_2', 'trash', 'dummy_3'],
                               'dNorth': [-11.1, -99.9, +8.9],
                               'dEast': [0.0, +99.9, 0.0]})  # verify safety of non-relevant rows.
    im.add_punches(df_punches)
    assert len(im.apertures) == 3
    this_ap = im.apertures['dummy_3']
    assert [this_ap.x_centroid, this_ap.y_centroid] == pytest.approx([505.53, 481.35], abs=0.005)
    results = im.results_from_aperture('dummy_3')
    assert results['y_centroid'] == this_ap.y_centroid
    assert results['fwhm'] == pytest.approx(7.05, abs=0.05)
    result_3 = im.results_from_aperture('dummy_3')
    ap_3 = im.apertures['dummy_3']
    assert result_3['x_centroid'] == ap_3.x_centroid


def test_fits__xy_from_radec():
    from photrix.util import RaDec
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory='$data_for_test',
                      filename='CE Aur-0001-V.fts')
    # All tests lack distortion corrections (as none available in FITS header),
    #    and so in real images calculated (x,y) values at edges will not quite line up with stars.
    radec_near_center = RaDec('06:55:21.25', '+46:17:33.0')
    x, y = fits.xy_from_radec(radec_near_center)
    assert list((x, y)) == pytest.approx([1557.6, 1005.8], abs=0.25)

    radec_upper_left = RaDec('06:56:10.6', '+46:02:27.1')
    x, y = fits.xy_from_radec(radec_upper_left)
    assert list((x, y)) == pytest.approx([229.8, 270.3], abs=0.25)

    radec_upper_right = RaDec('06:56:14.3', '+46:29:11.6')
    x, y = fits.xy_from_radec(radec_upper_right)
    assert list((x, y)) == pytest.approx([2567.6, 197.2], abs=0.25)

    radec_lower_left = RaDec('06:54:26.0', '+46:02:44.9')
    x, y = fits.xy_from_radec(radec_lower_left)
    assert list((x, y)) == pytest.approx([271.8, 1857.0], abs=0.25)

    radec_lower_right = RaDec('06:54:18.0', '+46:30:02.0')
    x, y = fits.xy_from_radec(radec_lower_right)
    assert list((x, y)) == pytest.approx([2658.7, 1946.6], abs=0.25)



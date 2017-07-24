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
    assert fits.image.shape == (2047, 3072)
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
    assert results['fwhm'] == pytest.approx(6.42, abs=0.01)
    assert set(results.index) == set(['r_disc', 'r_inner', 'r_outer', 'n_disc_pixels',
                                      'n_annulus_pixels', 'net_flux', 'net_flux_sigma',
                                      'x_centroid', 'y_centroid', 'fwhm'])

    # Aperture case: near image center, two punches:
    im.add_aperture('dummy_2', 1535, 979)
    df_punches = pd.DataFrame({'StarID': 'dummy_2',
                               'dNorth': [-10.0, 0.0],
                               'dEast': [+9.0, +3.4]})
    im.add_punches(df_punches)
    assert len(im.apertures) == 2




# def test_fits__xy_from_radec():
#     from photrix.util import RaDec
#     fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory='$data_for_test',
#                       filename='CE Aur-0001-V.fts')
#     radec1 = RaDec('06:55:58.22', '+46:14:23.4')
#     x1, y1 = fits.xy_from_radec(radec1)
#     assert list((x1, y1)) == pytest.approx([1276.6, 448.9], abs=0.25)



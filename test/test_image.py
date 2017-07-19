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

    # Open FITS file without known FITS extension (which FITS constructor must first find):
    given_filename = 'CE Aur-0001-V'
    fits = image.FITS(TEST_TOP_DIRECTORY, rel_directory='$data_for_test',
                      filename=given_filename)
    assert fits.is_valid
    assert fits.fullpath == os.path.join(TEST_TOP_DIRECTORY, test_rel_directory,
                                         given_filename +'.fts')
    assert fits.object == 'CE Aur'
    assert fits.airmass == pytest.approx(1.5263, abs=0.0001)





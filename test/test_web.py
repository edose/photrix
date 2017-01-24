import os
import pytest   # so we can have the @pytest... decorator
import ephem

from photrix.web import *
from photrix.user import Astronight

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTRIX_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEST_DATA_DIRECTORY = os.path.join(PHOTRIX_ROOT_DIRECTORY, "test", "$data_for_test")


@pytest.mark.webtest
def test_get_aavso_webobs_raw_table():  # ========================================================
    # Valid star case (for reasonableness only, as webobs data may become updated at any time).
    star_name = "ST Tri"
    num_obs = 60
    df = get_aavso_webobs_raw_table(star_name, num_obs)
    assert 25 <= len(df) <= MAX_WEBOBS_LINES   # webobs API always returns at least 25 lines.
    assert len(df[df['filter'] == "V"]) >= 1
    required_cols = ['target_name', 'jd', 'date_string', 'mag', 'error', 'filter', 'observer']
    assert all([colname in list(df.columns) for colname in required_cols])
    assert all(name.lower() == "st tri" for name in df['target_name'])
    assert all(10 < mag < 20 for mag in df['mag'])
    assert all(0 <= err < 1 for err in df['error'])
    assert 2457372 < df['jd'].max() < ephem.julian_date()  # reasonable JD values.
    assert df['jd'].min() > 2457372  # day before my own known obs date
    assert all([df['jd'].iloc[i] >= df['jd'].iloc[i+1] for i in range(len(df)-1)])  # sorted JDs.
    assert all(observer.isalpha() and (2 <= len(observer) <= 4) for observer in df['observer'])

    # Test invalid star (no data in AAVSO webobs):
    star_name = "Not a Star"
    df = get_aavso_webobs_raw_table(star_name)
    assert len(df) == 0   # empty dataframe
    required_cols = ['target_name', 'jd', 'date_string', 'mag', 'error', 'filter', 'observer']
    assert all([colname in list(df.columns) for colname in required_cols])  # but w/required cols.

    # Test valid star with start and end Julian Dates (but enough to cover num_results):
    star_name = "ST Tri"
    num_obs = 150
    jd_start = 2455190.87  # webobs has 192 obs between these two Julian Dates (1/20/2017).
    jd_end = 2457400
    df = get_aavso_webobs_raw_table(star_name, num_obs, jd_start, jd_end)
    assert len(df) == 150
    assert jd_start <= min(df['jd']) < max(df['jd']) < jd_end
    assert all([df['jd'].iloc[i] >= df['jd'].iloc[i+1] for i in range(len(df)-1)])  # sorted JDs.
    assert set(df['observer']) == set(['DERA'])
    data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_150.csv")
    df150 = pd.read_csv(data_fullpath, index_col=0)
    assert len(df) == len(df150)
    assert all([df['date_string'].iloc[irow] == df150['date_string'].iloc[irow]
               for irow in range(len(df))])



    # data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_150.csv")
    # df.to_csv(data_fullpath)

    # Test valid star with start and end Julian Dates (but NOT enough to cover num_results):
    star_name = "ST Tri"
    num_obs = 200
    jd_start = 2455190.87  # webobs has 192 obs between these two Julian Dates (1/20/2017).
    jd_end = 2457400
    df = get_aavso_webobs_raw_table(star_name, num_obs, jd_start, jd_end)
    assert len(df) == 192
    assert jd_start <= min(df['jd']) < max(df['jd']) < jd_end
    assert all([df['jd'].iloc[i] >= df['jd'].iloc[i+1] for i in range(len(df)-1)])  # sorted JDs.
    assert set(df['observer']) == set(['DERA', 'SX'])
    data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_192.csv")
    df192 = pd.read_csv(data_fullpath, index_col=0)
    assert len(df) == len(df192)
    assert all([df['date_string'].iloc[irow] == df192['date_string'].iloc[irow]
               for irow in range(len(df))])
    # data_fullpath = os.path.join(TEST_DATA_DIRECTORY, "ST_Tri_192.csv")
    # df.to_csv(data_fullpath)

    # # Test .most_recent_jd_mag():
    # filter_code = "V"
    # filter_bools = df['filter'].str.lower() == filter_code.lower()
    # df2 = df[filter_bools]
    # jd_target = df2['jd'].max()
    # jd_bools = df2['jd'] == jd_target
    # df3 = df2[jd_bools]
    # assert obj.most_recent_jd_mag(filter='V')[0] == jd_target
    # assert obj.most_recent_jd_mag(filter='V')[1] == df3['mag'].iloc[0]

    # # Test .days_gap_jd():
    # an_date_string = "20160910"
    # site_string = "BDO_Kansas"
    # an = Astronight(an_date_string, site_string)
    # days_target = an.local_middark_jd - jd_target  # from prev test section.
    # assert obj.days_gap_jd(an.local_middark_jd, filter='V') == days_target
    #
    # # Test .days_gap_an():
    # jd_an = an.local_middark_jd
    # assert obj.days_gap_an(an, filter='V') == jd_an - obj.most_recent_jd_mag(filter='V')[0]






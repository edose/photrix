import pytest   # so we can have the @pytest... decorator
import ephem
from photrix.web import *
from photrix.user import Astronight

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


@pytest.mark.webtest
def test_web_obs_aavso():  # =====================================================================
    star_name = "ST Tri"
    obj = WebObsAAVSO(star_name)
    df = obj.table

    # Test class constructor:
    assert len(df) == 50
    assert len(df[df['filter'] == "V"]) >= 1
    required_cols = ['target_name', 'jd', 'date_string', 'mag', 'error', 'filter', 'observer']
    assert all([colname in list(df.columns) for colname in required_cols])
    assert all(name.lower() == "st tri" for name in df['target_name'])
    assert all(10 < mag < 20 for mag in df['mag'])
    assert all(0 <= err < 1 for err in df['error'])
    assert ephem.julian_date() > df['jd'].max() > 2457372  # my own obs date
    assert df['jd'].min() > 2457372  # day before my own obs date
    assert all(observer.isalpha() and 2 <= len(observer) <= 4 for observer in df['observer'])

    # Test .most_recent_jd_mag():
    filter_code = "V"
    filter_bools = df['filter'].str.lower() == filter_code.lower()
    df2 = df[filter_bools]
    jd_target = df2['jd'].max()
    jd_bools = df2['jd'] == jd_target
    df3 = df2[jd_bools]
    assert obj.most_recent_jd_mag(filter='V')[0] == jd_target
    assert obj.most_recent_jd_mag(filter='V')[1] == df3['mag'].iloc[0]

    # Test .days_gap_jd():
    an_date_string = "20160910"
    site_string = "BDO_Kansas"
    an = Astronight(an_date_string, site_string)
    days_target = an.local_middark_jd - jd_target  # from prev test section.
    assert obj.days_gap_jd(an.local_middark_jd, filter='V') == days_target

    # Test .days_gap_an():
    jd_an = an.local_middark_jd
    assert obj.days_gap_an(an, filter='V') == jd_an - obj.most_recent_jd_mag(filter='V')[0]






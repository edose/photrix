from web_aavso import *
import pytest   # so we can have the @pytest... decorator

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


@pytest.mark.webtest
def test_web_obs_aavso():
    w = WebObsAAVSO("ST Tri")
    assert len(w.observations) >= 1

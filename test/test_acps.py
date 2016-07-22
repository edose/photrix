from photrix import acps  # don't include "acps." in calls to functions and classes.

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


def test_ACPS_observation():
    obs = acps.ACPS_observation("ST Tri", 34.555, 21.334)
    obs.add_imageset("ST Tri", 3, 120, 'V')
    obs.add_imageset("ST Tri I filter", 3, 30, 'I')
    txt = obs.rtml()
    print("\n" + txt + "\n")


def test_run():
    project = acps.ACPS_project("AN20160630-BOREA")  # for the whole night on one instrument.
    plan = project.make_plan("first plan")
    plan.horizon = 30  # override class defaults
    plan.priority = 4  # "
    obs = acps.ACPS_observation('obs_id', 37, -1)
    obs.add_imageset("", 5, 60, 'V')
    obs.add_imageset("name2", 2, 120, 'I')
    plan.add_observation(obs)
    project.add_plan(plan)
    print('\n'*2 + project.rtml())
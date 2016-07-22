from photrix.user import Instrument

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"


# package photrix.ACPS: generate RTML for ACP Scheduler
# Typical usage:
#
# from photrix.acps import *
# project = ACPS_project("AN20160630-BOREA")  # for the whole night on one instrument.
# for plan_source in plan_sources:
#     plan = project.make_plan(id)
#     plan.horizon = 30  # override class defaults (optional)
#     plan.priority = 4  #  "
#     obs = ACPS_observation(id, RA_deg, dec_deg)
#     obs.add_imageset(name1, count1, exposure1, filter1)
#     obs.add_imageset(name2, count2, exposure2, filter2)    ...
#     plan.add_observation(obs)
#     project.add_plan(plan)
# rtml_text = project.rtml()

class ACPS_project:
    """

    """
    def __init__(self, project_name=None, instrument_name="Borea",
                 user="Eric", email="astro@ericdose.com",
                 organization="Bois d'Arc Observatory, Kansas"):
        if project_name is None:
            raise ValueError('project name may not be null.')
        if len(project_name) <= 0:
            raise ValueError('project name may not be zero-length.')
        self.project_name = project_name
        self.instrument = Instrument(instrument_name)
        self.telescope = instrument_name
        self.user = user
        self.email = email
        self.organization = organization
        self.plan_list = []

    def make_plan(self, plan_id):
        plan = ACPS_plan(plan_id, self.project_name, self.telescope, self.user)
        plan.horizon = self.instrument.min_altitude
        plan.moon_distance = self.instrument.min_distance_full_moon
        plan.min_exposure = self.instrument.camera["shortest_exposure"]
        return plan

    def add_plan(self, new_acps_plan):
        if isinstance(new_acps_plan, ACPS_plan):
            self.plan_list.append(new_acps_plan)

    def rtml(self):
        rtml_text = '<?xml version="1.0"?>\n' + \
            '<RTML xmlns:xs="http://www.w3.org/2001/XMLSchema-instance" version="2.3">\n' + \
                    xml_open(1, "Contact") + \
                    xml_line(2, "User", self.user) + \
                    xml_line(2, "Email", self.email) + \
                    xml_line(2, "Organization", self.organization) + \
                    xml_close(1, "Contact")
        for plan in self.plan_list:
            rtml_text += plan.rtml()
        rtml_text += xml_close(0, "RTML")
        return rtml_text


class ACPS_plan:
    """
    ACPS plan is

    """
    def __init__(self, plan_id, project_name, telescope, user):
        if plan_id is None:
            raise ValueError('plan id may not be null.')
        if len(plan_id) <= 0:
            raise ValueError('plan id may not be zero-length.')
        self.plan_id = plan_id
        self.user_name = user
        self.description = ""
        self.reason = ""
        self.project_name = project_name
        self.airmass = None                # default
        self.airmass_range_minimum = None  # default
        self.airmass_range_maximum = None  # default
        self.horizon = 28                  # default, in degrees
        self.hour_angle_range_east = None  # default
        self.hour_angle_range_west = None  # default
        self.sky_condition = "Fair"        # default
        self.moon_distance = 50            # default, in degrees
        self.moon_width = 10               # default, in days
        self.earliest = None
        self.latest = None
        self.priority = 4            # default; NB: mean for all plans post-normalized to 0.5
        self.telescope = telescope
        self.min_exposure = 0              # default
        self.observation_list = []

    def add_observation(self, ACPS_observation_to_add):
        self.observation_list.append(ACPS_observation_to_add)

    def rtml(self):
        rtml_text = xml_open(1, "Request") + \
            xml_line(2, "ID", self.plan_id) + \
            xml_line(2, "UserName", self.user_name) + \
            xml_line(2, "Description", self.description) + \
            xml_line(2, "Reason", self.reason) + \
            xml_line(2, "Project", self.project_name) + \
            xml_open(2, "Schedule")
        # TODO: airmass vs airmass range vs horizon probably s/b handled with @property
        #       ...(so that setting one turns off the others).
        this_airmass_minimum = self.airmass_range_minimum
        this_airmass_maximum = self.airmass
        if (self.airmass is not None) and (self.airmass_range_maximum is not None):
            this_airmass_maximum = min(self.airmass, self.airmass_range_maximum)
        if this_airmass_minimum is not None:
            rtml_text += xml_open(3, "AirmassRange") + xml_line(4, "Minimum", this_airmass_minimum)
            if this_airmass_maximum is not None:
                rtml_text += xml_line(4, "Maximum", this_airmass_maximum)
            rtml_text += xml_close(3, "AirmassRange")
        else:
            if this_airmass_maximum is not None:
                rtml_text += xml_open(3, "Airmass") + \
                             xml_line(4, "Airmass", this_airmass_maximum) +\
                             xml_close(3, "Airmass")
        if self.horizon is not None:
            rtml_text += xml_line(3, "Horizon", str(self.horizon))
        if (self.hour_angle_range_east is not None) and (self.hour_angle_range_west is not None):
            rtml_text += xml_open(3, "HourAngleRange") + \
                xml_line(4, "East", "{0:.2f}".format(self.hour_angle_range_east)) + \
                xml_line(4, "West", "{0:.2f}".format(self.hour_angle_range_west)) + \
                xml_line(3, "HourAngleRange")
        rtml_text += xml_line(3, "SkyCondition", self.sky_condition)
        rtml_text += xml_open(3, "Moon") + \
            xml_line(4, "Distance", "{0:.2f}".format(self.moon_distance)) + \
            xml_line(4, "Width", "{0:.2f}".format(self.moon_width)) + \
            xml_close(3, "Moon")
        if (self.earliest is not None) or (self.latest is not None):
            rtml_text += xml_open(3, "TimeRange")
            if self.earliest is not None:
                rtml_text += xml_line(4, "Earliest", self.earliest)
            if self.latest is not None:
                rtml_text += xml_line(4, "Latest", self.latest)
            rtml_text += xml_close(3, "TimeRange")
        rtml_text += xml_line(3, "Priority", "{0:d}".format(self.priority))
        rtml_text += xml_close(2, "Schedule")
        rtml_text += xml_line(2, "Telescope", self.telescope)
        for observation in self.observation_list:
            rtml_text += observation.rtml()
        rtml_text += xml_close(1, "Request")
        return rtml_text


class ACPS_observation:
    """
    Object: holds one observation (RTML Target)
    Usage: obs = ACPS_observation('ST Tri', 34.555, +21.334)
           obs.add_imageset('ST Tri, 3, 120, 1, 'V')
    """
    def __init__(self, obs_id, RA_deg, dec_deg, autofocus=False, count=1):
        if obs_id is None:
            raise ValueError('observation id may not be null.')
        if len(obs_id) <= 0:
            raise ValueError('observation id may not be zero-length.')
        self.id = obs_id
        self.ra = RA_deg
        self.dec = dec_deg
        self.autofocus = autofocus
        self.count = count
        self.imageset_list = []

    def add_imageset(self, name, count, exposure, filter):
        # TODO: construct default name if given imageset name == "" or None
        imageset_to_add = ACPS_imageset(name, count, exposure, filter)
        self.imageset_list.append(imageset_to_add)

    def rtml(self):
        open_text = 'Target count=\"' + "{0:d}".format(self.count) + '\"'
        if self.autofocus is True:
            open_text += ' autofocus=\"True\"'
        rtml_text = xml_open(3, open_text)
        rtml_text += xml_line(3, "ID", self.id)
        rtml_text += xml_line(3, "Name", self.id)  # ID == Name, for now.
        rtml_text += xml_open(3, "Coordinates") + \
            xml_line(4, "RightAscension", "{0:.4f}".format(self.ra)) + \
            xml_line(4, "Declination", "{0:.4f}".format(self.dec)) + \
            xml_close(4, "Coordinates")
        for imageset in self.imageset_list:
            rtml_text += imageset.rtml()
        rtml_text += xml_close(2, 'Target')
        return rtml_text


class ACPS_imageset:
    """
    Pseudo-private class.
    Object: holds one image set (RTML Picture)
    Usage: is = ACPS_imageset('ST Tri', 3, 120, 'V')
           acps_obs.add_imageset(is)
    """
    def __init__(self, name, count, exposure, filter):
        self.name = name
        self.count = count
        self.exposure = max(0, exposure)
        self.binning = 1  # ALWAYS BINNING == 1
        self.filter = filter

    def rtml(self):
        rtml_text = xml_open(3, 'Picture count=\"' + "{0:d}".format(self.count) + '\"')
        rtml_text += xml_line(4, 'Name', self.name)
        rtml_text += xml_line(4, 'Description', self.name)  # Name == Description, for now.
        rtml_text += xml_line(4, 'ExposureTime', "{0:d}".format(self.exposure))
        rtml_text += xml_line(4, 'Binning', "{0:d}".format(self.binning))
        rtml_text += xml_line(4, 'Filter', self.filter)
        rtml_text += xml_line(4, 'Dither', '0')  # Dither always zero, for now.
        rtml_text += xml_close(3, 'Picture')
        return rtml_text


def xml_line(n_tabs, tag, content=''):
    return n_tabs * '\t' + '<' + tag + '>' + content + '</' + tag + '>\n'


def xml_open(n_tabs, tag):
    return n_tabs * '\t' + '<' + tag + '>\n'


def xml_close(n_tabs, tag):
    return n_tabs * '\t' + '</' + tag + '>\n'

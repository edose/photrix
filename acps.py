__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

class ACPS_project:
    def __init__(self, project_name=None, telescope="BOREA", user="Eric", email="astro@ericdose.com",
                 ):
        self.project_name = project_name
        self.telescope = telescope
        self.user = user
        self.email = email
        self.plan_list = []

    def make_plan(self):
        return ACPS_plan(self.project_name, self.telescope, self.user)

    def add_plan(self, new_acps_plan):
        if isinstance(new_acps_plan, ACPS_plan):
            self.plan_list.append(new_acps_plan)

    def rtml(self):
        rtml_text = '<?xml version="1.0"?>\n' + \
                    '<RTML xmlns:xs="http://www.w3.org/2001/XMLSchema-instance" version="2.3">\n' + \
                    xml_open(1, "Contact") + \
                    xml_line(2, "User", self.user) + \
                    xml_line(2, "Email", self.email) + \
                    xml_close(1, "Contact")
        for plan in self.plan_list:
            rtml_text += plan.rtml()
        return rtml_text


class ACPS_plan:
    def __init__(self, project_name, telescope, user):
        self.id = "This is the Request ID"
        self.user_name = user
        self.description = "This is the Description"
        self.reason = "This is the Reason"
        self.project_name = project_name
        self.airmass = None                # default
        self.airmass_range_minimum = None  # default
        self.airmass_range_maximum = None  # default
        self.horizon = None                # default, in degrees
        self.hour_angle_range_east = None  # default
        self.hour_angle_range_west = None  # default
        self.sky_condition = "Fair"        # default
        self.moon_distance = 60            # default, in degrees
        self.moon_width = 12               # default, in days
        self.earliest = None
        self.latest = None
        self.priority = 4            # default; NB: mean for all plans normalized to 0.5
        self.telescope = telescope
        self.observation_list = []

    def add_observation(self, ACPS_observation_to_add):
        pass

    def rtml(self):
        rtml_text = xml_open(1, "Request") + \
            xml_line(2, "ID", "This is the Request") + \
            xml_line(2, "UserName", self.user_name) + \
            xml_line(2, "Description", "This is the Description") + \
            xml_line(2, "Reason", "This is the Reason") + \
            xml_line(2, "Project", ) + \
            xml_open(2, "Schedule")
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
                rtml_text += xml_open(3, "Airmass") + xml_line(4, "Airmass", this_airmass_maximum) +\
                    xml_close(3, "Airmass")
        if self.horizon is not None:
            rtml_text += xml_line(3, "Horizon", self.horizon)
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
        rtml_text += xml_line(3, "Priority", self.priority)
        rtml_text += xml_close(2, "Schedule")
        rtml_text += xml_line(2, "Telescope", self.telescope)
        for observation in self.observation_list:
            rtml_text += observation.rtml()
        return rtml_text


class ACPS_observation:
    '''
    Object: holds one observation (RTML Target)
    Usage: obs = ACPS_observation('ST Tri', 34.555, +21.334)
           obs.add_imageset('ST Tri, 3, 120, 1, 'V')
    '''
    def __init__(self, id, RA_deg, dec_deg, autofocus, count):
        self.id = id
        self.ra = RA_deg
        self.dec = dec_deg
        self.autofocus = autofocus
        self.count = count
        self.imageset_list = []

    def add_imageset(self, name, count, exposure, binning, filter):
        imageset_to_add = _ACPS_imageset(name, count, exposure, binning, filter)
        self.imageset_list.append(imageset_to_add)

    def rtml(self):
        open_text = 'Target count=\"' + self.count + '\"'
        if self.autofocus==True:
            open_text += ' autofocus=\"True\"'
        rtml_text = xml_open(3, open_text)
        rtml_text += xml_line(3, "ID", self.id)
        rtml_text += xml_line(3, "Name", self.id)  # ID == Name, for now.
        rtml_text += xml_open(3, "Coordinates") + \
            xml_line(4, "RightAscension", self.ra) + \
            xml_line(4, "Declination", self.dec) + \
            xml_close(4, "Coordinates")
        for imageset in self.imageset_list:
            rtml_text += imageset.rtml()
        rtml_text += xml_close(3, 'Target')
        return rtml_text


class _ACPS_imageset:
    '''
    Pseudo-private class.
    Object: holds one image set (RTML Picture)
    Usage: is = ACPS_imageset('ST Tri', 3, 120, 1, 'V')
           acps_obs.add_imageset(is)
    '''
    def __init__(self, name, count, exposure, binning, filter):
        self.name = name
        self.count = count
        self.exposure = max(0, exposure)
        self.binning = binning
        self.filter = filter

    def rtml(self):
        rtml_text = xml_open(3, 'Picture count=\"' + self.count + '\"')
        rtml_text += xml_line(4, 'Name', self.name)
        rtml_text += xml_line(4, 'Description', self.name)  # Name == Description, for now.
        rtml_text += xml_line(4, 'ExposureTime', self.exposure)
        rtml_text += xml_line(4, 'Binning', self.binning)
        rtml_text += xml_line(4, 'Filter', self.filter)
        rtml_text += xml_close(3, 'Picture')
        return rtml_text


def xml_line (n_tabs, tag, content=''):
    return n_tabs * '\t' + '<' + tag + '>' + content + '</' + tag + '>\n'


def xml_open (n_tabs, tag):
    return n_tabs * '\t' + '<' + tag + '>\n'


def xml_close (n_tabs, tag):
    return n_tabs * '\t' + '</' + tag + '>\n'


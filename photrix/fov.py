from photrix.util import *

__author__ = "Eric Dose :: Bois d'Arc Observatory, Kansas"

PHOTOMETRY_PATH = "C:/Dev/Photometry/FOV/"


class FOV:
    """
    Object: holds info for one Field Of View (generally maps to one AAVSO sequence or chart).
    Usage: fov = FOV("ST Tri")
    """
    def __init__(self, fov_name):
        self.name = fov_name

        # For now (20160413), let's stick with std/legacy R/text file format.
        fov_relpath = PHOTOMETRY_PATH + fov_name + ".txt"
        fov_f = open(fov_relpath, 'r')
        lines = fov_f.readlines()
        fov_f.close()

        lines = [line.split(";")[0] for line in lines] # remove comments
        directive_lines = []
        for line in lines:
            # line = line.split(";")[0].strip()  # remove comments
            if line.startswith("#"):
                directive_lines.append(line)

        # Required directives here.
        self.sequence = self._directive_values(lines, "#SEQUENCE", 1)
        ra_str, dec_str = self._directive_values(lines, "#CENTER", 2)
        self.ra = ra_as_degrees(ra_str)
        self.dec = dec_as_degrees(dec_str)
        self.chart = self._directive_values(lines, "#CHART", 1)
        self.fov_date = self._directive_values(lines, "#DATE", 1)
        self.fov_type = self._directive_values(lines, "#TYPE", 1)
        acp_string = self._directive_values(lines, "#ACP", 1)
        if acp_string is not None:
            acp_words = acp_string.split("|")
            self.acp_comments = acp_words[len(acp_words) - 1]
            self.acp_directives = acp_words[:len(acp_words)-1]
        else:
            self.acp_directives = self._directive_values(lines, "#ACP_DIRECTIVES", 1)
            self.acp_comments = self._directive_values(lines, "#ACP_COMMENTS", 1)

        # Optional directives here: #PERIOD, #ACP_DIRECTIVES, #ACP_COMMENT, #STARE, #PUNCH.
        self.period = self._directive_values(lines, "#PERIOD", 1)
        self.stare = self._directive_values(lines, "#STAREFOR", 1)
        self.punches = self._get_punch_values(lines)

        # Get AAVSO star lines.
        self.aavso_stars = self._get_aavso_star_lines(lines)

    def _directive_values(self, lines, directive_string, n_values):
        default_value = None  # if directive absent from FOV file.
        for line in lines:
            if line.upper().startswith(directive_string):
                value_string = line[len(directive_string):].strip()
                if n_values == 1:
                    return value_string
                else:
                    value_strings = value_string.split()
                    if len(value_strings) >= n_values:
                        return [s.strip() for s in value_strings]
        return default_value

    def _get_punch_values(self, lines):
        punch_values = []
        for line in lines:
            if line.upper().startswith("#PUNCH"):
                value_string = line[len("#PUNCH"):].strip()
                # colon delimiter, as star ID may contain spaces (e.g., "ST Tri")
                star_id = (value_string.split(":")[0])
                terms = (value_string.split(":")[1]).split()
                if len(terms) == 2:
                    punch_item = [star_id.strip(), terms[0].strip(), terms[1].strip()]
                    punch_values.append(punch_item)
        return punch_values

    def _get_aavso_star_lines(self, lines):
        aavso_stars = []
        stars_line_found = False
        for line in lines:
            if line.upper().startswith("#STARS"):
                stars_line_found = True
            else:
                if stars_line_found:
                    aavso_star = AavsoSequenceStar(line)
                    if aavso_star.is_valid:
                        aavso_stars.append(aavso_star)
        return aavso_stars

    def __str__(self):
        return "FOV '" + self.name + "' with " + str(len(self.aavso_stars)) + " sequence stars."

    def available(self, astronight):
        """
        Returns Timespan object defining which this FOV may be observed during this astronight.
        Usage: ts = fov.available(astronight)
        """
        pass


class AavsoSequenceStar:
    """
    Object: holds parsed info for one star of an AAVSO sequence.
    Internal use only (by class FOV).
    """
    def __init__(self, star_string):
        words = star_string.strip().split("\t")   # AAVSO sequence star lines have tab delimiters.
        words = [word.strip() for word in words]  # strip whitespace from all items.
        self.star_id, self.star_type = ("", "")   # invalid defaults.

        if len(words) >= 7:
            # Get raw data strings from line.
            self.star_id = words[0]
            self.ra = float(words[1])   # already in degrees RA
            self.dec = float(words[2])  # already in degrees Dec

            # Extract star type.
            star_type_char = words[4]
            star_type_dict = {"C": "comp", "H": "check", "T": "target"}
            if star_type_char in star_type_dict.keys():
                self.star_type = star_type_dict[star_type_char]
            else:
                self.star_type = None

            # Extract magnitudes. (This code probably sub-optimal.)
            mag_string = words[3]
            self.magB, self.magV, self.magR, self.magI, self.magU = (None, None, None, None, None)
            if self.star_type != "target":
                mag_words = mag_string.split("|")
                mag_words = [mag_word.strip() for mag_word in mag_words]  # strip whitespace
                for mag_word in mag_words:
                    mag_split = mag_word.split("_")
                    this_mag = float(mag_split[1])
                    if this_mag != 0:
                        if mag_split[0] == "1":
                            self.magB = this_mag
                        if mag_split[0] == "2":
                            self.magV = this_mag
                        if mag_split[0] == "4":
                            self.magR = this_mag
                        if mag_split[0] == "8":
                            self.magI = this_mag
                        if mag_split[0] == "1024":
                            self.magU = this_mag

        # Finally, validate this object, or not.
        self.is_valid = (len(self.star_id) >= 2) & (self.star_type is not None)

    def __str__(self):
        return "AAVSO Sequence Star '" + self.star_id + "', is_valid=" + self.is_valid





# ---------------------------------------------------------------------

if __name__ == '__main__':
    f = FOV("ST Cnc")
    print(f.name, "\n", repr(f))

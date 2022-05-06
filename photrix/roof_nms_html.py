__author__ = "Eric Dose, Albuquerque"

""" This module: 
      
"""

# Python core:
import os

# Author's packages:
from roof_nms import monitor_nms_status_via_html


THIS_PACKAGE_ROOT_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INI_DIRECTORY = os.path.join(THIS_PACKAGE_ROOT_DIRECTORY, 'ini')


if __name__ == '__main__':
    # monitor_nms_status_via_e_mail()
    monitor_nms_status_via_html()
    # record_nms_status_image_name()

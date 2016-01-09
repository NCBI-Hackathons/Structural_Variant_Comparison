""" Dynamically generate summary for dbVar
"""

import ConfigParser, os
from jinja2 import Environment, FilesystemLoader

def generate_report(variable_dict):
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    j2_env = Environment(loader=FilesystemLoader)
    rpath = config.get('output', 'report_dir')
    vd = variable_dict


    out_html = html_string.format(**vd)
    f = open(rpath + 'report.html', 'w+')
    f.write(out_html)
    f.close()



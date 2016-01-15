""" Dynamically generate summary for dbVar
"""

import ConfigParser, os
from jinja2 import Template

def generate_report(variable_dict):
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    rpath = config.get('output', 'report_dir')

    template = Template(open('report_template.html', 'rU').read())
    #j2_env = Environment(loader=FilesystemLoader)
    vd = variable_dict
    f = open(rpath + 'report.html', 'w+')
    f.write(template.render(vd))
    f.close()

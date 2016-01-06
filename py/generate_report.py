""" Dynamically generate summary for dbVar
"""

import ConfigParser
import jinja2

def generate_report(variable_dict):
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    rpath = config.get('output', 'report_dir')
    vd = variable_dict
    html_string = '''
    <html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
    </head>
    <body>
    <h1> dbVar summary </h1>
    <h2> Variation Types</h2>
    {type_counts!s}
    <h3> Size distribution by type</h3>
    <img src='./copy number loss_size_distribution.png' />
    <img src='./copy number gain_size_distribution.png' />
    <img src='./inversion_size_distribution.png' />
    <h2> Study Statistics </h2>
    <h2> Unique Features </h2>
    </body>
    </html>
    '''

    out_html = html_string.format(**vd)
    f = open(rpath + 'report.html', 'w+')
    f.write(out_html)
    f.close()



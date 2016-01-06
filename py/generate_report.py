""" Dynamically generate summary for dbVar
"""

import ConfigParser

def generate_report(variable_dict):
    config = ConfigParser.RawConfigParser()
    config.read('../example.cfg')
    rpath = config.get('output', 'report_dir')
    vd = variable_dict
    html_string = '''
    <html>
    <head>
        <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">
        <style>body{ margin:0 100; background:whitesmoke; }</style>
    </head>
    <body>
    <h1> dbVar summary </h1>
    </body>
    </html>
    '''

    out_html = html_string.format()
    f.open(rpath + 'report.html')
    f.write(out_html)
    f.close()



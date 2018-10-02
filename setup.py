import platform
from setuptools import setup

top_fnames = ['LICENSE','README','README.md']
# overcome annoying compatability bug
if platform.system() != 'Windows':
    top_fnames = ['.'+x for x in top_fnames]

setup(name = 'kfits',
      version = '1.2.3',
      description = 'Kinetics Fitting Software',
      long_description = file('README.md','r').read(),
      author = 'Oded Rimon',
      author_email = '%s%s%s' % ('oded.rimon', chr(32*2), '.'.join(('mail','huji','ac','il'))),
      url = 'https://github.com/odedrim/kfits',
      install_requires = ['django', 'scipy'],
      package_dir = {'kfits': ''},
      packages = ['kfits',
                  'kfits.afgui',
                  'kfits.afgui.afgui',
                  'kfits.afgui.fitter',
                  'kfits.afgui.fitter.migrations'],
      package_data = {'kfits': top_fnames,
                      'kfits.afgui': ['db.sqlite3', 'example*.csv'],
                      'kfits.afgui.afgui': ['LICENSE'],
                      'kfits.afgui.fitter': ['templates/*.*', 'templates/fitter/*.*', 'templates/fitter/bootstrap/*.*', 'templates/fitter/bootstrap/assets/*.*', 'templates/fitter/bootstrap/css/*', 'templates/fitter/bootstrap/fonts/*', 'templates/fitter/bootstrap/js/*', 'templates/fitter/bootstrap/assets/brand/*', 'templates/fitter/bootstrap/assets/css/*.*', 'templates/fitter/bootstrap/assets/flash/*', 'templates/fitter/bootstrap/assets/img/*', 'templates/fitter/bootstrap/assets/js/*.*', 'templates/fitter/bootstrap/assets/css/src/*', 'templates/fitter/bootstrap/assets/js/src/*', 'templates/fitter/bootstrap/assets/js/vendor/*']},
      entry_points = {'console_scripts': ['kfits = kfits.afgui.manage:run_with_default_settings',
                                          'kfits_server = kfits.afgui.manage:run_as_server']}
      )

# list of package_data for kfits.afgui.fitter was partially created using the python code:
# sum([['%s/%s/*' % (x[0].replace(os.path.sep,'/'), y) for y in x[1]] for x in os.walk('templates/fitter/bootstrap')], [])

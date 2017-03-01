from setuptools import setup

setup(name = 'kfits',
      version = '0.1',
      description = 'Kinetics Fitting Software',
      long_description = file('README.md','r').read(),
      author = 'Oded Rimon',
      author_email = '%s%s%s' % ('oded.rimon', chr(32*2), '.'.join(('mail','huji','ac','il'))),
      url = 'https://github.com/odedrim/kfits',
      install_requires = {'django', 'scipy'},
      package_dir = {'kfits': ''},
      packages = ['kfits',
                  'kfits.afgui',
                  'kfits.afgui.afgui',
                  'kfits.afgui.fitter',
                  'kfits.afgui.fitter.migrations'],
      package_data = {'kfits': ['LICENSE', 'README'],
                      'kfits.afgui': ['db.sqlite3'],
                      'kfits.afgui.afgui': ['LICENSE'],
                      'kfits.afgui.fitter': ['templates/*.*', 'templates/fitter/*.*', 'templates/fitter/bootstrap/*', 'templates/fitter/bootstrap/assets/*', 'templates/fitter/bootstrap/css/*', 'templates/fitter/bootstrap/fonts/*', 'templates/fitter/bootstrap/js/*', 'templates/fitter/bootstrap/assets/brand/*', 'templates/fitter/bootstrap/assets/css/*', 'templates/fitter/bootstrap/assets/flash/*', 'templates/fitter/bootstrap/assets/img/*', 'templates/fitter/bootstrap/assets/js/*', 'templates/fitter/bootstrap/assets/css/src/*', 'templates/fitter/bootstrap/assets/js/src/*', 'templates/fitter/bootstrap/assets/js/vendor/*']},
      )

# list of package_data for kfits.afgui.fitter was partially created using the python code:
# sum([['%s/%s/*' % (x[0].replace(os.path.sep,'/'), y) for y in x[1]] for x in os.walk('templates/fitter/bootstrap')], [])

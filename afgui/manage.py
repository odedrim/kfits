#!/usr/bin/env python
import os
import sys

def main(args):
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "afgui.settings")
    try:
        from django.core.management import execute_from_command_line
    except ImportError:
        # The above import may fail for some other reason. Ensure that the
        # issue is really that Django is missing to avoid masking other
        # exceptions on Python 2.
        try:
            import django
        except ImportError:
            raise ImportError(
                "Couldn't import Django. Are you sure it's installed and "
                "available on your PYTHONPATH environment variable? Did you "
                "forget to activate a virtual environment?"
            )
        raise
    execute_from_command_line(args)

def run_with_default_settings():
    import inspect
    import subprocess
    import webbrowser
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    # THIS WAY DOESN'T WORK FOR SOME REASON
    #os.chdir(currentdir)
    #sys.path.insert(0, currentdir)
    #main(['manage.py', 'runserver'])
    # THIS DOES...
    subprocess.Popen(['python', 'manage.py', 'runserver'], cwd=currentdir)
    webbrowser.open_new_tab('http://127.0.0.1:8000/fitter')


if __name__ == "__main__":
    main(sys.argv)

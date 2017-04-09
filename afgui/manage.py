#!/usr/bin/env python
import os
import sys

def main(args, settings="afgui.settings"):
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", settings)
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
    # find the directory of `manage.py`
    import inspect
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    os.chdir(currentdir)
    sys.path.insert(0, currentdir)
    # run django server (noreload is necessary, because otherwise this command fails when not running from the same path as `manage.py`
    main(['', 'runserver', '--noreload'], "kfits.afgui.afgui.settings")
##    # ALTERNATIVE WAY
##    import subprocess
##    import webbrowser
##    subprocess.Popen(['python', 'manage.py', 'runserver'], cwd=currentdir)
##    webbrowser.open_new_tab('http://127.0.0.1:8000/fitter')

def run_as_server():
    # find the directory of `manage.py`
    import inspect
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    os.chdir(currentdir)
    sys.path.insert(0, currentdir)
    # get ip and port from user
    ip = raw_input("Enter IP to serve on (use 0.0.0.0 for all interfaces, 127.0.0.1 for local loop interface only): ")
    port = raw_input("Enter port to serve on (use 8000 if you don't have another preference): ")
    if not port.isdigit or int(port) < 1 or int(port) > 65535:
        print >> sys.stderr, "Port must be a number between 1-65535."
        port = raw_input("Try port again: ")
    addr = '%s:%s' % (ip, port)
    # set environment variable for new address
    if ip == '0.0.0.0':
        client_addr = '127.0.0.1:%s' % port
    else:
        client_addr = addr
    os.environ.setdefault("DJANGO_SERVER_ADDRESS", client_addr)
    # run django server (noreload is necessary, because otherwise this command fails when not running from the same path as `manage.py`
    main(['', 'runserver', '--noreload', addr], "kfits.afgui.afgui.settings")


if __name__ == "__main__":
    main(sys.argv)

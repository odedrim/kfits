"""
WSGI config for afgui project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.10/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "afgui.settings")

application = get_wsgi_application()

# get server addr
addr = os.environ.get("DJANGO_SERVER_ADDRESS", '127.0.0.1:8000')
# open browser
try:
    import webbrowser
    webbrowser.open_new_tab('http://%s/fitter' % addr)
except:
    import sys
    print >> sys.stderr, "Warning: Couldn't open web browser client"

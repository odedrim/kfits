import os
import mimetypes
import json
import time
import datetime
import random
import hashlib
import json
import inspect
# django
from django.shortcuts import render
from django import http
from django import template
from .forms import UploadFileForm
from .models import KfitsSession
# this module
import fbackend


TMP_FILENAME = 'kfits_%s.txt'


def _get_tmp_dir():
    return os.environ.get('TEMP', os.environ.get('TMP', '/tmp'))

# Create your views here.
def index(request):
    # create index
    tmplt = template.loader.get_template('fitter/index.htm')
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    return http.HttpResponse(tmplt.render(dict(tmp_file = '',
                                               example1 = os.path.join(parent_dir, 'example1.csv'),
                                               example2 = os.path.join(parent_dir, 'example2.csv'),
                                               example3 = os.path.join(parent_dir, 'example3.csv'),
                                               model_choice = [(model, fbackend.tfitter.FITTING_PAIRS[model][3]) for model in fbackend.tfitter.FITTING_PAIRS.keys()],
                                               width = fbackend.FIG_W,
                                               height = fbackend.FIG_H,
                                              ), request))

def test(request):
    tmplt = template.loader.get_template('fitter/test.htm')
    sessions = KfitsSession.objects.all()
    return http.HttpResponse(tmplt.render({'sessions': '<br>'.join(['%s (%s): %d bytes, %s' % (ses.sid, ses.timestamp.ctime(), len(ses.saved_output) if ses.saved_output else 0, ses.kinetic_params) for ses in sessions])}, request))

def bootstrap(request):
    response = file(os.path.join(os.getcwd(), 'fitter/templates', request.path.strip('/')),'rb').read()
    return http.HttpResponse(response, mimetypes.guess_type(request.path)[0])

def backend(request):
    if request.GET.has_key("function") and hasattr(fbackend, request.GET['function']):
        params = dict(request.GET)
        func = getattr(fbackend, params.pop('function')[0])
        # run function
        params.pop('_', None)
        res = func(**params)
        if res[0]:
            new_res = res[1:]
            if len(new_res) == 1:
                new_res = new_res[0]
            # return output
            return http.HttpResponse(json.dumps(new_res), "application/json")
        else:
            # return output
            response = http.HttpResponse(res[1], res[2])
            if len(res) > 3:
                response['Content-Disposition'] = 'attachment; filename="%s"' % res[3];
            return response
    else:
        return http.HttpResponseBadRequest()

def upload_text(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if form.is_valid():
            # get client ip
            forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
            client_ip = forwarded_for.split(',')[0] if forwarded_for else request.META.get('REMOTE_ADDR')
            # randomize session id
            full_session_details = '%s%.3f%.8f' % (client_ip, time.time(), random.random())
            session_id = hashlib.md5(full_session_details).hexdigest()
            # create session
            session = KfitsSession(sid = session_id)
            session.save()
            # delete old sessions
            old_sessions = KfitsSession.objects.filter(timestamp__lte=datetime.datetime.now()-datetime.timedelta(1,0,0))
            for oldses in old_sessions:
                # remove temporary file
                sfname = os.path.join(_get_tmp_dir(), TMP_FILENAME % oldses.sid)
                try:
                    os.remove(sfname)
                except OSError:
                    pass
                # delete entry from database
                oldses.delete()
            # make form
            fname = os.path.join(_get_tmp_dir(), TMP_FILENAME % session_id)
            f = file(fname, 'w')
            for chunk in request.FILES['fdata'].chunks():
                f.write(chunk)
            return http.HttpResponse(json.dumps(fname), "application/json")
    else:
        return render(request, 'upload.htm', {'form': UploadFileForm()})
    return http.HttpResponse(json.dumps(False), "application/json")

def get_past_results(request):
    if request.GET.has_key("sid"):
        sid = request.GET['sid']
        try:
            session = KfitsSession.objects.get(sid=sid)
        except Exception:
            return http.HttpResponseNotFound('<h1>Cannot find session %s</h1>It was either deleted or never existed.' % sid)
        if request.GET.has_key("download"):
            data = session.saved_output
            response = http.HttpResponse(data, 'text/plain')
            response['Content-Disposition'] = 'attachment; filename="%s.csv"' % sid
            return response
        else:
            if session.kinetic_params:
                tmplt = template.loader.get_template('fitter/get.htm')
                return http.HttpResponse(tmplt.render(dict(kparams = [dict(zip(['name','value'],param)) for param in sorted(json.loads(session.kinetic_params).items())],
                    sid = sid)))
            else:
                return http.HttpResponseNotFound('<h1>Cannot find results for session %s</h1>The session exists, but cleaning results were not saved by the user.' % sid)
    else:
        tmplt = template.loader.get_template('fitter/get_nosid.htm')
        return http.HttpResponse(tmplt.render())

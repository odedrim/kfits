import os
import mimetypes
import json
import inspect
# django
from django.shortcuts import render
from django import http
from django import template
from .forms import UploadFileForm
# this module
import fbackend


TMP_FILENAME = 'kfits.txt'


def _get_tmp_dir():
    return os.environ.get('TEMP', os.environ.get('TMP', '/tmp'))

# Create your views here.
def index(request):
    tmplt = template.loader.get_template('fitter/index.htm')
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
    return http.HttpResponse(tmplt.render(dict(tmp_file = os.path.join(_get_tmp_dir(), TMP_FILENAME),
                                               example1 = os.path.join(parent_dir, 'example1.csv'),
                                               example2 = os.path.join(parent_dir, 'example2.csv'),
                                               example3 = os.path.join(parent_dir, 'example3.csv'),
                                               model_choice = [(model, fbackend.tfitter.FITTING_PAIRS[model][3]) for model in fbackend.tfitter.FITTING_PAIRS.keys()],
                                              ), request))

def test(request):
    tmplt = template.loader.get_template('fitter/test.htm')
    return http.HttpResponse(tmplt.render({}, request))

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
            fname = os.path.join(_get_tmp_dir(), TMP_FILENAME)
            f = file(fname, 'w')
            for chunk in request.FILES['fdata'].chunks():
                f.write(chunk)
            return http.HttpResponse(json.dumps(fname), "application/json")
    else:
        return render(request, 'upload.htm', {'form': UploadFileForm()})
    return http.HttpResponse(json.dumps(False), "application/json")

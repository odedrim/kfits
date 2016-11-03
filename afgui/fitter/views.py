import os
import mimetypes
import json
# django
from django.shortcuts import render
from django import http
from django import template
# this module
import fbackend

# Create your views here.
def index(request):
    tmplt = template.loader.get_template('fitter/index.htm')
    return http.HttpResponse(tmplt.render({}, request))

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
        print func, params
        res = func(**params)
        if res[0]:
            new_res = res[1:]
            if len(new_res) == 1:
                new_res = new_res[0]
            # return output
            return http.HttpResponse(json.dumps(new_res), "application/json")
        else:
            # return output
            return http.HttpResponse(res[1], res[2])
    else:
        return http.HttpResponseBadRequest()

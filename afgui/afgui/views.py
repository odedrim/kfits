import os
# django
from django.shortcuts import redirect

MAIN_PROJECT = 'fitter/'

# Create your views here.
def index(request):
    return redirect(MAIN_PROJECT)

import os
import scipy
import time
import math
import json
# this package
import tfitter
from .models import KfitsSession


##########
# CONSTS #
##########
# dimensions of visualisation
FIG_W = 870
FIG_H = 580
SMALL_FIG_W = 520
SMALL_FIG_H = 370
# input file params
MAX_FILE_SIZE = 10*1024*1024
# noise reduction params
MAX_THRESHOLD_SEARCH_ITERATIONS = 10
MAX_NOISE_REMOVAL = 0.4
NOISE_MAX_CONCENTRATION = 1.0


def check_input_file(fnames):
    TO_JSON = True
    fname = fnames[0]
    if not os.access(fname, os.R_OK) or os.path.isdir(fname):
        return TO_JSON, False, "Permission Denied or Directory"
    if os.stat(fname).st_size > MAX_FILE_SIZE:
        return TO_JSON, False, "I am not allowed to read files bigger than %.1f MB" % (MAX_FILE_SIZE / 1048576.)
    data = file(fname,'rb').read()
    if 'XYDATA' not in data:
        return TO_JSON, False, "This file does not seem to contain aggregation data in CSV or TXT format."
    # if we got all the way here...
    return TO_JSON, True, "OK"

def _calibrate_threshold_points(points, w, h):
    w, h = float(w), float(h)
    res = []
    for x,y in points:
        res.append((x/w, 1.-y/h))
    return res

def plot_data(fnames, threshold_points=None, rev_threshold_points=None):
    # prepare return type
    TO_JSON = False
    RETTYPE = 'image/svg+xml'
    # parse parameters
    fname = fnames[0]
    threshold_points = _calibrate_threshold_points(eval(threshold_points[0]), FIG_W, FIG_H) if threshold_points else None
    rev_threshold_points = _calibrate_threshold_points(eval(rev_threshold_points[0]), FIG_W, FIG_H) if rev_threshold_points else None
    # run
    data = scipy.array(tfitter.parse_fluorometer_csv(fname, threshold_points, rev_threshold_points)).T
    return TO_JSON, tfitter.plot_to_svg(data[0], data[1], FIG_W, FIG_H), RETTYPE

def _fit_data(fnames, model, threshold_points=None, rev_threshold_points=None, approx_start=None, search_for_end=True):
    # parse parameters
    fname = fnames[0]
    threshold_points = _calibrate_threshold_points(eval(threshold_points[0]), FIG_W, FIG_H) if threshold_points else None
    rev_threshold_points = _calibrate_threshold_points(eval(rev_threshold_points[0]), FIG_W, FIG_H) if rev_threshold_points else None
    approx_start = float(approx_start[0])/FIG_W if approx_start else 0
    # run
    data = tfitter.parse_fluorometer_csv(fname, threshold_points, rev_threshold_points)
    sim_data, params = tfitter.fit_data(data, approx_start=approx_start*data[-1][0], search_for_end=search_for_end, model=model)
    return data, sim_data, params

def fit_data(fnames, model, threshold_points=None, rev_threshold_points=None, approx_start=None, search_for_end=None):
    # prepare return type
    TO_JSON = False
    RETTYPE = 'image/svg+xml'
    # run
    data, sim_data, params = _fit_data(fnames, model[0], threshold_points, rev_threshold_points, approx_start, search_for_end is not None)
    adata = scipy.array(data).T
    return TO_JSON, tfitter.plot_two_to_svg(adata[0], adata[1], [0.5*i for i in xrange(len(sim_data))], sim_data, FIG_W, FIG_H), RETTYPE

def _clean_data_short(data, sim_data, noise_threshold, noise_only_above=True):
    if noise_only_above:
        fltr = lambda (t,v): v-sim_data[int(t*2)] < noise_threshold
    else:
        fltr = lambda (t,v): abs(v-sim_data[int(t*2)]) < noise_threshold
    new_data = filter(fltr, data)
    adata = scipy.array(new_data).T
    return adata, new_data

def _clean_data(fnames, model, noise_threshold, threshold_points=None, rev_threshold_points=None, approx_start=None, noise_only_above=True, search_for_end=None):
    # parse parameters
    if isinstance(noise_threshold[0], int):
        noise_threshold = noise_threshold[0]
    elif noise_threshold[0].isdigit():
        noise_threshold = int(noise_threshold[0])
    elif noise_threshold[0].replace('.','').isdigit():
        noise_threshold = int(float(noise_threshold[0]))
    else:
        raise ValueError("Noise threshold has illegal value %r" % noise_threshold[0])
    # run fitting
    data, sim_data, params = _fit_data(fnames, model[0], threshold_points, rev_threshold_points, approx_start, search_for_end is not None)
    # find basal level
    baseline = tfitter.find_absolute_baseline(data)
    # clean data
    adata, new_data = _clean_data_short(data, sim_data, noise_threshold, noise_only_above)
    return adata, new_data, baseline, params

def clean_data(fnames, model, noise_threshold, threshold_points=None, rev_threshold_points=None, approx_start=None, noise_only_above=None, search_for_end=None):
    # prepare return type
    TO_JSON = True
    RETM = "OK"
    # run clean
    noise_only_above = noise_only_above is not None
    adata, new_data, baseline, params = _clean_data(fnames, model, noise_threshold, threshold_points, rev_threshold_points, approx_start, noise_only_above, search_for_end)
    # return clean figure and fitting parameters
    return TO_JSON, [tfitter.plot_to_svg(adata[0], adata[1], SMALL_FIG_W, SMALL_FIG_H, classes=['kfigure']), noise_threshold, params.pop('model'), params.pop('t1'), params.pop('t2'), params], RETM

def clean_data_optimise_noise_threshold(fnames, model, threshold_points=None, rev_threshold_points=None, approx_start=None, noise_only_above=None, search_for_end=None):
    # prepare return type
    TO_JSON = True
    RETM = "OK"
    # parse params
    noise_only_above = noise_only_above is not None
    # optimisation loop
    data, sim_data, params = _fit_data(fnames, model[0], threshold_points, rev_threshold_points, approx_start, search_for_end is not None)
    span = max(data, key=lambda (t,v):v)[1] - min(data, key=lambda (t,v):v)[1]
    noise_threshold = span / 3.
    ##last_reduction = 0
    for i in xrange(MAX_THRESHOLD_SEARCH_ITERATIONS):
        adata1, new_data1 = _clean_data_short(data, sim_data, noise_threshold, noise_only_above)
        adata2, new_data2 = _clean_data_short(data, sim_data, noise_threshold/2., noise_only_above)
        diff = set(new_data1) - set(new_data2)
        th2rem = set(data) - set(new_data2)
        # TODO: stop iteration if gain diminishes
        # test if the half threshold is too small
        if len(new_data2) < (1 - MAX_NOISE_REMOVAL) * len(data) or \
                tfitter.get_point_max_concentration(sorted(th2rem), 10)[2] > NOISE_MAX_CONCENTRATION:
            noise_threshold = noise_threshold * 1.5
        # test if halving the threshold is contributing or detrimental
        if len(diff) > 0.7 * MAX_NOISE_REMOVAL * len(data) or \
                tfitter.get_point_max_concentration(sorted(diff), 10)[2] > 0.7 * NOISE_MAX_CONCENTRATION:
            # probably detrimental - go half way
            noise_threshold *= 0.75
        else:
            noise_threshold *= 0.5
        if noise_threshold < 1.:
            # this has gone far enough...
            break
    # final run!
    noise_threshold = int(noise_threshold+0.5)
    adata, new_data, baseline, params = _clean_data(fnames, model, [noise_threshold], threshold_points, rev_threshold_points, approx_start, noise_only_above, search_for_end)
    # return clean figure and fitting parameters
    return TO_JSON, [tfitter.plot_to_svg(adata[0], adata[1], SMALL_FIG_W, SMALL_FIG_H, classes=['kfigure']), noise_threshold, params.pop('model'), params.pop('t1'), params.pop('t2'), params], RETM

def get_clean_data(fnames, model, noise_threshold, output_fnames, threshold_points=None, rev_threshold_points=None, approx_start=None, noise_only_above=None, search_for_end=None):
    # prepare return type
    TO_JSON = False
    RETTYPE = 'text/plain'
    # parse params
    noise_only_above = noise_only_above is not None
    # run clean
    adata, new_data, baseline, params = _clean_data(fnames, model, noise_threshold, threshold_points, rev_threshold_points, approx_start, noise_only_above, search_for_end)
    # prepare file
    new = ['XYDATA,values']
    for x,y in new_data:
        new.append('%f,%f' % (x,y-baseline))
    data = '\n'.join(new)
    # save output to database
    sid = fnames[0].rsplit('.',1)[0][-32:]
    try:
        session = KfitsSession.objects.get(sid=sid)
        session.saved_output = data
        session.kinetic_params = json.dumps(params)
        session.save()
    except Exception:
        pass
    # return clean figure and fitting parameters
    return TO_JSON, data, RETTYPE, output_fnames[0]


#########################
# Non-algorithmic Stuff #
#########################

def get_email_oded():
    return False, '.'.join(('oded','rimon123mail','huji','ac','il')), 'text/plain'

def get_email_dana():
    return False, '.'.join(('danare123mail','huji','ac','il')), 'text/plain'

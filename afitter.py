#!/usr/bin/env python
import os
import sys
import math
import scipy
from scipy.optimize import curve_fit
# for debugging
import inspect
# this package
import kmodels

#########
# UTILS #
#########

def _debug(d):
    last_func_name = _debug.__dict__.setdefault('last_func', None)
    func_name = inspect.getouterframes(inspect.currentframe(), 2)[1][3]
    if last_func_name != func_name:
        print >> sys.stderr
        print >> sys.stderr, func_name
    _debug.__dict__['last_func'] = func_name
    for k,v in d.iteritems():
        print >> sys.stderr, '\t',k,'=',repr(v)

def check_threshold(point, threshold_points):
    """
    @param point: Point to check
    @type point: 2-tuple
    @param threshold_points: Points on the graph defining thresholds for noise (optional, default: none)
    @type threshold_points: list of 2-tuples

    @return: Is the point above the line defined by the threshold points? (if it is neither above nor below - but left or right to the line - will return None)
    @rtype: bool or None
    """
    if point[0] <= threshold_points[0][0] or point[0] > threshold_points[-1][0]:
        return None
    for i in xrange(len(threshold_points)-1):
        x1, y1 = map(float, threshold_points[i])
        x2, y2 = map(float, threshold_points[i+1])
        if point[0] > x1 and point[0] <= x2:
            rel = (point[0] - x1) / (x2 - x1)
            return point[1] <= y1 + (y2 - y1) * rel
    raise RuntimeError("How did we get here?!")

def get_point_max_concentration(data, min_window=3):
    """
    @param data: Filtered aggregation data
    @type data: List of 2-tuples
    @param min_window: Minimum number of points on which we can calculate concentrations [default: 3]
    @type min_window: int

    @return: The maximum concentration of points found (at points/sec), and where they were found - (start, end, concentration)
    @rtype: tuple (int, int, float)
    """
    # check params
    if len(data) <= min_window:
        return (0, len(data), 0., 0.)
    # init
    top = (0, 0, 0., 0.)
    start = 0
    end = min_window
    cur_con = float(end) / (data[end-1][0] - data[start][0])
    cur_score = float(end)**0.5 * cur_con
    while end < len(data):
        # calculate scores for adding one point
        advcon = float(end-start+1) / (data[end][0] - data[start][0])
        # adjust score to number of points
        adv = float(end-start+1)**0.5 * advcon
        # should we advance the end index?
        if adv >= cur_score:
            # yes!
            cur_score = adv
            cur_con = advcon
            end += 1
        else:
            # no - but have we reached a top scoring set?
            if cur_score > top[3]:
                # new record!
                top = (start, end, cur_con, cur_score)
            end += 1
            start = end - min_window
    # have we reached a top scoring set at the end?
    if cur_score > top[3]:
        # new record!
        return (start, end, cur_con)
    else:
        return top[:3]


################
# FILE PARSERS #
################

def parse_fluorometer_csv(path, threshold_points=None, rev_threshold_points=None, debug=False):
    """
    Parse a CSV or TXT file created by the fluorometer.

    @param path: The path of the file
    @type path: str
    @param threshold_points: Points on the graph defining lower thresholds for noise (optional, default: none)
    @type threshold_points: list of 2-tuples
    @param rev_threshold_points: Points on the graph defining upper thresholds for noise (optional, default: none)
    @type rev_threshold_points: list of 2-tuples

    @return: Aggregation Data
    @rtype: List of 2-tuples
    """
    if not os.path.isfile(path):
        raise TypeError("%s is not a path to a regular file" % path)
    elif not os.access(path, os.R_OK):
        raise TypeError("Do not have reading permissions for %s" % path)
    else:
        data = file(path, 'r').read()
    if debug:
        _debug(dict(path=path, threshold_points=threshold_points, rev_threshold_points=rev_threshold_points))
    retval = []
    started = False
    ended = False
    ncomma, ntab, nlines = data.count(','), data.count('\t'), data.count('\n')
    if ncomma < ntab and ntab > nlines / 2.:
        sep = '\t'
    else:
        sep = ','
    for line in data.splitlines():
        if 'XYDATA' in line:
            started = True
        elif started and not ended:
            if line.strip() == "":
                ended = True
            else:
                retval.append(tuple(map(float, line.split(sep))))
    if not threshold_points and not rev_threshold_points:
        return retval
    else:
        tmp = scipy.array(retval)
        min_x, min_y = tmp.min(0)
        max_x, max_y = tmp.max(0)
        x_span = max_x - min_x
        y_span = max_y - min_y
        if debug:
            _debug({"points before filtering\n": len(retval)})
        if threshold_points:
            new_tps = []
            for x,y in threshold_points:
                new_tps.append((min_x+x*x_span, min_y+y*y_span))
            retval = filter(lambda point: check_threshold(point, new_tps) in (True, None), retval)
            if debug:
                _debug({"points after filtering by green line\n": len(retval)})
        if rev_threshold_points:
            new_tps = []
            for x,y in rev_threshold_points:
                new_tps.append((min_x+x*x_span, min_y+y*y_span))
            retval = filter(lambda point: check_threshold(point, new_tps) in (False, None), retval)
            if debug:
                _debug({"points after filtering by red line\n": len(retval)})
        return retval


################
# DATA FITTING #
################

def closest_t(data, t):
    """
    Return the index of the line with the closest t (time) value to the given t.
    Assumes the data is sorted by time. Returns None in case of failure.
    """
    if data is None or len(data) == 0:
        return
    cur = int(t*2)
    if cur >= len(data) - 1:
        cur = len(data) - 2
    if cur < 0:
        cur = 0
    inc = (-1,1)[data[cur][0] < t]
    for i in xrange(len(data)):
        if abs(data[cur][0] - t) < abs(data[cur-1][0] - t) and \
           abs(data[cur][0] - t) < abs(data[cur+1][0] - t):
            return cur
        else:
            cur += inc
            if cur == 0:
                return
            elif cur == len(data)-1:
                return cur
    return None

def find_absolute_baseline(data, debug=False):
    """
    Find the absolute baseline of the measurement.

    @param data: The aggregation data
    @type data: List of 2-tuples
    @param debug: Show debug information? (default: no)
    @type debug: bool

    @return: Basal value
    @rtype: float
    """
    # find first 30 seconds of measurement
    tmp_data = scipy.array(data[:60])
    mean = tmp_data.mean(0)[1]
    stdev = tmp_data.std(0)[1]
    # filter out "bad" data points
    filtered_data = scipy.array(filter(lambda (t,v): abs(v-mean) < stdev, tmp_data))
    better_mean = filtered_data.mean(0)[1]
    if debug:
        _debug(dict(mean=mean, stdev=stdev, better_mean=better_mean, better_stdev=filtered_data.std(0)[1]))
    return better_mean

def find_aggregation_initiation(data, approx=120, debug=False):
    """
    Find the beginning of the aggregation process, i.e. the point in time
    where the actual experiment began, i.e. the time when the substrate was
    added to the mix.

    @param data: The aggregation data
    @type data: List of 2-tuples
    @param approx: Approximate time to search around (default: 120 sec)
    @type approx: number
    @param debug: Show debug information? (default: no)
    @type debug: bool

    @return: Basal value and Time of aggregation beginning (baseline, t1)
    @rtype: tuple (float, float)
    """
    if debug:
        _debug(dict(data_len=len(data), beginning=data[:10], approx=approx))
    RANGE_WIDTH = 40
    # find last 30 seconds before approximate init
    end_idx = closest_t(data, approx)
    if end_idx is None:
        end_idx = len(data) - 1
    if debug:
        _debug({'end index':end_idx})
    tmp_data = scipy.array(data[end_idx - 60 if end_idx > 60 else 0:end_idx + 1])
    mean = tmp_data.mean(0)[1]
    stdev = tmp_data.std(0)[1]
    # filter out "bad" data points
    filtered_data = scipy.array(filter(lambda (t,v): abs(v-mean) < stdev, tmp_data))
    better_mean = filtered_data.mean(0)[1]
    better_stdev = filtered_data.std(0)[1]
    if debug:
        _debug(dict(mean=mean, stdev=stdev, better_mean=better_mean, better_stdev=better_stdev))
    # locate t1
    best_idx = end_idx
    range_width = RANGE_WIDTH if RANGE_WIDTH < end_idx else end_idx
    for i in xrange(end_idx - range_width, end_idx + range_width + 1):
        smean = sum([x[1] for x in data[i:i+10]]) / 10.
        if smean > better_mean + 1.5 * better_stdev:
            best_idx = i
            break
    if debug:
        _debug(dict(end_idx=end_idx, best_idx=best_idx, smean=smean))
    return better_mean, data[best_idx][0]

def find_apparent_maximum(data, range_to_fit=120, minimal_t2=350., debug=False):
    """
    Find the point where aggregation reached its apparent maximum.
    This is often long before the end of the trace, as larger aggregates
    cannot be measured in some methods and effectively appear as plateaux
    or even slight decrease in measured values (for example in scattering).

    @param data: The aggregation data
    @type data: List of 2-tuples
    @param range_to_fit: How many timepoints at the end of the trace to fit with line (default: 120, i.e. 2 min)
    @type range_to_fit: int
    @param minimal_t2: Minimal value for t2, important for low-aggregation samples (default: 350s)
    @type minimal_t2: number
    @param debug: Show debug information? (default: no)
    @type debug: bool

    @return: The maximal point (t2, v_t2)
    @rtype: tuple (float, float)
    """
    DEVIATIONS_FROM_LINE = 10. 
    # fit final minute with line
    to_fit = data[-range_to_fit:]
    popt, pcov = curve_fit(lambda x, a, b: a*x+b, [x[0] for x in to_fit], [x[1] for x in to_fit])
    a, b = popt
    if debug:
        _debug(dict(fit_a=a, fit_b=b))
    ##perr = scipy.sqrt(scipy.diag(pcov))
    # calculate standard distance from line
    diffs = scipy.array([abs((a*x[0]+b)-x[1]) for x in data])
    fit_dev = diffs[-range_to_fit:].std()
    # forecast backwards using line, find point of divergence from linearity
    bins = [sum(diffs[i:i+10])/10. for i in xrange(0,len(diffs),10)]
    point_of_divergence = None
    while point_of_divergence is None or data[point_of_divergence][0] < minimal_t2:
        latest_divergent_bin = 0
        for i in xrange(len(bins)-2,-1,-1):
            if bins[i] > DEVIATIONS_FROM_LINE*fit_dev:
                latest_divergent_bin = i
                break
        if debug:
            _debug(dict(latest_divergent_bin=latest_divergent_bin, num_bins=len(bins)))
        for i in xrange(9,-1,-1):
            idx = latest_divergent_bin*10 + i
            if diffs[idx] > DEVIATIONS_FROM_LINE*fit_dev:
                point_of_divergence = idx
        DEVIATIONS_FROM_LINE /= 2.
    if debug:
        _debug(dict(point_of_divergence=point_of_divergence, retval=data[point_of_divergence]))
    return data[point_of_divergence]

def fit_aggregation_kinetics(data, baseline, t1, t2, apparent_max, fit_func, init_values, debug=False):
    """
    Fit the main aggregation period data with (t*vmax)/(t+t1/2)) kinetics.

    @param data: The aggregation data
    @type data: List of 2-tuples
    @param baseline: The calculated baseline of the curve
    @type baseline: float
    @param t1: Time of the beginning of the aggregation kinetics
    @type t1: number
    @param t2: Time of the end of the aggregation kinetics
    @type t2: number
    @param apparent_max: The apparent maximum value
    @type apparent_max: float
    @param fit_func: The fitting function; receives t (time, shifted to start from 0) and any parameters that should be optimised
    @type fit_func: function
    @param init_values: Initial values for all parameters of fit_func except t
    @type init_values: tuple of numbers
    @param debug: Show debug information? (default: no)
    @type debug: bool

    @return: popt and perr
    @rtype: tuple, tuple
    """
    to_fit = data[closest_t(data,t1):closest_t(data,t2)+1]
    popt, pcov = curve_fit(fit_func, [x[0]-t1 for x in to_fit], [x[1]-baseline for x in to_fit], p0=init_values, maxfev=10000)
    if debug:
        _debug(dict(popt=popt, pcov=pcov))
    perr = scipy.sqrt(scipy.diag(pcov))
    ##return [x[0]-t1 for x in to_fit], [x[1]-baseline for x in to_fit]
    return popt, perr

def choose_best_model(data, baseline, approx_start, t1, t2, apparent_max, fit_pairs, debug=False):
    """
    Fit the main aggregation period data with (t*vmax)/(t+t1/2)) kinetics.

    @param data: The aggregation data
    @type data: List of 2-tuples
    @param baseline: The calculated baseline of the curve
    @type baseline: float
    @param approx_start: Approximate time of beginning of kinetics, estimated by the user
    @type approx_start: number
    @param t1: Time of the beginning of the aggregation kinetics
    @type t1: number
    @param t2: Time of the end of the aggregation kinetics
    @type t2: number
    @param apparent_max: The apparent maximum value
    @type apparent_max: float
    @param fit_pairs: Dictionary of names pointing to 2-tuples, each comprised of a fitting function and initial values
    @type fit_pairs: dict (str:tuple)
    @param debug: Show debug information? (default: no)
    @type debug: bool

    @return: The best fit_pair's name
    @rtype: str
    """
    to_fit = data[closest_t(data,t1):closest_t(data,t2)+1]
    x_vals = [x[0]-t1 for x in to_fit]
    y_vals = [x[1]-baseline for x in to_fit]
    # fitting loop
    best = None
    best_score = None
    for name, (fit_func, init_values, param_names, human_name) in fit_pairs.iteritems():
        init_values = init_values(apparent_max, baseline)
        try:
            popt, pcov = curve_fit(fit_func, x_vals, y_vals, p0=init_values, maxfev=10000)
        except Exception, e:
            print >> sys.stderr, "Encountered error trying to fit data with %r: %s" % (name, e)
            continue
        score = sum([(y_vals[i] - fit_func(x_vals[i], *popt))**2 for i in xrange(len(y_vals)) if to_fit[i][0]>approx_start])
        if debug:
            _debug(dict(approx_start=approx_start, model_name=name, score=score, old_score=sum([abs(y_vals[i] - fit_func(x_vals[i], *popt)) for i in xrange(len(y_vals))])))
        if best is None or score < best_score:
            best = name
            best_score = score
    return best


#####################
# FITTING FUNCTIONS #
#####################

FITTING_PAIRS = kmodels.get_models()


#############
# INTERFACE #
#############

def fit_data(data, model='auto', approx_start=120, search_for_end=True, debug=False):
    # find initiation
    baseline, t1 = find_aggregation_initiation(data, approx=approx_start, debug=debug)
    # find plateau
    if search_for_end:
        t2, apparent_max = find_apparent_maximum(data, minimal_t2=(data[-1][0]+approx_start)/2, debug=debug)
    else:
        t2, apparent_max = data[-1]
    # choose model
    if model == 'auto':
        model = choose_best_model(data, baseline, approx_start, t1, t2, apparent_max, FITTING_PAIRS, debug=True)
        if model is None:
            raise RuntimeError("None of the models successfully fit the data!")
    fit_pair = FITTING_PAIRS[model]
    # fit kinetics
    popt, perr = fit_aggregation_kinetics(data, baseline, t1, t2, apparent_max, fit_pair[0], fit_pair[1](apparent_max, baseline), debug=debug)
    if debug:
        d = dict(zip(inspect.getargspec(fit_pair[0])[0][2:], ['%.3f +- %.3f' % (popt[i], perr[i]) for i in xrange(len(popt))]))
        _debug(d)
    # fit rest of data with straight line
    end_of_agg_curve = baseline + fit_pair[0](t2-t1, *popt)
    to_fit = data[closest_t(data,t2):]
    a, pcov = curve_fit(lambda x, a: a*x, [x[0]-t2 for x in to_fit], [x[1]-end_of_agg_curve for x in to_fit])
    a = float(a)

    # return fitted (simulated) data
    simulated_data = [baseline for x in xrange(int(t1*2))] + \
                     [baseline + fit_pair[0](x/2., *popt) for x in xrange(int((t2-t1)*2))] + \
                     [(a*((x/2.)-t2)+end_of_agg_curve) for x in range(int(t2*2.)-1, int(2*max(data)[0]))]
    # prepare return value
    params = dict(model = model, t1 = t1, t2 = t2)
    for i, param_name in enumerate(fit_pair[2]):
        params[param_name] = '%.3f &plusmn; %.3f' % (popt[i], perr[i])
    return simulated_data, params


##########################
# COMMAND-LINE INTERFACE #
##########################

def main(args):
    if len(args) < 1 or len(args) > 2:
        print >> sys.stderr, "Usage: %s <csv file path> (<approx start>)" % sys.argv[0]
        return 1
    try:
        from matplotlib import pyplot
    except ImportError, e:
        print >> sys.stderr, "The command line interface requires installation of the python matplotlib package."
        return 1
    data = parse_fluorometer_csv(args[0], debug=True)
    _debug({'data length':len(data)})
    if len(args) > 1:
        sim_data, params = fit_data(data, approx_start=int(args[1]), debug=True)
    else:
        sim_data, params = fit_data(data, debug=True)
    for item in params.iteritems():
        print '%s: %s' % item
    pyplot.plot([x[0] for x in data], [x[1] for x in data])
    pyplot.plot([0.5*i for i in xrange(len(sim_data))], sim_data)
    pyplot.show()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

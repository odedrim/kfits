#!/usr/bin/env python
import os
import sys
import afitter as fitter


DEFAULT_THRESHOLD = 25


def main(args):
    if len(args) < 1:
        print >> sys.stderr, "Usage: %s <csv file path> (<threshold>)" % sys.argv[0]
        return 1
    threshold = DEFAULT_THRESHOLD if len(args) < 2 else int(args[1])
    # read & fit data
    data = fitter.parse_fluorometer_csv(args[0])
    sim_data = fitter.fit_data(data)
    # find basal level
    baseline = fitter.find_absolute_baseline(data)
    # run fitting
    sim_data, params = fitter.fit_data(data)
    # clean data
    print 'threshold:', threshold
    new_data = filter(lambda (t,v): v-sim_data[int(t*2)] < threshold, data)
    # write to new file
    new_path = args[0].rsplit('.',1)[0] + '_clean.csv'
    w = file(new_path, 'w')
    w.write('Time (s),Light Scattering (A.U.)\n')
    w.writelines(['%.1f,%f\n' % (point[0], point[1]-baseline) for point in new_data])
    w.close()
    print "\nWRITTEN CLEAN DATA TO %s" % new_path
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

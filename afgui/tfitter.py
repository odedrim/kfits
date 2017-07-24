import os, sys
try:
    from kfits.afitter import *
except ImportError:
    # in case the package is not properly installed, try to work locally
    import inspect
    currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
    parentdir = os.path.dirname(currentdir)
    sys.path.insert(0,parentdir)
    from afitter import *

#########
# UTILS #
#########

SVG_TEMPLATE = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<svg xmlns="http://www.w3.org/2000/svg" width="%(width)d" height="%(height)d" %(classes)s viewBox="0 0 %(width)d %(height)d" preserveAspectRatio="none">
    <defs>
        <style type="text/css" />
    </defs>
    <g>
        %(gdata)s
    </g>
</svg>"""
SVG_PATH_TEMPLATE = """<path style="fill:none; stroke:%(color)s; stroke-width:1px"
    d="%(path)s" />"""

def plot_to_svg(xx, yy, w, h, color="#000000", classes=[]):
    min_x, max_x = min(xx), max(xx)
    min_y, max_y = min(yy), max(yy)
    norm_xx = (scipy.array(xx) - min_x) / (max_x - min_x) * w
    norm_yy = h - ((scipy.array(yy) - min_y) / (max_y - min_y) * h)
    retval = ["M %f %f" % (norm_xx[0], norm_yy[0])]
    for i in xrange(1, len(norm_xx)):
        retval.append('L %f %f' % (norm_xx[i], norm_yy[i]))
    path = SVG_PATH_TEMPLATE % dict(color = color,
                                    path = ' '.join(retval))
    return SVG_TEMPLATE % dict(width = w,
                               height = h,
                               gdata = path,
                               classes = ('class="%s"' % ' '.join(classes)) if classes else '')

def plot_two_to_svg(xx1, yy1, xx2, yy2, w, h, color1="#000000", color2="#ff4414"):
    min_x, max_x = min(min(xx1),min(xx2)), max(max(xx1),max(xx2))
    min_y, max_y = min(min(yy1),min(yy2)), max(max(yy1),max(yy2))
    norm_xx1 = (scipy.array(xx1) - min_x) / (max_x - min_x) * w
    norm_yy1 = h - ((scipy.array(yy1) - min_y) / (max_y - min_y) * h)
    norm_xx2 = (scipy.array(xx2) - min_x) / (max_x - min_x) * w
    norm_yy2 = h - ((scipy.array(yy2) - min_y) / (max_y - min_y) * h)
    retval = ["M %f %f" % (norm_xx1[0], norm_yy1[0])]
    for i in xrange(1, len(norm_xx1)):
        retval.append('L %f %f' % (norm_xx1[i], norm_yy1[i]))
    path1 = SVG_PATH_TEMPLATE % dict(color = color1,
                                     path = ' '.join(retval))
    retval = ["M %f %f" % (norm_xx2[0], norm_yy2[0])]
    for i in xrange(1, len(norm_xx2)):
        retval.append('L %f %f' % (norm_xx2[i], norm_yy2[i]))
    path2 = SVG_PATH_TEMPLATE % dict(color = color2,
                                     path = ' '.join(retval))
    return SVG_TEMPLATE % dict(width = w,
                               height = h,
                               gdata = path1 + '\n' + path2,
                               classes = '')

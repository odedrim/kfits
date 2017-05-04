import math
import scipy


FIT_BASIC_INIT = lambda apparent_max, baseline: (apparent_max-baseline, 12)
FIT_BASIC_PARAM_NAMES = ('v<sub>max</sub>', 't<sub>&#189;</sub>')
def fit_basic(t, vmax, thalf):
    return (t * vmax) / (t + thalf)

FIT_CAMBRIDGE_INIT = lambda apparent_max, baseline: (apparent_max-baseline, 1.1, 1000)
FIT_CAMBRIDGE_PARAM_NAMES = ('v<sub>max</sub>', 'n<sub>c</sub>', 'k')
class CambridgeFit(object):
    def __init__(self, m_total):
        self.mtot = m_total

    def fit_cambridge(self, t, vmax, nc, k):
        if nc < 0 or nc > 10 or k < 0 or k > 1e12:
            return 0
        ins = (nc * k * self.mtot**nc)**0.5
        try:
            outs = scipy.vectorize(math.cosh)(ins * t)
        except OverflowError, e:
            return 0
        return vmax * (1 - outs ** (-2./nc))


def get_models():
    return dict(basic = (fit_basic, \
                         FIT_BASIC_INIT, \
                         FIT_BASIC_PARAM_NAMES, \
                         'Basic One-Site Binding'),
                nucleation_elongation = (CambridgeFit(75e-6).fit_cambridge, \
                                         FIT_CAMBRIDGE_INIT, \
                                         FIT_CAMBRIDGE_PARAM_NAMES, \
                                         'Nucleation Elongation'),
                )


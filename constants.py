"""Just some global constants."""

import math

GeV2nb = 389379.;  # GeV -> nb conversion factor
GeV2fb = GeV2nb * 1e6

# SM

sw2 = 0.234
sw = math.sqrt(sw2)
cw2 = 1 - sw2
cw = math.sqrt(cw2)
v = 246
alpha = 1/127.9
mh = 125.09
GF = 1.1166e-5
mW = 80.385
mZ = 91.188
mt = 173.34

# BSM

mH = 750
alphas = 0.1 # at MH/2 or so
Kfactor = 2 # for ggF

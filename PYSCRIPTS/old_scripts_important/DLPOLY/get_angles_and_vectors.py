from __future__ import print_function, division
#import numpy as np
import sys
sys.path.append("/home/gadelmeier/Python/myPYMODULES/DLPOPLY_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/OTHER_MODULES/")
sys.path.append("/home/gadelmeier/Python/myPYMODULES/MATH_MODULES/")
import DL_HISTORY_class8_new_approaches_8_4 as dlhc

config, field = sys.argv

o = dlhc.ConfigFile(field, config)
o.read_config()
a, b, c, alpha, beta, gamma = o.config.box.convert_vectors()
print("\{{} {} {} {} {} {}\}".format(
    a, b, c, alpha, beta, gamma)
)

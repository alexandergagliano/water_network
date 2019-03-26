import numpy as np
import pandas as pd
import matplotlib.pyplot as pylab

colors = ["blue", "brown", "magenta", "purple"]
logZ = [-6., -4., -2., 0.]

for i in range(len(logZ)):
    filenamefH2 = 'logfH2_v_logN_[Z]=%d_2005.csv' % logZ[i]
    n = []
    fH2 = []
    with open(filenamefH2, 'r') as datafH2:
        lines = datafH2.readlines()
        for line in lines[1:]:
            p = line.split(',')
            n.append(float(p[0]))
            fH2.append(float(p[1]))
    
    pylab.plot(n, fH2, label='log (Z/Z$_{\odot}$) = %d' % logZ[i], color=colors[i])
    #pylab.axhline(y=0, color='r', linestyle='-')
    #print(n)
    #print(fH2)
    
    n = []
    fH2 = []
    if i == 0:
        pylab.ylim(ymin = -5, ymax = 0.5)
        pylab.xlim(xmin = -1.5, xmax=12)
        pylab.xlabel('number density log n$_{H}$ (cm$^{-3}$)')
        pylab.ylabel('H$_2$ fraction log f$_{H_{2}}$')
   
    ## OPEN OMUKAI 2005 FOR COMPARISON
    with open('logfH2_vs_logn_[ZH]=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
        lines = omukaiData.readlines()
        Omukai_n = []
        Omukai_fH2 = []
        for line in lines[6:]:
            p = line.split(',')
            Omukai_n.append(float(p[0]))
            Omukai_fH2.append(float(p[1]))
        pylab.plot(Omukai_n, Omukai_fH2, label='log (Z/Z$_{\odot}$) = %d, Omukai' % logZ[i], color=colors[i], linestyle=':')

        Omukai_n = []
        Omukai_fH2 = []

pylab.tight_layout()
#leg = pylab.legend(fancybox = True, labelspacing=0.0, loc='lower right')
#leg.get_frame().set_alpha(0.5)
pylab.savefig('logn-logfH2_GrackleWater_2005.png')

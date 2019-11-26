import numpy as np
import matplotlib.pyplot as plt


T, HYD, HI, HII, HM, H2I, H2II = np.loadtxt("abundances.txt", unpack=True, skiprows=1)

'''loglog(MET, C, label='C')
loglog(MET, O, label='O')
loglog(MET, CO, label='CO')
loglog(MET, OH, label='OH')
loglog(MET, H2O, label='H2O')
loglog(MET, O2, label='O2')
loglog(MET, CH, label='CH')
xlim(1.0-3,1.0)
ylim(1.0e-10, 1.0e-3)
legend(loc='best')
show()'''

fig, ax1 = plt.subplots()
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Fraction (Z/H)')
ax1.semilogy(T,HYD/HYD, label='Baseline')
ax1.semilogy(T,HI, label='HI')
ax1.semilogy(T,HII, label='HII')
ax1.semilogy(T,HM, label='HM')
ax1.semilogy(T,H2I, label='H2I')
ax1.semilogy(T,H2II,label='H2II')

ax2 = ax1.twinx()

ax2.set_ylabel('Total hydrogen number density (cm**-3)')
ax2.semilogy(T,HYD,'k--', label='Hydrogen')

fig.legend(loc='best')
fig.tight_layout()
plt.show()

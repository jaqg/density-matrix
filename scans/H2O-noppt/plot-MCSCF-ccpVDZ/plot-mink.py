# +---------------------------------------------+
# | Author: Jose Antonio Quiñonero Gris         |
# | Creation date: Saturday 15:24:39 05-11-2022 |
# +---------------------------------------------+

#
# Libreries
#
import matplotlib
from matplotlib.lines import MarkerStyle
import numpy as np
import matplotlib.pyplot as plt
import string
import os

matplotlib.use("pgf")

colores=['#56B4E9', '#E69F00', '#009E73', '#0072B2', '#D55E00', '#CC79A7', '#F0E442']
#
# Style
#
style_file = 'mine.mplstyle'
style_file = os.path.dirname(__file__)+'/{}'.format(style_file)
plt.style.use(style_file)

# Links
# https://jwalton.info/Embed-Publication-Matplotlib-Latex/
# https://jwalton.info/Matplotlib-latex-PGF/

#
# Input data
#
data_file = 'ord-1-minkowski'
data_file = os.path.dirname(__file__)+'/{}'.format(data_file)
imol, m1_els, m1_ebbc1, m1_ebbc2, m1_ebbc3, m1_ebbc3m = np.loadtxt(data_file,
                                                                        unpack=True,
                                                                        skiprows=1)
data_file = 'ord-2-minkowski'
data_file = os.path.dirname(__file__)+'/{}'.format(data_file)
imol, m2_els, m2_ebbc1, m2_ebbc2, m2_ebbc3, m2_ebbc3m = np.loadtxt(data_file,
                                                                        unpack=True,
                                                                        skiprows=1)
data_file = 'ord-3-minkowski'
data_file = os.path.dirname(__file__)+'/{}'.format(data_file)
imol, m3_els, m3_ebbc1, m3_ebbc2, m3_ebbc3, m3_ebbc3m = np.loadtxt(data_file,
                                                                        unpack=True,
                                                                        skiprows=1)
data_file = 'ord-4-minkowski'
data_file = os.path.dirname(__file__)+'/{}'.format(data_file)
imol, m4_els, m4_ebbc1, m4_ebbc2, m4_ebbc3, m4_ebbc3m = np.loadtxt(data_file,
                                                                        unpack=True,
                                                                        skiprows=1)
delta_d = 0.5  # angstrom
bd = delta_d * imol  # bond distances

#
# Plot
#
# fig, axs = plt.subplots(1,3,figsize=(7.2*1.5,4.45))  # For the pdf
fig, axs = plt.subplots(2,2,figsize=(6.33,2*3.45), sharex=True)        # For the pgf
fig.subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)

fig.suptitle(r'Minkowski distances for CASSCF(8,8)/cc-pVDZ with 1 inactive orbital for $\mathrm{H_2O}$')

# =============================================================================
# Energy plot
# -----------------------------------------------------------------------------
#
ax=axs[0,0]
#
ax.plot(bd, m1_els, label=r'LS')
ax.plot(bd, m1_ebbc1, label=r'BBC1')
ax.plot(bd, m1_ebbc2, label=r'BBC2')
ax.plot(bd, m1_ebbc3, label=r'BBC3')
ax.plot(bd, m1_ebbc3m, label=r'BBC3M')
ax.set(
        # title=r'Energy, $\Delta E$',
        # xlabel=r'Bond distance, $d$ (\AA)',
        ylabel=r'Minkowski 1st-order distance, $d^{(1)}$'
        )
ax.legend()

# =============================================================================
# Errors plot
# -----------------------------------------------------------------------------
#
ax=axs[0,1]
#
ax.plot(bd, m2_els, label=r'LS')
ax.plot(bd, m2_ebbc1, label=r'BBC1')
ax.plot(bd, m2_ebbc2, label=r'BBC2')
ax.plot(bd, m2_ebbc3, label=r'BBC3')
ax.plot(bd, m2_ebbc3m, label=r'BBC3M')
ax.set(
        # title=r'Energy, $\Delta E$',
        # xlabel=r'Bond distance, $d$ (\AA)',
        ylabel=r'Minkowski 2nd-order distance, $d^{(2)}$'
        )
ax.legend()

# =============================================================================
# Minkowski distances plot
# -----------------------------------------------------------------------------
#
ax=axs[1,0]
#
ax.plot(bd, m3_els, label=r'LS')
ax.plot(bd, m3_ebbc1, label=r'BBC1')
ax.plot(bd, m3_ebbc2, label=r'BBC2')
ax.plot(bd, m3_ebbc3, label=r'BBC3')
ax.plot(bd, m3_ebbc3m, label=r'BBC3M')
ax.set(
        # title=r'Energy, $\Delta E$',
        xlabel=r'Bond distance, $d$ (\AA)',
        ylabel=r'Minkowski 3rd-order distance, $d^{(3)}$'
        )
ax.legend()
# =============================================================================
# Minkowski distances plot
# -----------------------------------------------------------------------------
#
ax=axs[1,1]
#
ax.plot(bd, m4_els, label=r'LS')
ax.plot(bd, m4_ebbc1, label=r'BBC1')
ax.plot(bd, m4_ebbc2, label=r'BBC2')
ax.plot(bd, m4_ebbc3, label=r'BBC3')
ax.plot(bd, m4_ebbc3m, label=r'BBC3M')
ax.set(
        # title=r'Energy, $\Delta E$',
        xlabel=r'Bond distance, $d$ (\AA)',
        ylabel=r'Minkowski 4th-order distance, $d^{(4)}$'
        )
ax.legend()
# =============================================================================

# Enumerate the subplots
# for n, ax in enumerate(axs):
#     ax.text(-0.1, 1.1, string.ascii_uppercase[n], transform=ax.transAxes,
#             size=20, weight='bold')

nombre_grafica = os.path.basename(__file__).replace(".py", ".pdf")
plt.savefig(nombre_grafica, transparent='True', bbox_inches='tight')
#
nombre_grafica_2 = os.path.basename(__file__).replace(".py", ".pgf")
plt.savefig(nombre_grafica_2, format='pgf')

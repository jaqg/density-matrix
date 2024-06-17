# +-------------------------------------------+
# | Author: Jose Antonio Quinonero Gris       |
# | Creation date: Sunday 23:48:05 16/06/2024 |
# +-------------------------------------------+

#
# Libreries
#
import matplotlib
from matplotlib.lines import MarkerStyle
import numpy as np
import matplotlib.pyplot as plt
import string
import os

# matplotlib.use("pgf")

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
fichero_datos = 'all_energies.dat'
fichero_datos = os.path.dirname(__file__)+'/{}'.format(fichero_datos)
molnumb, ed1d2, els, bbc1, bbc2, bbc3, bbc3m = np.loadtxt(fichero_datos,
                                                          unpack=True, skiprows=1)

#
# Plot
#
fig, ax = plt.subplots()

ax.plot(molnumb, ed1d2, label=r'$E_{\text{exact}}$')
ax.plot(molnumb, els, label=r'$E_{\text{ELS}}$')
ax.plot(molnumb, bbc1, label=r'$E_{\text{BBC1}}$')
ax.plot(molnumb, bbc2, label=r'$E_{\text{BBC2}}$')
ax.plot(molnumb, bbc3, label=r'$E_{\text{BBC3}}$')
ax.plot(molnumb, bbc3m, label=r'$E_{\text{BBC3M}}$')

ax.set(
        title=r'H$_2$O',
        xlabel=r'Scan coordinate ($\mathrm{\AA}$)',
        ylabel=r'Energy, $E$ ($\mathrm{Ha}$)'
        )

ax.legend()
#
# ax.grid(False)
# plt.show()
nombre_grafica = os.path.basename(__file__).replace(".py", ".pdf")
plt.savefig(nombre_grafica, transparent='True', bbox_inches='tight')
# #
# nombre_grafica_2 = os.path.basename(__file__).replace(".py", ".pgf")
# plt.savefig(nombre_grafica_2, format='pgf')
#

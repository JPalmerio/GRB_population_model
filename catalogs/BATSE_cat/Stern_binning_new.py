import sys
import platform
if platform.system() == 'Linux':
    sys.path.insert(0,'/nethome/palmerio/Dropbox/Plotting_GUI/Src')
elif platform.system() == 'Darwin': 
    sys.path.insert(0,'/Users/palmerio/Dropbox/Plotting_GUI/Src')

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import plotting_functions as pf
from pathlib import Path
import seaborn as sns
root_dir = Path('/Users/palmerio/Science_projects/')


def efficiency_correction_Stern(pflx, c_e0=0.097, nu=2.34, norm=0.7):
    """
        The efficiency function of BATSE for detecting GRBs as a function of peak flux, derived by Stern+01
        c_e0 is in [counts/s/cm2]
        pflx is in [ph/s/cm2]
    """
    c_e = pflx * 0.75  # the conversion factor from counts to pflx comes from the Stern+01 paper as well, figure 7.
    return norm * (1.0 - np.exp(-(c_e/c_e0)**2) )**nu


T_min = 2.05   # secondes
P23_min = 0.01  # ph/cm2/s
N_bin_min = 1  # Nb of bins with flx >= 50% of pflx -> more robust than T90>2s because of uncertainties on T90 measurement (see Stern+01)
verbose = 1
file_Stern = str(root_dir/'GRB_population_model/observational_constraints/lognlogp.stern.dat')
text_size = 22
font = {'fontname':'Serif', 'size':text_size}


#################### Stern #######################
filename_S = 'Stern_cat.txt'
name_S = pf.read_column(filename_S, 0, dtype=str, splitter='\t')
T90_S = pf.read_column(filename_S, 7, splitter='\t')
N_bin_S = pf.read_column(filename_S, 8, splitter='\t')
peak_flux_1024_S = pf.read_column(filename_S, 3, splitter='\t')
final_long_peak_flux_1024_S = peak_flux_1024_S[N_bin_S > N_bin_min]
print("Summary for data read from Stern catalog:")
print(f"filename: {filename_S}")
print(f"N_GRB total: {len(name_S)}")
print(f"N_GRB with N_bin >= 2 (i.e. long GRBs): {len(final_long_peak_flux_1024_S)}")
print("pflx min, median, max (ph/s/cm2):"
      + f" {final_long_peak_flux_1024_S.min()} {np.median(final_long_peak_flux_1024_S)}, {final_long_peak_flux_1024_S.max()}")


# Stern
bins = pf.read_column(file_Stern, 0, array=False)
hist_Stern = 10**pf.read_column(file_Stern, 1)
hist_Stern_err = 10**pf.read_column(file_Stern, 2)
bins.append(1.8)
bins = 10**np.asarray(bins)/0.75
hist_S, bins_S = np.histogram(final_long_peak_flux_1024_S, bins=bins)
hist_S = np.asarray(hist_S, dtype=float)
hist_S_err = np.sqrt(hist_S)

bins2 = []
for i in range(len(bins)-8):
    bins2.append(bins[i])
bins2.append(16.)
bins2.append(20.)
bins2.append(28.)
bins2.append(50.)
bins2 = np.asarray(bins2, dtype=float)

hist_S_v2, __unused = np.histogram(final_long_peak_flux_1024_S, bins=bins2)

ln_oi = 0.
for i, val in enumerate(hist_S_v2):
    ln_oi += val*np.log(val) - val
print(f"ln(o_i!) = {ln_oi} before efficiency correction")

hist_S_v2 = np.asarray(hist_S_v2, dtype=float)
hist_S_v2_err = np.sqrt(hist_S_v2)

plt.hist(np.log10(final_long_peak_flux_1024_S), bins=np.log10(bins2),  # density=True,
         label='Original raw data from Stern+01 with new binning', histtype='step', linewidth=1.5, color='C1', )
plt.hist(np.log10(final_long_peak_flux_1024_S), bins=np.log10(bins),  # density=True,
         label='Original raw data from Stern+01 with old binning', histtype='step', linewidth=1.5, color='C2', )
# sns.kdeplot(np.log10(final_long_peak_flux_1024_S), color='k', ax=plt.gca(), label='KDE')
plt.gca().axhline(10, label='Minimum of 10 objects per bin')
plt.legend()

plt.gca().set_ylabel('N')
plt.gca().set_xlabel('log pflx 50-300 keV (ph/s/cm2)')
plt.gca().set_yscale('log')
plt.gca().set_ylim(0,400)
plt.show()

write_to_file = True
if write_to_file:
    outputfile_Stern_rebinned = 'Stern_lognlogp_rebinned_latest.txt'
    with open(outputfile_Stern_rebinned, 'w') as f:
        f.write("# Histogram of peak fluxes from Stern+01 catalog, rebinned to avoid low number (<10) bins.\n")
        f.write("# We removed the smallest bin as it was plagued with significant uncertainties.\n")
        f.write("# Peak count rate was NOT converted to peak flux since we used the raw data from Stern+01.\n")
        f.write("# The last bin ends at {}, and there are 26 bins.\n".format(bins2[-1]))
        f.write("# Nb_GRBs/year/log10(pflx) is NOT corrected for efficiency.\n")
        f.write("# 1. Peak flux (left bin edge) [ph/cm2/s]  2. Nb_GRBs  3. Error = sqrt(N) \n")
        for i in range(len(hist_S_v2)):
            if i == 0:
                pass
            else:
                f.write("{:e} {:e} {:e}\n".format(bins2[i], hist_S_v2[i], np.sqrt(hist_S_v2[i])))


exit()
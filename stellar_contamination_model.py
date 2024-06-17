import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits, ascii
import spectres
import sys

'''
This code uses PHOENIX models to assess the correction factor that must be applied to a transmission spectrum due to unocculted cold spots. The inputs are:

    1) percent photometric variability (usually want to model a few)
    2) spot temperatures (usually want to model a few)

The correction to the transit depth is given by:

delta_d/d=delta_f_lambda_0*(1-F_spot/F_photosphere)/(1-F_spot_lambda_0/F_photosphere_lambda_0)        (1)

where delta_f_lambda_0 is the measured or estimated variability at a given wavelength (lambda_0), F_spot is the flux of the spot, F_photosphere is the
flux of the photosphere, and the final term on the right-hand side is a normalization factor to the lambda_0 considered.

See Sing et al. 2011 for more details. 

The outputs of this code are:
1) Reduced chi-squared values for your inputs, informing you if the fit is good or not.
2) A plot of your data+stellar contamination models. 

'''


#-----------------------------------Inputs----------------------------------------------------------------

'''
This section contains all the inputs you will need to change. Be sure to read through it! 
'''
#Path to save figure
savepath='/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/spot_contamination_model/tests/commented_document_test/'

#Get PHOENIX models loaded in - I used https://www.stsci.edu/hst/instrumentation/reference-data-for-calibration-and-tools/astronomical-catalogs/phoenix-models-available-in-synphot
#For most MS stars, assume [M/H]=0.0 and log g=5.0. 
#For a photosphere=3300 K:
spectra_teff_3300=fits.open('/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/spot_contamination_model/phoenix_models/phoenixm00_3300.fits')
#Test spots between 2600-3200 K: 
spectra_teff_3000=fits.open('/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/spot_contamination_model/phoenix_models/phoenixm00_3000.fits')
spectra_teff_2800=fits.open('/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/spot_contamination_model/phoenix_models/phoenixm00_2800.fits')
spectra_teff_2600=fits.open('/Users/katiebennett/Documents/LTT1445Ab/UVIS_Zafar/spot_contamination_model/phoenix_models/phoenixm00_2600.fits')


#Pull in transmission spectrum data in .txt form with headers 'wvl', 'wvl_width', 'depth', and 'err'. 
data_all=ascii.read('/Users/katiebennett/Documents/LTT1445Ab/results/LTT1445Ab_spectra_combined_final_apr24_poseidon_w_header.txt')
wvl_all=data_all['wvl']
wvl_width_all=data_all['wvl_width']
depth_all=data_all['depth']*1e6 #convert to ppm
err_all=data_all['err']*1e6


#Pick a couple of delta_f_lambda_0 to try:
delta_f_lo=0.03 
delta_f_hi=0.1

#Temperatures you are probing: 
T_phot=3300
T_spot1=3000
T_spot2=2800
T_spot3=2600

#Based on the PHOENIX model you use, you will need to find the index corresponding to your desired lambda_0. To do this, we will find the closest
#wavelength value to lambda_0, but the width of wavelengths you search will depend on the resolution of your grid. 
#For example, for a fine grid, if we want to set lambda_0=8000 Angstroms, search around:
wvl_lo=7999.9 
wvl_hi=8000.1 
#If, however, your grid is coarser, you might need to set wvl_lo=7999 and wvl_hi=8001. Only one index should be found - the code will crash if 
#more than one is found, and then you know you need to adjust these search parameters. 



d=1914 #Mean transit depth in ppm (need to multiply right-side side of equation (1) by d, the mean transit depth, to get the transit correction factor in ppm)
k=3 #Number of free parameters in PHOENIX models (usually Teff, log g, and [M/H])


'''
Other things you MAY need to change in the code:

1) If you want to test more spot temperatures, you will need to add those lines accordingly. (Can ignore the fact that the variables I label correspond to certain temperatures.)
----> I will eventually make this code more automated so it does not matter your initial input number of spots. 
2) Change wavelength resolution in line 144 if your model grids are much different resolution than mine
'''



#-----------------------------------Set-up------------------------------------------------------------

#Get wavelengths from FITS file in Angstroms
wvls=spectra_teff_3300[1].data['WAVELENGTH']

#Pull intensity data from FITS file - may need to adjust this depending on format of FITS file you are using
flx_teff_3300=spectra_teff_3300[1].data['g50']
flx_teff_3000=spectra_teff_3000[1].data['g50']
flx_teff_2800=spectra_teff_2800[1].data['g50']
flx_teff_2600=spectra_teff_2600[1].data['g50']

'''
#UNCOMMENT TO PLOT HERE
#Double check that stellar models look good (you may want to )

plt.plot(wvls, flx_teff_3000, color='goldenrod', label='3000K')
plt.plot(wvls, flx_teff_2800, color='green', label='2800K')
plt.plot(wvls, flx_teff_2600, color='blue', label='2600K')
plt.xlim(2000, 20000)
plt.legend()
plt.show()
'''

#Get flux of spot model divided by flux of photosphere
F_3000_F_3300=flx_teff_3000/flx_teff_3300
F_2800_F_3300=flx_teff_2800/flx_teff_3300
F_2600_F_3300=flx_teff_2600/flx_teff_3300

#Determine index that corresponds to lambda_0
index=[]
for i in range(len(wvls)):
    if wvls[i]<wvl_hi and wvls[i]>wvl_lo:
        index.append(i)
#print('Index found closest lambda_0 is:', wvls[index])

#Calculate spot correction using Equation 4 from Sing 2011
delta_d_d_3000_lo=delta_f_lo*(1-F_3000_F_3300)/(1-flx_teff_3000[index]/flx_teff_3300[index])
delta_d_d_2800_lo=delta_f_lo*(1-F_2800_F_3300)/(1-flx_teff_2800[index]/flx_teff_3300[index])
delta_d_d_2600_lo=delta_f_lo*(1-F_2600_F_3300)/(1-flx_teff_2600[index]/flx_teff_3300[index])

delta_d_d_3000_hi=delta_f_hi*(1-F_3000_F_3300)/(1-flx_teff_3000[index]/flx_teff_3300[index])
delta_d_d_2800_hi=delta_f_hi*(1-F_2800_F_3300)/(1-flx_teff_2800[index]/flx_teff_3300[index])
delta_d_d_2600_hi=delta_f_hi*(1-F_2600_F_3300)/(1-flx_teff_2600[index]/flx_teff_3300[index])

#Bin up spectrum - how much you need to bin will depend on how fine your grid is! 
wvl_even=[]
for i in range(len(wvls)):
    if i % 2:
        wvl_even.append(wvls[i])
wvl_even=np.array(wvl_even)

wvl_even2=[]
for i in range(len(wvl_even)):
    if i % 2:
        wvl_even2.append(wvl_even[i])
wvl_even2=np.array(wvl_even2)

wvl_even3=[]
for i in range(len(wvl_even2)):
    if i % 2:
        wvl_even3.append(wvl_even2[i])
wvl_even3=np.array(wvl_even3)

wvl_even4=[]
for i in range(len(wvl_even3)):
    if i % 2:
        wvl_even4.append(wvl_even3[i])
wvl_even4=np.array(wvl_even4)

wvl_even5=[]
for i in range(len(wvl_even4)):
    if i % 2:
        wvl_even5.append(wvl_even4[i])
wvl_even5=np.array(wvl_even5)

wvl_even6=[]
for i in range(len(wvl_even5)):
    if i % 2:
        wvl_even6.append(wvl_even5[i])
wvl_even6=np.array(wvl_even6)

wvl_even7=[]
for i in range(len(wvl_even6)):
    if i % 2:
        wvl_even7.append(wvl_even6[i])
wvl_even7=np.array(wvl_even7)

#Bin up spectrum to a reasonable resolution  
#wvl_even5 is coarser than wvl_even4, which is coarser than wvl_even3, etc. 
delta_d_3000_bin_lo=spectres.spectres(wvl_even4, wvls, delta_d_d_3000_lo*d)
delta_d_2800_bin_lo=spectres.spectres(wvl_even4, wvls, delta_d_d_2800_lo*d)
delta_d_2600_bin_lo=spectres.spectres(wvl_even4, wvls, delta_d_d_2600_lo*d)

delta_d_3000_bin_hi=spectres.spectres(wvl_even4, wvls, delta_d_d_3000_hi*d)
delta_d_2800_bin_hi=spectres.spectres(wvl_even4, wvls, delta_d_d_2800_hi*d)
delta_d_2600_bin_hi=spectres.spectres(wvl_even4, wvls, delta_d_d_2600_hi*d)



#-----------------Calculate statistics for different models---------------------------------

def chi2v(obs, mod, err, k):
    answer=np.sum((obs-mod)**2/err**2)
    answer_red=answer/(len(obs)-k)
    print('good fit is between:', 1-np.sqrt(2/(len(obs)-k)), 'and', 1+np.sqrt(2/(len(obs)-k)))
    print('chi2v=', answer_red)
    if answer_red <=1+np.sqrt(2/(len(obs)-k)) and answer_red >=1-np.sqrt(2/(len(obs)-k)):
        print('chi2v indicates a good fit!')
    else: 
        print('chi2v does not show a good fit - rule this model out')
    return answer_red

#Apply an offset to get correction spectrum to line up with data at level of lambda_0
new_index=[]
for i in range(len(wvl_even4)):
    if wvl_even4[i]<wvl_hi+1 and wvl_even4[i]>wvl_lo-1:
        new_index.append(i)
#print('New index found closest lambda_0 is:', wvl_even4[new_index])
diff_hi=d-delta_d_2800_bin_hi[new_index]
diff_lo=d-delta_d_2800_bin_lo[new_index]

#Apply this offset and bin the spectrum to the data to be able to calculate chi2_nu
delta_d_3000_hi_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_3000_bin_hi+diff_hi))
delta_d_2800_hi_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_2800_bin_hi+diff_hi))
delta_d_2600_hi_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_2600_bin_hi+diff_hi))

delta_d_3000_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_3000_bin_lo+diff_lo))
delta_d_2800_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_2800_bin_lo+diff_lo))
delta_d_2600_bin_to_data=spectres.spectres(wvl_all, wvl_even4/1e4, (delta_d_2600_bin_lo+diff_lo))

print('Statistics for High Variability scenario:')
print('Chi2v for', T_spot1, 'K spots:')
chi2v_3000=chi2v(depth_all, delta_d_3000_hi_bin_to_data, err_all, k)
print('Chi2v for', T_spot2, 'K spots:')
chi2v_2800=chi2v(depth_all, delta_d_2800_hi_bin_to_data, err_all, k)
print('Chi2v for', T_spot3, 'K spots:')
chi2v_2600=chi2v(depth_all, delta_d_2600_hi_bin_to_data, err_all, k)
print()
print('Statistics for Low Variability Scenario:')
print('Chi2v for', T_spot1, 'K spots:')
chi2v_3000=chi2v(depth_all, delta_d_3000_bin_to_data, err_all, k)
print('Chi2v for', T_spot2, 'K spots:')
chi2v_2800=chi2v(depth_all, delta_d_2800_bin_to_data, err_all, k)
print('Chi2v for', T_spot3, 'K spots:')
chi2v_2600=chi2v(depth_all, delta_d_2600_bin_to_data, err_all, k)



#----------------------------------------Plot------------------------------------------------


plt.figure(figsize=[10,6])

plt.errorbar(wvl_all, depth_all, xerr=wvl_width_all, yerr=err_all, fmt='.', alpha=1, color='navy')

plt.plot(wvl_even4/1e4, delta_d_3000_bin_hi+diff_hi, alpha=0.5, color='red', label=r'$T_{\rm spots}=$'+str(T_spot1)+' K at 10% variability')
plt.plot(wvl_even4/1e4, delta_d_3000_bin_lo+diff_lo, color='red', alpha=0.5, ls='dotted', label=r'$T_{\rm spots}=$'+str(T_spot1)+' K at 2.5% variability')

plt.plot(wvl_even4/1e4, delta_d_2800_bin_hi+diff_hi, alpha=0.5, color='darkgreen', label=r'$T_{\rm spots}=$'+str(T_spot2)+' K')
plt.plot(wvl_even4/1e4, delta_d_2600_bin_hi+diff_hi, alpha=0.5, color='deepskyblue', label=r'$T_{\rm spots}=$'+str(T_spot3)+' K')

plt.plot(wvl_even4/1e4, delta_d_2800_bin_lo+diff_lo, color='darkgreen', alpha=0.5, ls='dotted')
plt.plot(wvl_even4/1e4, delta_d_2600_bin_lo+diff_lo, color='deepskyblue', alpha=0.5, ls='dotted')
plt.axhline(y=1907, ls='--', color='grey')

plt.legend()
plt.xlabel(r'Wavelength ($\rm \mu m$)', fontsize=14)
plt.ylabel(r'Transit Depth (ppm)', fontsize=14)
plt.xlim((0.3, 1.7))
#plt.ylim((1350, 3100))
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.ylim((1550, 2300))
plt.legend()
plt.savefig(savepath+'stellar_contamination_fwd_model.pdf')
plt.show()

#--------------------------------------------------------------------------------------------

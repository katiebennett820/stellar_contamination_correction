# stellar_contamination_correction
This code uses PHOENIX models to assess the correction factor that must be applied to a transmission spectrum due to unocculted cold spots. It is based on the work done in Sing et al. 2011. The inputs are:


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

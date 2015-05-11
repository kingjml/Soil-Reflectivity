soilReflectivity<-function(freq, theta, tSnowK, epiSnow, tSoilK, mvSoil, sigSoil, epiSoil){
	c = 299792458;
	theta_r = theta * (pi / 180);

	#Addapted from DMRT-ML
	costheta0 = cos(theta_r);
	n=sqrt(epiSoil/epiSnow);
	b = 1-(1-costheta0^2)/(n^2);
	costheta1=sqrt(b);
	fresnel_v=abs( (n*costheta0-costheta1)/(n*costheta0+costheta1) )^2;
	fresnel_h=abs( (costheta0-n*costheta1)/(costheta0+n*costheta1) )^2;

	effTheta = asin(sin(theta_r)/sqrt(epiSnow)) #Effective theta in snow
 	effWave = (c/sqrt(epiSnow))/freq*1e-7 #Effective wavelength (cm) in snow
 	effK = 2*pi/effWave #Effective wavenumber (rad/cm) in snow

      kSigma = effK * sigSoil

  	if (kSigma<0.07|kSigma>27.5){
		cat("kSigma out of range!")
	}

	r_h_mod = fresnel_h*exp(-kSigma^(sqrt(0.1*cos(effTheta))))
  	r_v_mod = r_h_mod*cos(effTheta )^0.655 #TODO Build in beta tuning function (0.655 currently)

      return(list(r_h_mod,r_v_mod))
}
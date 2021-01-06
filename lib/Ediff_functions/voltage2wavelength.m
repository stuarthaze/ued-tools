function eWavl = voltage2wavelength(volts)
h = 6.626069E-34;
e = 1.602176E-19;
me = 9.10938E-31;
c = 2.99792458E8;
term1 = h/sqrt(2*me*e*volts);
term2 = 1/sqrt(1 + (e*volts/(2*me*c*c)));
eWavl = term1*term2;
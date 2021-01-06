function v = volt2vel(volt)
   e = 1.602176e-19;
   me = 9.10938e-31;
   c = 2.997925e8;
   v = sqrt(1-(me./(me + e*volt./c^2)).^2);
return

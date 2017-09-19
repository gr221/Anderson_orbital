function a = fermi(e,Temp)

%a = 1./(1+exp(e/T));
a = 1/2*(1-tanh(e./(2*Temp)));
end
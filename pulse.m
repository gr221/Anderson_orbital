function v = pulse(t)

global e_d U Gamma  T mu0 w D A t0 t_del

%v = A*sin(2*pi*w*(t-t0+3)).*exp(-(t-t0).^2/D^2);
%v = A*cos(2*pi*w*(t-t0)).*exp(-(t-t0).^2/D^2);
%v = A/2*(1-fermi(t-t0,D));
v = A*exp(-(t-t0).^2/D^2);
%v = zeros(size(t));

%v = A*(exp(-(t-t0).^2/D^2) - exp(-(t-t0-t_del).^2/D^2));
end
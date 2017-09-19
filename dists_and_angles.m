function [d_ben, d_io, d_No, d_Ni, d_CH, d_NH, d_e, d_NM, r_p, phi, psi, gamma, gamma2] = dists_and_angles
%
% [d_ben, d_io, d_No, d_Ni, d_CH, d_NH, d_e, d_NM, r_p, phi, psi, gamma, gamma2] = dists_and_angles
% 
% 20.02.2015: ALL values taken from Liao et al. J. Chem. Phys. 114, 9780 (2001)
% 
% 
% 
% Pc Unit Cell
%
%	    C12
%	  /		 \	
%	 /		  \
%	/		   \
%	N11			N10
%	|			  \
% 	|			   \        H7
% 	|				\       |
% 	|				C9     C7     H6
% 	|			   /  \  /    \  /
% 	|			  /    C8      C6
%	0------------N1     |       |
% 				  \    C3      C5
% 				   \  /  \    /  \
% 					C2     C4     H5
% 							|
% 						   H4
% 
% 
% phi     = 109.1;				% Winkel C9-N1-C2
% psi		= 73.6;					% Winkel x-Achse - C9-N10
% gamma	= 16.7;					% Winkel x-Achse - C9-C8
% gamma2	= 90-psi;				% Winkel x-Achse N10-C12 (next unit cell)
% 
% d_ben = 1.41;		% zwischen 2 C Atomen im äußeren benzene (a im PH)
% d_io  = 1.46;		% zwischen C innerer Ring zu C im benzene (b im PH)
% d_No  = 1.33;		% zwischen C und außen liegendem N (c im PH)
% d_Ni  = 1.38;		% zwischen C und innen liegendem N (d im PH)
% d_CH  = 1.09;		% Abstand C - H (R_C-H in LiS)
% d_NH  = 1.0;		% Abstand N - H bei H2Pc (gemittelt) (aus DaW)
% d_NM  = 1.98;		% Abstand Metall - Stickstoff
% d_e   = 2*d_NM;		% Abstand von Atom 1 - Atom 21, 3.65 < e < 4.50()
% r_p   = 5.29;		% p-state radius Stickstoff



d_ben = 1.41;		% zwischen 2 C Atomen im äußeren benzene (a im PH)
d_io  = 1.46;		% zwischen C innerer Ring zu C im benzene (b im PH)
d_No  = 1.33;		% zwischen C und außen liegendem N (c im PH)
d_Ni  = 1.38;		% zwischen C und innen liegendem N (d im PH)
d_CH  = 1.09;		% Abstand C - H (R_C-H in LiS)
d_NH  = 1.0;		% Abstand N - H bei H2Pc (gemittelt) (aus DaW)
d_NM  = 1.98;		% Abstand Metall - Stickstoff
d_e   = 2*d_NM;		% Abstand von Atom 1 - Atom 21, 3.65 < e < 4.50()
r_p   = 5.29;		% p-state radius Stickstoff

phi     = 109.1;				% Winkel C9-N1-C2
psi		= 73.6;					% Winkel x-Achse - C9-N10
gamma	= 16.7;					% Winkel x-Achse - C9-C8
gamma2	= 90-psi;				% Winkel x-Achse N10-C12 (next unit cell)

















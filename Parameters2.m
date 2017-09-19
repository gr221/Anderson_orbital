function [eps_H1s, eps_C2s, eps_C2p, eps_N2s, eps_N2p, eps_M4s, eps_M3d,r_d, Valence_electrons] = Parameters2(M)
% [eps_H1s, eps_C2s, eps_C2p, eps_N2s, eps_N2p, eps_M4s, eps_M3d,r_d, Valence_electrons] = Parameters2(M)
%

% the values for the energies are taken from the HF calculations from Mann

Ry = -13.605693;

eps_H1s = Ry;					%[eV] H 1s 
eps_C2s = 1.424123*Ry;			%[eV] C 2s 
eps_C2p = .8137992*Ry;			%[eV] C 2p 
eps_N2s = 1.927336*Ry;			%[eV] N 2s 
eps_N2p = 1.017308*Ry;			%[eV] N 2p

switch M
	case 'Mn'
    eps_M4s = .5026975*Ry;		%[eV]  Mn 4s
    eps_M3d = 1.122340*Ry;		%[eV]  Mn 3d
    r_d=0.799;					% Mn d-state radius
    Valence_electrons = 191;	%Mn

	case 'Fe'
    eps_M4s = .5202591*Ry;		%[eV]  Fe 4s 
    eps_M3d = 1.215704*Ry;		%[eV]  Fe 3d
    r_d=0.744;					% Fe d-state radius 
    Valence_electrons = 192;	%Fe
	
	case 'Co'
    eps_M4s = .5371718*Ry;		%[eV]  Co 4s 
    eps_M3d = 1.306100*Ry;		%[eV]  Co 3d
    r_d=0.696;					% Co d-state radius  
    Valence_electrons = 193;	%Co

	case 'Ni'
    eps_M4s = .5535508*Ry;		%[eV]  Ni 4s 
    eps_M3d = 1.394161*Ry;		%[eV]  Ni 3d
    r_d=0.652;					% Ni d-state radius   
    Valence_electrons = 194;	%Ni

	case 'Cu'
    eps_M4s = .4769883*Ry;		%[eV] Cu 4s -6.49
    eps_M3d = .9824537*Ry;		%[eV] Cu 3d -20.26
    r_d    = 0.688;				% Cu d-state radius   0.688  
    Valence_electrons = 195;	%Cu
end


end
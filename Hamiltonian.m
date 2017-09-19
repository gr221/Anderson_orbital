function [H_sigma,H_pi]=Hamiltonian(Metal) 
% [H_sigma,H_pi]=Hamiltonian(Metal)
%
% Calculates sigma and pi Hamiltonian of MPc
% Metal = Cu Mn Fe Co Ni
% 
% date: 02.02.2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[E_H_1s,E_C_2s,E_C_2p,E_N_2s,E_N_2p,E_M_4s,E_M_3d,r_d,~] = Parameters2(Metal);

splitting = zeros(5,1);       % ligand field induced splitting of metallic d levels
Delta0 = 0;                   % mean splitting fur Cu2+

splitting(1) = -0.43 * Delta0;           % z^2
splitting(2) =  1.23 * Delta0;           % x^2-y^2
splitting(3) =  0.23 * Delta0;           % xy
splitting(4) = -0.51 * Delta0;           % xz
splitting(5) = -0.51 * Delta0;           % yz

%% hopping Parameter eta_( s/p , s/p , [s]igma/[p]i )
%   aus W.A. Harrison - New Tight Binding Parameters...

eta_sss = -1.32; % -1.32
eta_sps =  1.42; % -1.42
eta_pps =  2.22; %  2.22
eta_ppp = -0.63; % -0.63
eta_sds = -3.16; % -3.16
eta_pds = -3*sqrt(15)/(2*pi);	% -3*sqrt(15)/(2*pi) from ??? , -2.95 from Pure&Appl Chem. 61, 2161 (1989), arxiv:cond-mat/0405067
eta_pdp = 3*sqrt(5)/(2*pi);		% 3*sqrt(5)/(2*pi) from ???, 1.36 from arxiv:cond-mat/0405067

%% lengths (in Angstrom) and angles 

[d_ben, d_io, d_No, d_Ni, d_CH, ~, ~, d_NM, r_p, phi, psi, gamma, gamma2]   = dists_and_angles;

%% interatomic hopping integrals V_(...)_( s/p , s/p , [s]igma/[p]i )
%   Faktor W = 7.62 kommt von hbar^2/m_e, Einheit eV/Angstrom^2

W=7.62;

% hopping zwischen 2 C Atomen im äußeren benzene
V_ben_sss = W*eta_sss/d_ben^2;
V_ben_sps = W*eta_sps/d_ben^2;
V_ben_pps = W*eta_pps/d_ben^2;
V_ben_ppp = W*eta_ppp/d_ben^2;

% hopping zwischen C innerer Ring zu C im benzene
V_io_sss  = W*eta_sss/d_io^2;
V_io_sps  = W*eta_sps/d_io^2;
V_io_pps  = W*eta_pps/d_io^2;
V_io_ppp  = W*eta_ppp/d_io^2;

% hopping zwischen C und außen liegendem N
V_No_sss  = W*eta_sss/d_No^2;
V_No_sps  = W*eta_sps/d_No^2;
V_No_pps  = W*eta_pps/d_No^2;
V_No_ppp  = W*eta_ppp/d_No^2;

% hopping zwischen C und innen liegendem N
V_Ni_sss  = W*eta_sss/d_Ni^2;
V_Ni_sps  = W*eta_sps/d_Ni^2;
V_Ni_pps  = W*eta_pps/d_Ni^2;
V_Ni_ppp  = W*eta_ppp/d_Ni^2;

% hopping zwischen C und H im benzene
V_CH_sss  = W*eta_sss/d_CH^2;
V_CH_sps  = W*eta_sps/d_CH^2;

% hopping zwischen N und H bei H2Pc
% % V_NH_sss  = W*eta_sss/d_NH^2;
% % V_NH_sps  = W*eta_sps/d_NH^2;

% hopping zwischen Cu und N bei MPc
V_NM_sss  = W*eta_sss/d_NM^2;
V_NM_sps  = W*eta_sps/d_NM^2;
V_NM_sds  = W*eta_sds*sqrt(r_d^3/d_NM^7);
V_NM_pds  = W*eta_pds*sqrt(r_p*r_d^3)/d_NM^4;
V_NM_pdp  = W*eta_pdp*sqrt(r_p*r_d^3)/d_NM^4;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Sigma Hamiltonian

%% on site energies %%%

H_us = zeros(34);

% Atom 1 (N)
H_us(1,1)   = E_N_2s;
H_us(2,2)   = E_N_2p;
H_us(3,3)   = E_N_2p;

% Atom 2 (C)
H_us(4,4)   = E_C_2s;
H_us(5,5)   = E_C_2p;
H_us(6,6)   = E_C_2p;
         
% Atom 3 (C)
H_us(7,7)   = E_C_2s;
H_us(8,8)   = E_C_2p;
H_us(9,9)   = E_C_2p;

% Atom 4 (C)
H_us(10,10) = E_C_2s;
H_us(11,11) = E_C_2p;
H_us(12,12) = E_C_2p;

% H Atom 4
H_us(13,13) = E_H_1s;

% Atom 5 (C)
H_us(14,14) = E_C_2s;
H_us(15,15) = E_C_2p;
H_us(16,16) = E_C_2p;

% H Atom 5
H_us(17,17) = E_H_1s;

% Atom 6 (C)
H_us(18,18) = E_C_2s;
H_us(19,19) = E_C_2p;
H_us(20,20) = E_C_2p;

% H Atom 6
H_us(21,21) = E_H_1s;

% Atom 7 (C)
H_us(22,22) = E_C_2s;
H_us(23,23) = E_C_2p;
H_us(24,24) = E_C_2p;

% H Atom 7
H_us(25,25) = E_H_1s;

% Atom 8 (C)
H_us(26,26) = E_C_2s;
H_us(27,27) = E_C_2p;
H_us(28,28) = E_C_2p;

% Atom 9 (C)
H_us(29,29) = E_C_2s;
H_us(30,30) = E_C_2p;
H_us(31,31) = E_C_2p;

% Atom 10 (N)
H_us(32,32) = E_N_2s;
H_us(33,33) = E_N_2p;
H_us(34,34) = E_N_2p;

%% Hopping Matrix V_us(i,j) für Sigma - System %%%
%     zunächst nur für i<j 
%     => komplette Matrix ist (V_us)+(V_us)' (transponiert)

%% Atom 1 (N) (Indizes 1,2,3) entsprechen (s, p_x, p_y) Orbitalen

V_us = zeros(34);

V_us(1,4)   =   V_Ni_sss;                     
V_us(1,5)   =   V_Ni_sps*cosd(phi/2);
V_us(1,6)   = - V_Ni_sps*cosd(90-phi/2);

V_us(2,4)   = - V_Ni_sps*cosd(phi/2);
V_us(2,5)   =   V_Ni_pps*cosd(phi/2)^2 + V_Ni_ppp*sind(phi/2)^2;
V_us(2,6)   = -(V_Ni_pps - V_Ni_ppp)*cosd(phi/2)*sind(phi/2);

V_us(3,4)   =   V_Ni_sps*cosd(90-phi/2);
V_us(3,5)   = -(V_Ni_pps - V_Ni_ppp)*cosd(phi/2)*sind(phi/2);
V_us(3,6)   =   V_Ni_pps*cosd(90-phi/2)^2 + V_Ni_ppp*sind(90-phi/2)^2;

V_us(1,29)  =   V_Ni_sss;
V_us(1,30)  =   V_Ni_sps*cosd(phi/2);
V_us(1,31)  =   V_Ni_sps*cosd(90-phi/2);

V_us(2,29)  =  -V_Ni_sps*cosd(phi/2);
V_us(2,30)  =   V_Ni_pps*cosd(phi/2)^2 + V_Ni_ppp*sind(phi/2)^2;
V_us(2,31)  =  (V_Ni_pps - V_Ni_ppp)*cosd(phi/2)*sind(phi/2);

V_us(3,29)  =  -V_Ni_sps*cosd(90-phi/2);
V_us(3,30)  =  (V_Ni_pps - V_Ni_ppp)*cosd(phi/2)*sind(phi/2);
V_us(3,31)  =   V_Ni_pps*cosd(90-phi/2)^2 + V_Ni_ppp*sind(90-phi/2)^2;

%% Atom 2 (C) (4,5,6)

V_us(4,7)   =   V_io_sss;
V_us(4,8)   =   V_io_sps*cosd(gamma);
V_us(4,9)   =   V_io_sps*cosd(90-gamma);

V_us(5,7)   =  -V_io_sps*cosd(gamma);
V_us(5,8)   =   V_io_pps*cosd(gamma)^2 + V_io_ppp*sind(gamma)^2;
V_us(5,9)   =  (V_io_pps - V_io_ppp)*sind(gamma)*cosd(gamma);

V_us(6,7)   =  -V_io_sps*cosd(90-gamma);
V_us(6,8)   =  (V_io_pps - V_io_ppp)*sind(gamma)*cosd(gamma);
V_us(6,9)   =   V_io_pps*cosd(90-gamma)^2 + V_io_ppp*sind(90-gamma)^2;

% V_us(4,7)   =   V_io_sss;
% V_us(4,8)   =   V_io_sps*cosd(15);
% V_us(4,9)   =   V_io_sps*cosd(75);
% 
% V_us(5,7)   =  -V_io_sps*cosd(15);
% V_us(5,8)   =   V_io_pps*cosd(15)^2 + V_io_ppp*sind(15)^2;
% V_us(5,9)   =  (V_io_pps - V_io_ppp)*sind(15)*cosd(15);
% 
% V_us(6,7)   =  -V_io_sps*cosd(75);
% V_us(6,8)   =  (V_io_pps - V_io_ppp)*sind(15)*cosd(15);
% V_us(6,9)   =   V_io_pps*cosd(75)^2 + V_io_ppp*sind(75)^2;

%% Atom 3 (C) (7,8,9)

V_us(7,10)  =   V_ben_sss;
V_us(7,11)  =   V_ben_sps*cosd(30);
V_us(7,12)  =  -V_ben_sps*cosd(60);

V_us(8,10)  =  -V_ben_sps*cosd(30);
V_us(8,11)  =   V_ben_pps*cosd(30)^2 + V_ben_ppp*sind(30)^2;
V_us(8,12)  = -(V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);

V_us(9,10)  =   V_ben_sps*cosd(60);
V_us(9,11)  = -(V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);
V_us(9,12)  =   V_ben_pps*cosd(60)^2 + V_ben_ppp*sind(60)^2;

V_us(7,26)  =   V_ben_sss;
V_us(7,27)  =   0;
V_us(7,28)  =   V_ben_sps;

V_us(8,26)  =   0;
V_us(8,27)  =   V_ben_ppp;
V_us(8,28)  =   0;

V_us(9,26)  =  -V_ben_sps;
V_us(9,27)  =   0;
V_us(9,28)  =   V_ben_pps;

%% Atom 4 (C) (10,11,12) (H Atom 4: 13)

V_us(10,13) =   V_CH_sss;

V_us(10,14) =   V_ben_sss;
V_us(10,15) =   V_ben_sps*cosd(30);
V_us(10,16) =   V_ben_sps*cosd(60);

V_us(11,13) =   0;

V_us(11,14) =  -V_ben_sps*cosd(30);
V_us(11,15) =   V_ben_pps*cosd(30)^2 + V_ben_ppp*sind(30)^2;
V_us(11,16) =  (V_ben_pps - V_ben_ppp)*cosd(30)*sind(30); 

V_us(12,13) =   V_CH_sps;

V_us(12,14) =  -V_ben_sps*cosd(60);
V_us(12,15) =  (V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);
V_us(12,16) =   V_ben_pps*cosd(60)^2 + V_ben_ppp*sind(60)^2;

%% Atom 5 (C) (14,15,16) (H Atom 5: 17)

V_us(14,17) =   V_CH_sss;

V_us(14,18) =   V_ben_sss;
V_us(14,19) =   0;
V_us(14,20) =   V_ben_sps;

V_us(15,17) =  -V_CH_sps*cosd(30);

V_us(15,18) =   0;
V_us(15,19) =   V_ben_ppp;
V_us(15,20) =   0;

V_us(16,17) =   V_CH_sps*cosd(60);

V_us(16,18) =  -V_ben_sps;
V_us(16,19) =   0;
V_us(16,20) =   V_ben_pps;

%% Atom 6 (C) (18,19,20) (H Atom 6: 21)

V_us(18,21) =   V_CH_sss;

V_us(18,22) =   V_ben_sss;
V_us(18,23) =  -V_ben_sps*cosd(30);
V_us(18,24) =   V_ben_sps*cosd(60);

V_us(19,21) =  -V_CH_sps*cosd(30);

V_us(19,22) =   V_ben_sps*cosd(30);
V_us(19,23) =   V_ben_pps*cosd(30)^2 + V_ben_ppp*sind(30)^2;
V_us(19,24) = -(V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);

V_us(20,21) =  -V_CH_sps*cosd(60);

V_us(20,22) =  -V_ben_sps*cosd(60);
V_us(20,23) = -(V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);
V_us(20,24) =   V_ben_pps*cosd(60)^2 + V_ben_ppp*sind(60)^2;

%% Atom 7 (C) (22,23,24) (H Atom 7: 25)

V_us(22,25) =   V_CH_sss;

V_us(22,26) =   V_ben_sss;
V_us(22,27) =  -V_ben_sps*cosd(30);
V_us(22,28) =  -V_ben_sps*cosd(60);

V_us(23,25) =   0;

V_us(23,26) =   V_ben_sps*cosd(30);
V_us(23,27) =   V_ben_pps*cosd(30)^2 + V_ben_ppp*sind(30)^2;
V_us(23,28) =  (V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);

V_us(24,25) =  -V_CH_sps;

V_us(24,26) =   V_ben_sps*cosd(60);
V_us(24,27) =  (V_ben_pps - V_ben_ppp)*cosd(30)*sind(30);
V_us(24,28) =   V_ben_pps*cosd(60)^2 + V_ben_ppp*sind(60)^2;

%% Atom 8 (C) (26,27,28)

V_us(26,29) =   V_io_sss;
V_us(26,30) =  -V_io_sps*cosd(gamma);
V_us(26,31) =   V_io_sps*cosd(90-gamma);

V_us(27,29) =   V_io_sps*cosd(gamma);
V_us(27,30) =   V_io_pps*cosd(gamma)^2 + V_io_ppp*sind(gamma)^2;
V_us(27,31) = -(V_io_pps - V_io_ppp)*cosd(gamma)*sind(gamma);

V_us(28,29) =  -V_io_sps*cosd(90-gamma);
V_us(28,30) =  -(V_io_pps - V_io_ppp)*cosd(gamma)*sind(gamma);
V_us(28,31) =   V_io_pps*cosd(90-gamma)^2 + V_io_ppp*sind(90-gamma)^2;

% V_us(26,29) =   V_io_sss;
% V_us(26,30) =  -V_io_sps*cosd(15);
% V_us(26,31) =   V_io_sps*cosd(75);
% 
% V_us(27,29) =   V_io_sps*cosd(15);
% V_us(27,30) =   V_io_pps*cosd(15)^2 + V_io_ppp*sind(15)^2;
% V_us(27,31) = -(V_io_pps - V_io_ppp)*cosd(15)*sind(15);
% 
% V_us(28,29) =  -V_io_sps*cosd(75);
% V_us(28,30) =  -(V_io_pps - V_io_ppp)*cosd(15)*sind(15);
% V_us(28,31) =   V_io_pps*cosd(75)^2 + V_io_ppp*sind(75)^2;

%% Atom 9 (C) (29,30,31) [ Atom 10 (N) (32,33,34) ]

V_us(29,32) =   V_No_sss;
V_us(29,33) =  -V_No_sps*cosd(psi);
V_us(29,34) =   V_No_sps*cosd(90-psi);

V_us(30,32) =   V_No_sps*cosd(psi);
V_us(30,33) =   V_No_pps*cosd(psi)^2 + V_No_ppp*sind(psi)^2;
V_us(30,34) = -(V_No_pps - V_No_ppp)*cosd(psi)*sind(psi);

V_us(31,32) =  -V_No_sps*cosd(90-psi);
V_us(31,33) = -(V_No_pps - V_No_ppp)*cosd(psi)*sind(psi);
V_us(31,34) =   V_No_pps*cosd(90-psi)^2 + V_No_ppp*sind(90-psi)^2;
 
%% V_us zusammenbauen:

V_us = V_us + V_us';

%% Hamiltonian fürs Metall

H_M      = zeros(4);

H_M(1,1) = E_M_4s;                    % s orbital
H_M(2,2) = E_M_3d + splitting(1);     % z^2 orbital
H_M(3,3) = E_M_3d + splitting(2);     % x^2-y^2 orbital
H_M(4,4) = E_M_3d + splitting(3);     % xy orbital

%% hopping zwischen Metall und Unit Cell

% H_M1_u1      = zeros(4,34); % zwischen Cu und unit cell 1
% H_M1_u2      = zeros(4,34); % zwischen Cu und unit cell 2
% H_M1_u3      = zeros(4,34); % zwischen Cu und unit cell 3
% H_M1_u4      = zeros(4,34); % zwischen Cu und unit cell 4

H_M2_u1      = zeros(4,34); % zwischen Cu und unit cell 1
H_M2_u2      = zeros(4,34); % zwischen Cu und unit cell 2
H_M2_u3      = zeros(4,34); % zwischen Cu und unit cell 3
H_M2_u4      = zeros(4,34); % zwischen Cu und unit cell 4


%-----------------------------------------------------------%
%                                                           %
% --- !!! WICHTIG !!! NOCH ZU KLÄREN !!! ---                %
% wegen p_x oder p_y orbitalen in UC 2 und 4 und            %
% vorzeichen wegen verschiedener orientierungen !!!         %
%                                                           %
%-----------------------------------------------------------%


% 2 versionen: H_M1 (meine) und H_M2

%% s orbital mit s orbital

% H_M1_u1(1,1) = V_NM_sss;     % mit s orbital von N
% H_M1_u2(1,1) = V_NM_sss;     % mit s orbital von N
% H_M1_u3(1,1) = V_NM_sss;     % mit s orbital von N
% H_M1_u4(1,1) = V_NM_sss;     % mit s orbital von N


H_M2_u1(1,1) = V_NM_sss;     % mit s orbital von N
H_M2_u2(1,1) = V_NM_sss;     % mit s orbital von N
H_M2_u3(1,1) = V_NM_sss;     % mit s orbital von N
H_M2_u4(1,1) = V_NM_sss;     % mit s orbital von N

%% s orbital mit p orbital

% H_M1_u1(1,2) = V_NM_sps;     % mit p_x  
% 
% H_M1_u2(1,3) = V_NM_sps;     % mit_py  
% 
% H_M1_u3(1,2) = - V_NM_sps;   % mit p_x (minus wegen orientierung) 
%                             
% H_M1_u4(1,3) = - V_NM_sps;   % mit p_y (minus wegen orientierung)  
                            

H_M2_u1(1,2) = V_NM_sps;     % mit p_x  

H_M2_u2(1,2) = V_NM_sps;     % mit_py
                            
H_M2_u3(1,2) = V_NM_sps;    % mit p_x
                            
H_M2_u4(1,2) = V_NM_sps;    % mit p_y
                           
%% z^2 orbital mit s orbital

% H_M1_u1(2,1) = -.5*V_NM_sds;
% H_M1_u2(2,1) = -.5*V_NM_sds;
% H_M1_u3(2,1) = -.5*V_NM_sds;
% H_M1_u4(2,1) = -.5*V_NM_sds;

H_M2_u1(2,1) = -.5*V_NM_sds;
H_M2_u2(2,1) = -.5*V_NM_sds;
H_M2_u3(2,1) = -.5*V_NM_sds;
H_M2_u4(2,1) = -.5*V_NM_sds;

%% z^2 orbital mit p orbital
% 
% H_M1_u1(2,2) =  .5*V_NM_pds; % mit p_x 
% 
% H_M1_u2(2,3) =  .5*V_NM_pds; % mit p_y  
%                             
% H_M1_u3(2,2) = -.5*V_NM_pds; % mit p_x (minus wegen orientierung) 
%                             
% H_M1_u4(2,3) = -.5*V_NM_pds; % mit p_y (minus wegen orientierung)


H_M2_u1(2,2) =  .5*V_NM_pds; % mit p_x 

H_M2_u2(2,2) =  .5*V_NM_pds; % mit p_y
                            
H_M2_u3(2,2) =  .5*V_NM_pds; % mit p_x
                            
H_M2_u4(2,2) =  .5*V_NM_pds;  % mit p_y   
                                                                                 
%% x^2-y^2 orbital mit s orbital
% 
% H_M1_u1(3,1) =  sqrt(3)/2*V_NM_sds;
% 
% H_M1_u2(3,1) = -sqrt(3)/2*V_NM_sds;
% 
% H_M1_u3(3,1) =  sqrt(3)/2*V_NM_sds;
% 
% H_M1_u4(3,1) = -sqrt(3)/2*V_NM_sds;


H_M2_u1(3,1) =  sqrt(3)/2*V_NM_sds;

H_M2_u2(3,1) = -sqrt(3)/2*V_NM_sds;

H_M2_u3(3,1) =  sqrt(3)/2*V_NM_sds;

H_M2_u4(3,1) = -sqrt(3)/2*V_NM_sds;

%% x^2-y^2 orbital mit p orbital

% H_M1_u1(3,2) = -sqrt(3)/2*V_NM_pds;  % mit p_x
% 
% H_M1_u2(3,3) =  sqrt(3)/2*V_NM_pds;  % mit p_y
% 
% H_M1_u3(3,2) =  sqrt(3)/2*V_NM_pds;  % mit p_x 
%                                     
% H_M1_u4(3,3) = -sqrt (3)/2*V_NM_pds; % mit p_y  
                                    

H_M2_u1(3,2) = -sqrt(3)/2*V_NM_pds;  % mit p_x

H_M2_u2(3,2) =  sqrt(3)/2*V_NM_pds;  % mit p_y

H_M2_u3(3,2) = -sqrt(3)/2*V_NM_pds;  % mit p_x
                                    
H_M2_u4(3,2) =  sqrt(3)/2*V_NM_pds;  % mit p_y                                   
                                   
%% xy orbital mit p orbital

% H_M1_u1(4,3) = -V_NM_pdp;    % mit p_y
% 
% H_M1_u2(4,2) = -V_NM_pdp;    % mit p_x
% 
% H_M1_u3(4,3) =  V_NM_pdp;    % mit p_y  
%                             
% H_M1_u4(4,2) =  V_NM_pdp;    % mit p_x  
                            
                            
H_M2_u1(4,3) = -V_NM_pdp;    % mit p_y

H_M2_u2(4,3) = V_NM_pdp;    % mit p_x
                            
H_M2_u3(4,3) = -V_NM_pdp;    % mit p_y
                            
H_M2_u4(4,3) =  V_NM_pdp;      % mit p_x                          
                            
%%  hopping zwischen unit cells: V_uu - im Sigma - System: %%
%    hopping zwischen 10. und 2. Atom
%    WICHTIG: Atom 2 jetzt um 90 Grad verdreht => x <-> y
%    heißt: x zeigt nach oben, y nach links
%
%                   x
%                   ^
%                   |
%                   |
%          y <------o

V_uu = zeros(34);

V_uu(32,4) =   V_No_sss;                    
V_uu(32,5) =   V_No_sps*sind(gamma2);      
V_uu(32,6) =   V_No_sps*cosd(gamma2);        
 
V_uu(33,4) =   V_No_sps*cosd(gamma2);         
V_uu(33,5) = -(V_No_pps - V_No_ppp)*cosd(gamma2)*sind(gamma2);
V_uu(33,6) = -(V_No_pps*cosd(gamma2)^2 + V_No_ppp*sind(gamma2)^2);

V_uu(34,4) =  -V_No_sps*sind(gamma2);
V_uu(34,5) =   V_No_pps*sind(gamma2)^2 + V_No_ppp*cosd(gamma2)^2;
V_uu(34,6) =  (V_No_pps - V_No_ppp)*cosd(gamma2)*sind(gamma2);

% V_uu(32,4) =   V_No_sss;                    
% V_uu(32,5) =   V_No_sps*cosd(90-Gamma);      
% V_uu(32,6) =   V_No_sps*cosd(Gamma);        
%  
% V_uu(33,4) =   V_No_sps*cosd(Gamma);         
% V_uu(33,5) = -(V_No_pps - V_No_ppp)*cosd(Gamma)*sind(Gamma);
% V_uu(33,6) = -(V_No_pps*cosd(Gamma)^2 + V_No_ppp*sind(Gamma)^2);
% 
% V_uu(34,4) =  -V_No_sps*cosd(90-Gamma);
% V_uu(34,5) =   V_No_pps*cosd(90-Gamma)^2 + V_No_ppp*sind(90-Gamma)^2;
% V_uu(34,6) =  (V_No_pps - V_No_ppp)*cosd(90-Gamma)*sind(90-Gamma);

%% Sigma Hamiltonian zusammenbauen

% Unit Cell Hamiltonian (onsite energies + off-diagonal Elemente)

H_us = H_us + V_us;

%% kompletter Hamiltonian fürs Sigma System

% % H_sig1 = [  H_us        ,   V_uu        ,   zeros(34)   ,   V_uu'       ,   H_M1_u1'     ;   
% %             V_uu'       ,   H_us        ,   V_uu        ,   zeros(34)   ,   H_M1_u2'     ;    
% %             zeros(34)   ,   V_uu'       ,   H_us        ,   V_uu        ,   H_M1_u3'     ;   
% %             V_uu        ,   zeros(34)   ,   V_uu'       ,   H_us        ,   H_M1_u4'     ;
% %             H_M1_u1     ,   H_M1_u2     ,   H_M1_u3     ,   H_M1_u4     ,   H_M         ];
        

H_sig2 = [  H_us        ,   V_uu        ,   zeros(34)   ,   V_uu'       ,   H_M2_u1'     ;   
            V_uu'       ,   H_us        ,   V_uu        ,   zeros(34)   ,   H_M2_u2'     ;    
            zeros(34)   ,   V_uu'       ,   H_us        ,   V_uu        ,   H_M2_u3'     ;   
            V_uu        ,   zeros(34)   ,   V_uu'       ,   H_us        ,   H_M2_u4'     ;
            H_M2_u1     ,   H_M2_u2     ,   H_M2_u3     ,   H_M2_u4     ,   H_M         ];
        
%%

H_sigma = H_sig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pi - Hamiltonian

%% unit cell Hamiltonian ohne inneres H Atom
%   nummeriert nach sites

H_up = zeros(10);

H_up(1,1)   = E_N_2p;
H_up(1,2)   = V_Ni_ppp;
H_up(1,9)   = V_Ni_ppp;

H_up(2,1)   = V_Ni_ppp;
H_up(2,2)   = E_C_2p;
H_up(2,3)   = V_io_ppp;

H_up(3,2)   = V_io_ppp;
H_up(3,3)   = E_C_2p;
H_up(3,4)   = V_ben_ppp;
H_up(3,8)   = V_ben_ppp;

H_up(4,3)   = V_ben_ppp;
H_up(4,4)   = E_C_2p;
H_up(4,5)   = V_ben_ppp;

H_up(5,4)   = V_ben_ppp;
H_up(5,5)   = E_C_2p;
H_up(5,6)   = V_ben_ppp;

H_up(6,5)   = V_ben_ppp;
H_up(6,6)   = E_C_2p;
H_up(6,7)   = V_ben_ppp;

H_up(7,6)   = V_ben_ppp;
H_up(7,7)   = E_C_2p;
H_up(7,8)   = V_ben_ppp;

H_up(8,7)   = V_ben_ppp;
H_up(8,8)   = E_C_2p;
H_up(8,9)   = V_io_ppp;
H_up(8,3)   = V_ben_ppp;

H_up(9,8)   = V_io_ppp;
H_up(9,9)   = E_C_2p;
H_up(9,10)  = V_No_ppp;
H_up(9,1)   = V_Ni_ppp;

H_up(10,9)  = V_No_ppp;
H_up(10,10) = E_N_2p;

%% Hamiltonian für hopping zwischen unit cells
%   inter UC hopping immer zwischen Stickstoff an Platz 10
%   mit Kohlenstoff am Platz 2

H_up_hop = zeros(10);

H_up_hop(10,2) = V_No_ppp; 

%% Hamiltonian fürs Metall

xi = 0;

H_M = zeros(2);

H_M(1,1) = E_M_3d + splitting(4);     % xz

H_M(2,2) = E_M_3d + splitting(5);     % yz

H_M(1,2) = -1i*xi;

H_M(2,1) = +1i*xi;

%% hopping Metall - unitcell

H_M_u1 = zeros(2,10);
H_M_u2 = zeros(2,10);
H_M_u3 = zeros(2,10);
H_M_u4 = zeros(2,10);

%% xz orbital

H_M_u1(1,1) = -V_NM_pdp;
H_M_u3(1,1) =  V_NM_pdp;

%% yz orbital

H_M_u2(2,1) = -V_NM_pdp;
H_M_u4(2,1) =  V_NM_pdp;

%% Pi-Hamiltonian zusammenbauen

H_pi = [ H_up      , H_up_hop  , zeros(10) , H_up_hop' , H_M_u1'
         H_up_hop' , H_up      , H_up_hop  , zeros(10) , H_M_u2'
         zeros(10) , H_up_hop' , H_up      , H_up_hop  , H_M_u3'
         H_up_hop  , zeros(10) , H_up_hop' , H_up      , H_M_u4';
         H_M_u1    , H_M_u2    , H_M_u3    ,H_M_u4     , H_M  ];
     
%%


































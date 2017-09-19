function out = orbital(Metal,MO,X,Y,Z,repr)
addpath('../functions','../parameters')
%
% out = orbital(Metal,MO,X,Y,Z,repr)
%
% repr = 'imag', 'real'
%
% calculates the the molecular orbital MO of Metal 'Metal',
% at point X Y Z
% uses slater orbitals


%---------------- TO BE DELETED LATER ----------------------
% clear all
% Metal = 'Cu';
% finesse = 50;
% MO = 60;
% n_occ = 2;
%-----------------------------------------------------------

% finesse = length(X);

q_M3 = 8.3;  % 10
q_M4 = 2.7;
q_L_p = 3;
q_L_s = 3;
q_L_p_N = 3.6;
q_L_s_N = 3.6;

Q_pi = [q_L_p_N q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p_N ];
Q_pi = [Q_pi Q_pi Q_pi Q_pi];

Q_si_s	= [ q_L_s_N q_L_s q_L_s q_L_s q_L_s q_L_s q_L_s q_L_s q_L_s q_L_s_N ];
Q_si_s	= [Q_si_s Q_si_s Q_si_s Q_si_s];

Q_si_p	= [ q_L_p_N q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p q_L_p_N ];
Q_si_p	= [Q_si_p Q_si_p Q_si_p Q_si_s];

[H_sigma,H_pi] = Hamiltonian(Metal);

[psi_pi,E_pi] = eig(H_pi);

E_pi = [diag(E_pi),ones(length(E_pi),1),(1:length(E_pi))'];


[psi_sigma,E_sigma] = eig(H_sigma);
E_sigma = [diag(E_sigma),zeros(length(E_sigma),1),(1:length(E_sigma))'];

E = [E_sigma;E_pi];

E = [sortrows(E,1),(1:length(E))'];


if E(MO,2) == 0             % pi or sigma orbital? => sigma
    
	display('sigma')
	
    orb = E(MO,3);
   
    
    Psi = psi_sigma(:,orb);
    
    ind_1s = [13 17 21 25];
    ind_1s = [ind_1s 34+ind_1s 68+ind_1s 102+ind_1s];
    
    ind_2s = [1 4 7 10 14 18 22 26 29 32];
    ind_2s = [ind_2s 34+ind_2s 68+ind_2s 102+ind_2s];
    
    ind_2x = [2 5 8 11 15 19 23 27 30 33];
    ind_2x = [ind_2x 34+ind_2x 68+ind_2x 102+ind_2x];
    
    ind_2y = [3 6 9 12 16 20 24 28 31 34];
    ind_2y = [ind_2y 34+ind_2y 68+ind_2y 102+ind_2y];
    
    Rot_0=eye(34);
    a=fliplr(diag([-1,1]));
% 	a=fliplr(diag([1,-1]));
% 	a=fliplr(diag([1,1]));
% 	a = eye(2);
    Rot_m90=blkdiag(1,a,1,a,1,a,1,a,1,1,a,1,1,a,1,1,a,1,1,a,1,a,1,a);
    Rot_180=Rot_m90^2;
    Rot_90=Rot_m90^3;

    Rot=blkdiag(Rot_0,Rot_m90,Rot_180,Rot_90);
   
    Psi_4s    = Psi(end-3);
    Psi_3z2   = Psi(end-2);
    Psi_3x2y2 = Psi(end-1);
    Psi_3xy   = Psi(end);
    
    Psi = Psi(1:end-4);
    
    Psi = Rot*Psi;
   
    V    = zeros(size(X));

    C_L = coordinates('lig');
    C_H = coordinates('H');
    
    x_a = [C_L(:,1);0;0];
    y_a = [C_L(:,2);0;0];
    
    x_h = [C_H(:,1);0;0];
    y_h = [C_H(:,2);0;0];
    
%     if orb==75  % just for testing
%        
%         Psi = zeros(size(Psi));
%         Psi_4s = 0;
%         Psi_3z2 = 0;
%         Psi_3xy = 0;
%         
%     end
    
    for i=1:length(ind_2s)
        x = X - x_a(i);
        y = Y - y_a(i);
        z = Z;
        
        V = V + Psi( ind_2s(i) ).*slater_orb(x,y,z,Q_si_s(i),'2s') ...
              + Psi( ind_2x(i) ).*slater_orb(x,y,z,Q_si_p(i),'2px') ...
              + Psi( ind_2y(i) ).*slater_orb(x,y,z,Q_si_p(i),'2py');
    end
    
    V = V + Psi_4s.*slater_orb(X,Y,Z,q_M4,'4s') ...
          + Psi_3z2.*slater_orb(X,Y,Z,q_M3,'3z2') ...
          + Psi_3x2y2.*slater_orb(X,Y,Z,q_M3,'3x2y2') ...
          + Psi_3xy.*slater_orb(X,Y,Z,q_M3,'3xy');
    
    for i=1:length(ind_1s)
        x = X - x_h(i);
        y = Y - y_h(i);
        z = Z;
		V = V + Psi( ind_1s(i) ).*slater_orb(x,y,z,1,'1s');
    end
    
    out = V;
end

if E(MO,2) == 1             % pi or sigma orbital? => pi
    display('pi')
    orb = E(MO,3);
    
    if strcmp(repr,'imag')
        if orb==24
            Psi = -1/sqrt(2)*( psi_pi(:,24) + 1i * psi_pi(:,25) );
        elseif orb==25
            Psi = +1/sqrt(2)*( psi_pi(:,24) - 1i * psi_pi(:,25) );
        else
            Psi = psi_pi(:,orb);
        end
    elseif strcmp(repr,'real')
        Psi = psi_pi(:,orb);
    else
        return
    end
 
    
    C = coordinates('lig');
    
    x_a = [C(:,1);0;0];
    y_a = [C(:,2);0;0];
    
    V = zeros(size(X));
    
    % calculate volume matrix of wavefunction
    
    for i=1:length(Psi)-2
        x = X - x_a(i);
        y = Y - y_a(i);
        z = Z;
		V = V + Psi(i).*slater_orb(x,y,z,Q_pi(i),'2pz');
    end
    
    V = V + Psi(41).*slater_orb(X,Y,Z,q_M3,'3xz') ...
          + Psi(42).*slater_orb(X,Y,Z,q_M3,'3yz');
    
    out = V;
   
    
end















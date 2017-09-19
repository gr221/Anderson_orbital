function drho = anderson_orbital_pulse_L(t, rho)

global e_l U Gamma_sub Gamma_tip Gamma_tip2 B_x B_y B_z T mu0 A w D t0 t_del alpha
drho = zeros(6,1);    % a column vector P0 P1u P1d P2 C1ud C1du

%mu(1) is chemical potential in the tip and mu(2) in the substrate
mu(1)= mu0 + alpha*pulse(t);
mu(2) = mu0 - (1-alpha)*pulse(t);

Pauli_x = [0 1; 1 0];
Pauli_y = [0 -1i;1i 0];
Pauli_z = [1 0; 0 -1];

H1 = e_l*eye(2) - B_z*Pauli_z/2;
[V,H1d] = eig(H1); 
% Useful f+ and f- components of the Liouvillean

for i_LR = 1:2
    %This block is E_1l - E_0
    %(1) refers to the tip and (2) to the substrate
    fp{i_LR} = V*diag(fermi(diag(H1d)-mu(i_LR),T(i_LR)))*V';
    fm{i_LR} = eye(2)-fp{i_LR};
    p{i_LR} = V*diag(pp(diag(H1d)-mu(i_LR),T(i_LR)))*V';
    
    %This is E_2 - E_1l
    fUp{i_LR} = V*diag(fermi(2*e_l+U-mu(i_LR)-diag(H1d),T(i_LR)))*V';
    fUm{i_LR} = eye(2)-fUp{i_LR};
    pU{i_LR} = V*diag(pp(2*e_l+U-mu(i_LR)-diag(H1d),T(i_LR)))*V';
   
end

rho_0 = rho(1);
rho_1 = [rho(2),rho(5);rho(6),rho(3)]; 
rho_2 = rho(4);

drho_0 = 0;
drho_1 = -1i*(H1*rho_1 - rho_1*H1);
drho_2 = 0;

drho_0 = -trace(Gamma_sub*fp{2}+Gamma_tip*fp{1})*rho_0 + trace(Gamma_sub*fm{2}*rho_1) +...
            0.5*trace(fm{1}*(Gamma_tip*rho_1+rho_1*Gamma_tip))+ ...
            1i/(2*pi)*trace(p{1}*(Gamma_tip*rho_1 - rho_1*Gamma_tip));
for i_LR = 1:2
    if i_LR == 1
        Gamma = Gamma_tip;
        Gamma2 = Gamma_tip2;
    else
        Gamma = Gamma_sub;
        Gamma2 = Gamma_sub;
    end
    
    drho_1 = drho_1 - Gamma/2*(fm{i_LR} - 1i/pi *p{i_LR})*rho_1 - ...
                0.5*rho_1'*(fm{i_LR}' + 1i/pi*p{i_LR}')*Gamma' -...
                0.5*Gamma2*(fUp{i_LR} - 1i/pi*pU{i_LR})*rho_1 -...
                0.5*rho_1'*(fUp{i_LR}' + 1i/pi*pU{i_LR})*Gamma2' +...
                0.5*(fp{i_LR}*Gamma+Gamma*fp{i_LR} + 1i/pi*(p{i_LR}*Gamma - Gamma*p{i_LR}))*rho_0 +...
                0.5*(fUm{i_LR}*Gamma2+Gamma2*fUm{i_LR} + 1i/pi*(pU{i_LR}*Gamma2 - Gamma2*pU{i_LR}))*rho_2;
                
end       

drho_2 = -trace(Gamma_sub*fUm{2}+Gamma_tip2*fUp{1})*rho_2 + ...
        trace(Gamma_sub*fUp{2}*rho_1) + ...
        0.5*trace(fUp{1}*(Gamma_tip2*rho_1 + rho_1*Gamma_tip2)) +...
        i/(2*pi) * trace(pU{1}*(Gamma_tip2*rho_1 - rho_1*Gamma_tip2));
    
drho(1) = drho_0;
drho(2:3) = diag(drho_1);
drho(4) = drho_2;
drho(5) = drho_1(1,2);
drho(6) = drho_1(2,1);
    
end

clear all
close all

global e_l U Gamma_sub Gamma_tip Gamma_tip2 B_x B_y B_z theta  T mu0 A w D t0 alpha t_del
% Bath parameters
Gamma_0 = 0.001; %bare bath-impurity tunnelling rates (meV)
T=[0.1,0.1]; %Temperature of the tip and substrate
e_l = 0.2 %onsite energy
U = 1 %Interaction 
B_x = 0; % magnetic field components
B_y = 0; % quantization axis of the dot orthogonal to the one of the leads
B_z = 0.001;
B = sqrt(B_x^2 + B_y^2 + B_z^2); % strength of the magnetic field
mu0 = e_l+U/2; %Chemical potential

%The gamma Matrix for the substrate
Gamma_sub = Gamma_0*eye(2);

%We use the Homo, Lumo1 and Lumo2 for Copper by the orbital
%function written by Benni
% for x=-12:12
%     for y = -12:12
%         for z = -8:8
% r_tip = [x,y,z];
r_tip = [-11,-12,-8];
HOMO = orbital('Cu',98,r_tip(1),r_tip(2),r_tip(3),'imag')
LUM1 = orbital('Cu',99,r_tip(1),r_tip(2),r_tip(3),'imag')
LUM2 = orbital('Cu',100,r_tip(1),r_tip(2),r_tip(3),'imag')

Gamma_tip = Gamma_0*[conj(LUM1)*LUM1, conj(LUM1)*LUM2; ...
    conj(LUM2)*LUM1, conj(LUM2)*LUM2];
Gamma_tip2 = Gamma_0*[conj(LUM1)*LUM1, -conj(LUM1)*LUM2; ...
    -conj(LUM2)*LUM1, conj(LUM2)*LUM2];

%The Pauli matrices
Pauli_x = [0 1; 1 0];
Pauli_y = [0 -1i;1i 0];
Pauli_z = [1 0; 0 -1];

A = 0.7; % amplitude of the pulse

w = 0.1*Gamma_0; % frequency of the pulse
D = 1/Gamma_0; % width of the pulse
t0 = 5/Gamma_0; % center of the pulse

N = [0 1 1 2]; %possible occupation numbers on the metal
t_max = t0 + 30*D;

%Calculation of the analytical stationary states
        
H_S = blkdiag(0, e_l*eye(2) - B_z*Pauli_z/2, 2*e_l + U)

alpha = 0.85; %fraction of bias drop at the tip

Current = [];
t_vec = [];

%statistical density matrix
rho_Stat = expm(-1/T(1)*(H_S - mu0*diag(N)));       
rho_Stat = rho_Stat/trace(rho_Stat);
% the initial condition is written as rho_0 rho_1uu rho_1dd rho_2 rho_1ud rho_1du
rho_In = [diag(rho_Stat)',rho_Stat(2,3),rho_Stat(3,2)] 

options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[t_vec,rho] = ode45(@anderson_orbital_pulse_L,[0 t_max],rho_In,options);

t_del_vec = 0;
i_td = 1;

for i_t= 1:length(t_vec)
    t = t_vec(i_t);
    rho_0 = rho(i_t,1);
    rho_1 = [rho(i_t,2),rho(i_t,5);rho(i_t,6),rho(i_t,3)];
    rho_2 = rho(i_t,4);
    
    
    mu(1)= mu0 + alpha*pulse(t);
    mu(2) = mu0 - (1-alpha)*pulse(t);
    
    
    
    H1 = e_l*eye(2) - B_z*Pauli_z/2;
    [V,H1d] = eig(H1);

    % Useful f+ and f- components of the Liouvillean
    
    for i_LR = 1:2
        %This block is E_1l - E_0
        fp{i_LR} = V*diag(fermi(diag(H1d)-mu(i_LR),T(i_LR)))*V';
        fm{i_LR} = eye(2)-fp{i_LR};
        p{i_LR} = V*diag(pp(diag(H1d)-mu(i_LR),T(i_LR)))*V';

        %This is E_2 - E_1l
        fUp{i_LR} = V*diag(fermi(2*e_l+U-mu(i_LR)-diag(H1d),T(i_LR)))*V';
        fUm{i_LR} = eye(2)-fUp{i_LR};
        pU{i_LR} = V*diag(pp(2*e_l+U-mu(i_LR)-diag(H1d),T(i_LR)))*V';
    end   
            
    %first line fUP oder fUm{1}
    drho_2 = -trace(Gamma_sub*fUm{2}+Gamma_tip2*fUp{1})*rho_2 + ...
    trace(Gamma_sub*fUp{2}*rho_1) + ...
    0.5*trace(fUp{1}*(Gamma_tip2*rho_1 + rho_1*Gamma_tip)) +...
    1i/(2*pi) * trace(pU{1}*(Gamma_tip2*rho_1 - rho_1*Gamma_tip2));

    % drho_1 = -1i*(H1*rho_1-rho_1*H1);
        
    for i_LR = 1:2
        if i_LR == 1
            Gamma = Gamma_tip;
            Gamma2 = Gamma_tip2;

            % drho_2 = -trace(Gamma_tip*fUp{1})*rho_2 + ...
            % 0.5*trace(fUp{1}*(Gamma_tip*rho_1 + rho_1*Gamma_tip)) +...
            % 1i/(2*pi) * trace(pU{1}*(Gamma_tip*rho_1 - rho_1*Gamma_tip));
        else
            Gamma = Gamma_sub;
            Gamma2 = Gamma_sub;
            % drho_2 = -trace(Gamma_sub*fUm{2})*rho_2 + ...
            % trace(Gamma_sub*fUp{2}*rho_1);
        end

        drho_1 =  Gamma/2*(fm{i_LR} - 1i/pi *p{i_LR})*rho_1 - ...
            0.5*rho_1'*(fm{i_LR}' + 1i/pi*p{i_LR}')*Gamma' -...
            0.5*Gamma2*(fUp{i_LR} - 1i/pi*pU{i_LR})*rho_1 -...
            0.5*rho_1'*(fUp{i_LR}' + 1i/pi*pU{i_LR})*Gamma2' +...
            0.5*(fp{i_LR}*Gamma+Gamma*fp{i_LR} + 1i/pi*(p{i_LR}*Gamma - Gamma*p{i_LR}))*rho_0 +...
            0.5*(fUm{i_LR}*Gamma2+Gamma2*fUm{i_LR} + 1i/pi*(pU{i_LR}*Gamma2 - Gamma2*pU{i_LR}))*rho_2;
        % drho_1 = drho_1 - Gamma/2*(fm{i_LR} - i/pi *p{i_LR})*rho_1 - ...
        % 0.5*rho_1'*(fm{i_LR}' + i/pi*p{i_LR}')*Gamma' -...
        % 0.5*Gamma*(fUp{i_LR} - i/pi*p{i_LR})*rho_1 -...
        % 0.5*rho_1'*(fUp{i_LR}' + i/pi*p{i_LR})*Gamma' +...
        % 0.5*(fp{i_LR}*Gamma+Gamma*fp{i_LR} + i/pi*(p{i_LR}*Gamma - Gamma*p{i_LR}))*rho_0 +...
        % 0.5*(fUm{i_LR}*Gamma+Gamma*fUm{i_LR} + i/pi*(pU{i_LR}*Gamma - Gamma*pU{i_LR}))*rho_2;


        Current(i_t,i_LR) = trace(eye(2)*drho_1) + 2*drho_2;
        % Current = trace(eye(2)*drho_1) + 2*drho_2
    end
end

Q(i_td,:) = squeeze(trapz(t_vec,Current)) % # of electrons that enetered (>0) of left (<0)
% the impurity. Q_ij with i = L,R and j = u,d



if length(t_del_vec) == 1
    Nf=20; 
    % figure(1)
    % subplot(2,2,1)
    % plot(t_vec*Gamma_0,rho(:,1:4),'Linewidth',2);
    % hold on
    %
    % plot(t_vec*Gamma_0,repmat(diag(rho_Stat)',length(t_vec),1),'k-')
    %
    % legend('P_0','P_{1\uparrow}','P_{1\downarrow}','P_2')
    
    %%%%% Graphical makeup %%%%%%%%%%%%%%%%%%%%%%%%%%
    % ylim([0 1])
    % xlabel('time [1/\gamma]','Fontsize',Nf)
    % ylabel('Populations','Fontsize',Nf)
    % set(gca,'Fontsize',Nf)
    
    %    subplot(2,2,2)
    %     plot(t_vec*Gamma_0,real(rho(:,5:6)))
    %     hold on
    %     plot(t_vec*Gamma_0,imag(rho(:,5:6)),'--')
    %
    %     %%%%% Graphical makeup %%%%%%%%%%%%%%%%%%%%%%%%%%
    %     xlabel('time [1/\gamma]','Fontsize',Nf)
    %     ylabel('Coherences','Fontsize',Nf)
    %     set(gca,'Fontsize',Nf)
    
    % subplot(2,2,3)
    %
    % plot(t_vec*Gamma_0,[alpha*pulse(t_vec)';(alpha-1)*pulse(t_vec)'],'Linewidth',2)
    %
    % legend('\Delta\mu_L','\Delta\mu_R')
    
    %%%%% Graphical makeup %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %title('Pulse','Fontsize',Nf)
    % xlabel('time [1/\gamma]','Fontsize',Nf)
    % ylabel('\Delta\mu_\alpha','Fontsize',Nf)
    % set(gca,'Fontsize',Nf)
    %
    % hold on
    % rho_Stat_tilde = [rho_Stat(1,1), rho_Stat(2,2)+rho_Stat(3,3),rho_Stat(4,4)];
    % if max(rho_Stat_tilde) - rho_Stat_tilde(1) == 0
    %     plot(t_vec*Gamma_0,(e_l-B/2-mu0)*ones(1,length(t_vec)),'g--')
    %     plot(t_vec*Gamma_0,(e_l+B/2-mu0)*ones(1,length(t_vec)),'g--')
    % end
    % if max(rho_Stat_tilde) - rho_Stat_tilde(2) == 0
    %     plot(t_vec*Gamma_0,(e_l + B/2 +U-mu0)*ones(1,length(t_vec)),'g--')
    %     plot(t_vec*Gamma_0,(e_l - B/2 +U-mu0)*ones(1,length(t_vec)),'g--')
    %     plot(t_vec*Gamma_0,(e_l-B/2-mu0)*ones(1,length(t_vec)),'r--')
    %     plot(t_vec*Gamma_0,(e_l+B/2-mu0)*ones(1,length(t_vec)),'r--')
    % end
    % if max(rho_Stat_tilde) - rho_Stat_tilde(3) == 0
    %     plot(t_vec*Gamma_0,(e_l+B/2 + U-mu0)*ones(1,length(t_vec)),'r--')
    %     plot(t_vec*Gamma_0,(e_l-B/2 + U-mu0)*ones(1,length(t_vec)),'r--')
    % end
    % hold off
    %
    %
    % subplot(2,2,4)
    
    % for i_t = 1:length(t_vec)
    %     rho_t = [rho(i_t,2),rho(i_t,5);rho(i_t,6),rho(i_t,3)];
    %     P1(i_t,:) = sort(eig(rho_t),'descend');
    %     S_x(i_t) = 0.5*trace(Pauli_x*rho_t);
    %     S_y(i_t) = 0.5*trace(Pauli_y*rho_t);
    %     S_z(i_t) = 0.5*trace(Pauli_z*rho_t);
    % end
    %
    % plot(t_vec*Gamma_0,[rho(:,1),P1,rho(:,4)],'Linewidth',2)
    %
    % %%%%% Graphical makeup %%%%%%%%%%%%%%%%%%%%%%%%%%
    % legend('P_0','P_{1a}','P_{1b}','P_2')
    % xlabel('time [1/\gamma]','Fontsize',Nf)
    % ylim([0 1])
    % ylabel('Eigenvalues of \rho','Fontsize',Nf)
    % set(gca,'Fontsize',Nf)
    %
    %
    %
    % subplot(2,2,2)
    % plot(Gamma_0*t_vec,[S_x;S_y;S_z],'Linewidth',2)
    % %%%%% Graphical makeup %%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % legend('\langleS_x\rangle ','\langleS_y\rangle ','\langleS_z\rangle ')
    % xlabel('time [1/\gamma]','Fontsize',Nf)
    % ylabel('\langleS_\alpha\rangle','Fontsize',Nf)
    % set(gca,'Fontsize',Nf)
    
    
    
    
    
    %%
    % G_rel = sum(Gamma(:))*exp(-U/(2*T(1)));
    % Sz_eq = 0.5*tanh(B/(2*T(1)));
    
    % G_rel = sum(Gamma(:))/2*0.5*(1 - fermi(e_l + B/2-mu0,T(1)) + ...
    %     1 - fermi(e_l - B/2-mu0,T(1)) + ...
    %     fermi(e_l +U - B/2-mu0,T(1)) +...
    %     fermi(e_l +U + B/2-mu0,T(1)) );
    % Sz_eq = -0.5*(1 - fermi(e_l - B/2-mu0,T(1)) -  ...
    %     1 + fermi(e_l + B/2-mu0,T(1)) + ...
    %     fermi(e_l +U + B/2-mu0,T(1)) -...
    %     fermi(e_l +U - B/2-mu0,T(1)))...
    %     /...
    %     (1 - fermi(e_l + B/2-mu0,T(1)) + ...
    %     1 - fermi(e_l - B/2-mu0,T(1)) + ...
    %     fermi(e_l +U - B/2-mu0,T(1)) +...
    %     fermi(e_l +U + B/2-mu0,T(1)) );
    % domega = sum(Gamma(:))/2*(pp(e_l - B/2-mu0,T(1)) - pp(e_l + B/2-mu0,T(1))...
    %     + pp(e_l +U + B/2-mu0,T(1)) - pp(e_l +U - B/2-mu0,T(1)))/(2*pi);
    %
    % omega = B + domega;
    
    % Current
    figure(2)  % Current(i_t,i_LR,i_ud)
    K = 1.6/4.1*10^5; %(e/hbar in \nC/s*eV)
    plot(Gamma_0*t_vec,K*Current(:,1))
    hold on
    plot(Gamma_0*t_vec,K*Current(:,2))
    legend('I_{L\uparrow}','I_{R\uparrow}','I_{L\downarrow}','I_{R\downarrow}')
    xlabel('time [1/\gamma]','Fontsize',Nf)
    ylabel('Current [nA]','Fontsize',Nf)
    waitfor(figure(2))
    hold off
    
else
    %%
    figure(2)
    subplot(1,2,1)
    plot(Gamma_0*t_del_vec,Q(:,:,1))
    hold on
    plot(Gamma_0*t_del_vec,Q(:,:,2))
    xlabel('delay time [1/\gamma]','Fontsize',Nf)
    ylabel('Transf. charge/cycle [e]','Fontsize',Nf)
    legend('Q_{L\uparrow}','Q_{R\uparrow}','Q_{L\downarrow}','Q_{R\downarrow}')
    hold off
    
    subplot(1,2,2)
    plot(Gamma_0*t_del_vec,Q(:,1,1) + Q(:,1,2))
    hold on
    plot(Gamma_0*t_del_vec,Q(:,2,1) + Q(:,2,2))
    xlabel('delay time [1/\gamma]','Fontsize',Nf)
    ylabel('Transf. charge/cycle [e]','Fontsize',Nf)
    legend('Q_{L}','Q_{R}')
    hold off
    
    
    
end
% end
% end
% end

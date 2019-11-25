%%%% Transport amplitude %%%
clearvars;
close all;

n_k = 200;  % step of wave vector
n_p = 200;

%%% Hopping
V = 10;  % hopping lead
gamma = 0.25;    % hopping lead and conductor
delta = 0.0;    % hopping diammer
E_a = 0.1;     % on site of phi_a
E_b = -0.1;     % on site of phi_b

%%% Magnetic field
psi_list = linspace(-pi,pi,n_p);

%%% Wave vector space
E_list = linspace(-0.5,0.5,n_k);
k_list = linspace(-0.1,0.1,n_k);

res_vec = zeros(4,n_k,n_p);
for i = 1:n_k
    for j =1:n_p
        E = E_list(i);
        a = -E/V + 1i*sqrt(1-E^2/V^2);
        ad = -E/V - 1i*sqrt(1-E^2/V^2);
        gamma = gamma*exp(1i*psi_list(j));
    
        vec = [V/2; 0; conj(gamma)*ad; gamma*ad];
    
        transport_matrix = [[-V/2,0,-gamma,-conj(gamma)];...
                        [0,-V/2,-conj(gamma),-gamma];...
                        [-conj(gamma)*a,-gamma*a,E_a+(V/2)*(a+ad),-delta];...
                        [-gamma*a,-conj(gamma)*a,-delta,E_b+(V/2)*(a+ad)]];
    
    
        res_vec(:,i,j)=transport_matrix\vec;
    end
end

fig1 = figure();
pcolor(psi_list,E_list,squeeze(abs(res_vec(1,:,:)).^2));
colorbar
shading flat
ylabel('k[\pi]');
xlabel('e^{i\psi} [\psi:\pi]');
title('|r|^2');

fig2 = figure();
pcolor(psi_list/pi,E_list,squeeze(abs(res_vec(2,:,:)).^2));
colorbar
shading flat
ylabel('k[\pi]');
xlabel('e^{i\psi} [\psi:\pi]');
title('|t|^2');

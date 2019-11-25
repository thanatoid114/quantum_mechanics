%%%% Transport amplitude %%%
clearvars;
close all;

n_k = 400;  % step of wave vector

%%% Hopping
V = 10;  % hopping lead
gamma = 0.25;    % hopping lead and conductor
delta = 0.1;    % hopping diammer
E_a = 0.;     % on site of phi_a
E_b = 0.;     % on site of phi_b

%%% Magnetic field
psi = 0.0;    % Effective phase
gamma = gamma*exp(1i*psi);

%%% Wave vector space
E_list = linspace(-0.5,0.5,n_k);

res_vec = zeros(4,n_k);
for i =1:n_k
%     k = k_list(i);
    E = E_list(i);
    a = -E/V + 1i*sqrt(1-E^2/V^2);
    ad = -E/V - 1i*sqrt(1-E^2/V^2);
    
    vec = [V/2; 0; conj(gamma)*ad; gamma*ad];
    
    transport_matrix = [[-V/2,0,-gamma,-conj(gamma)];...
                        [0,-V/2,-conj(gamma),-gamma];...
                        [-conj(gamma)*a,-gamma*a,E_a+(V/2)*(a+ad),-delta];...
                        [-gamma*a,-conj(gamma)*a,-delta,E_b+(V/2)*(a+ad)]];
    
    
    res_vec(:,i)=transport_matrix\vec;
end

fig1 = figure();
subplot(2,1,1)
hold on
% plot(k_list/(pi),abs(res_vec(1,:)).^2./(abs(res_vec(1,:)).^2+abs(res_vec(2,:)).^2));
% plot(k_list/(pi),abs(res_vec(2,:)).^2./(abs(res_vec(1,:)).^2+abs(res_vec(2,:)).^2));
plot(E_list/V,abs(res_vec(1,:)).^2);
plot(E_list/V,abs(res_vec(2,:)).^2);
legend('|r|^2','|t|^2');
xlabel('k[\pi]')
hold off

subplot(2,1,2)
hold on
% plot(k_list/(pi),abs(res_vec(3,:)).^2./(abs(res_vec(3,:)).^2+abs(res_vec(4,:)).^2));
% plot(k_list/(pi),abs(res_vec(4,:)).^2./(abs(res_vec(3,:)).^2+abs(res_vec(4,:)).^2));
plot(E_list/(V),abs(res_vec(3,:)).^2);
plot(E_list/(V),abs(res_vec(4,:)).^2);
legend('|\phi_a|^2','|\phi_b|^2');
xlabel('k[\pi]')
hold off

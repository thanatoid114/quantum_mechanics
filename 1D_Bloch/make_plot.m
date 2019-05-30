%%%%% Make the plot %%%%%
clearvars
close all

load parameters
load result

fig1 = figure();
for i = 1:N_bands
    hold on
plot(kx_list/kxmax,real(eig_v(:,i)));
% plot(kx_list/kxmax,imag(eig_v(:,i)),'red');
xlim([-1 1])
end

hold off

fig2 = figure();
for i = 1:N_bands
    hold on
% plot(kx_list/kxmax,real(eig_v(:,i)),'blue');
plot(kx_list/kxmax,imag(eig_v(:,i)));
xlim([-1 1])
end

hold off


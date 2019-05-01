%%%%% Make the plot %%%%%
clearvars
close all

load parameters
load result

fig1 = figure();
for i = 1:5
hold on
plot(x_coord,eig_v(i)+100*(abs(eig_mode(:,i).^2)));
hold off
end

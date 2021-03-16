clearvars
close all

n = 3; % number of lattice
t = 1; % hopping of lattice in one unit cell
dt = 1; % hopping between the cell
dl = 1.2; % hopping of between lead and lattice left side
dr = 1.1; % hopping of between lead and lattice right side

E=1.3;
eta = 0.1*1i;


% Hamilton of cross lattice for one unit cell
H_D = [[0,t];[t,0]];
% Hamiltion with the lattice
Sigma_L = 
% isolate green function for the cross lattice
% g_isol = (E+eta)*one(2) - 


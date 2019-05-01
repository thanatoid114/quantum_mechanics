%%%%%% Continus mode of 1D quantum well
% particle in a infinite deep quantum well; using Dirichlet boundary condition

clearvars;
close all;

% Constant
hbar = 1;   % planck's constant
a = 1;  % lattice periodic
m = 1; % effective mass of particle
ku = pi/a; % half the reciprocal lattice vector
Eu = hbar*ku^2/(2*m);   % energy unit

% Mesh of system
nx = 128;
dx = a/nx;
x_coord = 0:dx:a;  % perodic direction [0,1]
Nx = length(x_coord);

save parameters;

[partial_x,partial_xx] = partial(Nx,dx); % Partial operator of [X,Y]
bloch_matrix = full(-(1/2)*(partial_xx));
[mode,vals] = eig(bloch_matrix);
vals = diag(vals);
[~,order] = sort(real(vals));
vals = vals(order);
eig_v = vals/Eu;
eig_mode = mode(:,order);

save('result.mat','eig_v','eig_mode','-v7.3');

make_plot;
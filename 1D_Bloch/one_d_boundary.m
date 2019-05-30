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

% lattice system
v1 = 5*1i*Eu; % height of barrier
a1 = a/4; % width of barrier
a2 = 3*a/4; % width of barrier
lattice = potential(x_coord,v1,a1,a2);
lattice = lattice - (v1/2)*ones(size(lattice));

% K space
Nk = 64;
kxmax = pi/a;
dk = 2*kxmax/Nk;
kx_list = -kxmax:dk:kxmax-dk;

N_bands = 2;

save parameters;

eig_v = zeros(Nk,N_bands);
eig_mode = zeros(Nk,Nx,N_bands);
lattice = reshape(lattice,Nx,1);
[partial_x,partial_xx] = partial(Nx,dx); % Partial operator of [X,Y]

for i = 1:Nk
    kx = kx_list(i);
    bloch_matrix = spdiags(lattice,0,Nx,Nx)-(1/2)*((-kx^2)*speye(Nx,Nx)...
        + partial_xx + 2*1i*(kx*partial_x));
    [mode,vals] = eigs(bloch_matrix,N_bands,'sm');
    vals = diag(vals);
    [~,order] = sort(real(vals));
    vals = vals(order);
    eig_v(i,:) = vals/Eu;
    eig_mode = mode(:,order);
end
save('result.mat','eig_v','eig_mode','-v7.3');

make_plot;
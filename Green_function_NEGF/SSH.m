%%%% Green function approach for SSH Chain
clearvars
close all

dE = 0.4;   % Energy round 
En = 801;   % Energy list for calculation
iter_lead = 100;    % Iteration times for leads to calculate the Surface Green's function
eta = 0.000001;  % Imaginary value for Green function
%%%% for contract cell
w = 1.0;    % inter cell hopping
v = 1.0;    % intral cell hopping
V0 = 0.25;    % on site energy a site (V0,-V0)
NC = 100;     % Number of contact cell

%%%% for lead
wld = 1.0;   % inter cell hopping for lead
vld = 1.0;   % intra cell hopping for lead
V0ld = 0.25;   % on site energy for lead (V0ld,-V0ld)


%%%% Define the element matrix to construct the Hamiltonian
s0 = eye(2);    % diagonal matrix
sx = [0 1; 1 0];    % Sigma matrix 1
sy = [0 -1i; 1i 0]; % Sigma matrix 2
sz = [1 0; 0 -1];   % Sigma matrix 3

%%%% Hamiltonian for single cell and the hopping to the next unit cell
Hld_intra = sz*V0ld - vld*sx;     % on site and intra cell hopping
Hld_inter = (-1/2)*(sx-1i*sy)*wld;  % for the inter cell hopping
Hsys_intra = sz*V0 - v*sx;
Hsys_inter = (-1/2)*(sx-1i*sy)*w;
%%%%% Hamiltonian for contact in k space %%%%%%%%%%
Nk = 101;   % number of points in k space for plot
klist = linspace(-pi,pi,Nk);    % k list
eig_val = zeros(2,Nk);
for i = 1:Nk
    k = klist(i);
    Hk = [V0 -w-v*exp(1i*k); -w-v*exp(-1i*k) -V0];
    [~,vals] = eig(Hk);
    eig_val(:,i) = diag(vals);
end
% hold on
% plot(klist/pi,eig_val(1,:),'-o');
% plot(klist/pi,eig_val(2,:),'-o');
% hold off

% Find the min and max of the band; where is at k=0;
Emin = min(eig_val(1,:))-dE;
Emax = max(eig_val(2,:))+dE;

%%%% Define the energy list for green function
Elist = linspace(Emin,Emax,En);
Glist = zeros(size(Elist));
%%%% For Surface Green's function
%%%% Based on Eq:4.60 Give the definition for the matrix block of Green
%%%% function;   The iteration is based on Eq:4.67
for i = 1:En
    en = (Elist(i)+1i*eta)*s0;
    % The Lift side G11
    % For initial term
    d1 = en - Hld_intra;  % terms for the surface of lead column (cell)
    D1 = en - Hld_intra;  % terms for hoppings for a single isolated column (cell)
    A1 = Hld_inter;  % terms for hoppings for neighbour column (cell)
    B1 = Hld_inter'; % terms for hoppings for neighbour column (cell) H.C. term
    % By iteration calculate the green function g11 Eq.4.67
    for n_iter = 1:iter_lead
        d2 = d1 - (A1/D1)*B1;
        D2 = D1 - (A1/D1)*B1 - (B1/D1)*A1;
        A2 = (A1/D1)*A1;
        B2 = (B1/D1)*B1;
        
        d1 = d2;
        A1 = A2;
        B1 = B2;
        D1 = D2;
    end
    g11_iso = s0/d1;    % Green's function surface term Eq.4.69 
    SigL = -d1 + en - Hld_intra;    % Self energy Eq.4.56
    GammaL =  1i*(SigL - SigL');    % Eq.4.20
    
    % The Right side (swap A1 and B1) GNN
    % For initial term
    d1 = en - Hld_intra;  % terms for the surface of lead column (cell)
    D1 = en - Hld_intra;  % terms for hoppings for a single isolated column (cell)
    A1 = Hld_inter';  % terms for hoppings for neighbour column (cell)
    B1 = Hld_inter; % terms for hoppings for neighbour column (cell) H.C. term
    for n_iter = 1:iter_lead
        d2 = d1 - (A1/D1)*B1;
        D2 = D1 - (A1/D1)*B1 - (B1/D1)*A1;
        A2 = (A1/D1)*A1;
        B2 = (B1/D1)*B1;
        
        d1 = d2;
        A1 = A2;
        B1 = B2;
        D1 = D2;
    end
    gNN_iso = s0/d1;    % Green's function surface term Eq.4.69 
    SigR = -d1 + en - Hld_intra;    % Self energy Eq.4.56
    GammaR =  1i*(SigR - SigR');    % Eq.4.20
    
    % Computing the all green function
    G11 = g11_iso;  % The initial term
    G10 = g11_iso;  % The initial term
    for j = 1:NC
        gii_iso = s0/(en-Hsys_intra);   % green's function for isolate term
        G22 = (s0/(s0-gii_iso*Hsys_inter*G11*Hsys_inter'))*gii_iso;   % Eq in P57
        G21 = G22*Hsys_inter*G10;
        G11 = G22;      % For interation
        G10 = G21;      % For iteration
    end
    
    % For the last term right lead
    G22 = (s0/(s0 - gNN_iso*Hsys_inter*G11*Hsys_inter'))*gNN_iso;
    G21 = G22*Hsys_inter*G10;
    g_t = real(trace(GammaL*G21'*GammaR*G21));        % Transmission coefficient Eq.4.44
    Glist(1,i) = g_t;
end

figure();
subplot(1,2,1)
plot(Glist,Elist,'.-k');
title('Transmission coefficient');
ylabel('Energy');
subplot(1,2,2)
hold on
plot(klist/pi,eig_val(1,:),'.-b');
plot(klist/pi,eig_val(2,:),'.-b');
title('band structure for contract cell');
xlabel('wave vector');
hold off
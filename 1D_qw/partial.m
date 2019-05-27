function [ partial_x, partial_xx ] = partial( Nx,dx )
% Generate the partial function for 1D problem;

e1 = ones(Nx,1);
partial_x = spdiags([e1,-1*e1],[1,-1],Nx,Nx)/(2*dx);
% partial_x(1,:) = 0;     % boundary condition
% partial_x(end,:) = 0;   % boundary condition

partial_xx = spdiags([e1,e1,-2*e1],[1,-1,0],Nx,Nx)/(dx*dx);    % partial operator partial_xx;
% partial_xx(1,:) = 0;    % periodic boundary condition;
% partial_xx(end,:) = 0; % periodic boundary condition;
% partial_xx(:,1) = 0;
% partial_xx(:,end) = 0;
end


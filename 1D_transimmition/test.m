% test
clearvars
close all
time=0:0.002:100;
xmax=5;
Nx=256;
dx=2*xmax/Nx;
x=(-Nx/2:1:Nx/2-1)'*dx;
dk=pi/xmax;
kx=(-Nx/2:1:Nx/2-1)'*dk;
Ek=kx.^2/2000;
%   plot(kx,Ek); stop

ki = -60;  
xi = -2;

  
psi_ini=exp(-1i*ki*x).*exp(-(x-xi).^2);
psi_ini=psi_ini/max(abs(psi_ini));

v_fun = @(x1,x2,v) ((x>x1)&(x<x2)).*ones(size(x))*v;
vp = v_fun(1,1+2*dx,2.);

odefcn = @(t,psi) ...
    -1i*(fftshift(ifft2(ifftshift(Ek.*fftshift(fft2(ifftshift(psi))))))...
    +vp.*psi);



%  plot(x,abs(psi_ini)); stop


[t,psi]=ode45(@(t,psi) ...
    odefcn(t,psi)...
    ,time,psi_ini);


save parameters

%%%% Plot
h1=figure(1);
plot(t,sum(abs(psi),2),'b-o');
xlabel('time');
ylabel('Amplitude');






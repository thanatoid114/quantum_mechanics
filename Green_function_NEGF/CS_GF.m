loop=1;

if loop==1

eta=0.000001;
d0=1.0; d0c=1.0; t0=1.0; em=0.0; ep=0.0; 
V0=1.0; tl=t0+0.1; eml=em; gam=1.0; gamp=0.0;
Nx=100;
Vd=0.0;


enmin=-6;
enmax=5;

ielen=1000;
iNlist=[Nx];
iNlen=1;
iNda=1;

s0=[1 0;0 1];sx=[0 1;1 0];sy=[0 -1i;1i 0];sz=[1 0;0 -1];

% for leads
hl0=eml*sz-tl*sx;
hlx=-V0*(s0+sx);

Hlx=hl0;
Blx=hlx;

% for system
h0=em*sz-t0*sx;
hx=-d0*s0-d0c*sx;

Hx=h0;
Bx=hx;    

Yy=eye(2);iclen=1;

Yy=zeros(2,2); Yy(1:1,1:1)=eye(1); 
Yy(end:end,end:end)=eye(1); 
    
listu=zeros(1,ielen);
listumean=zeros(1,ielen);    
listuimag=zeros(1,ielen);    
listth=zeros(1,ielen-1);    

for ie=1:ielen
        ie;
        en=enmin+(enmax-enmin)*(ie-1)/(ielen-1)+1i*eta;% choose one energy for lead! (E-i\eta)

        d1=en*Yy-Hlx;A1=-Blx; B1=-Blx'; D1=en*Yy-Hlx; % Eq.4.60  !!!!!!!!! A and B  have minus sign???
        for i=1:1:100
            d2=d1-A1/D1*B1;             % Eq.4.67(1)
            D2=D1-A1/D1*B1-B1/D1*A1;    % Eq.4.67(2)
            A2=A1/D1*A1;                % (3)
            B2=B1/D1*B1;                % (4)

            d1=d2;A1=A2;B1=B2;D1=D2;
        end
        SigL=-d1+en*Yy-Hlx; %Self energy from the Lead??????????????????????????? Eq.4.56
        GammL=1i*(SigL-SigL'); %Gamma for left lead. will be used for later.
        g11=Yy/d1; %Green's function for the first strip.  

        d1=en*Yy-Hlx;A1=-Blx'; B1=-Blx; D1=en*Yy-Hlx;
        for i=1:1:100
            d2=d1-A1/D1*B1;
            D2=D1-A1/D1*B1-B1/D1*A1;
            A2=A1/D1*A1;
            B2=B1/D1*B1;

            d1=d2;A1=A2;B1=B2;D1=D2;
        end
        gNN=Yy/d1;%Green's function for the last strip.
        SigR=-d1+en*Yy-Hlx; %Self energy from the right lead.????????????????????
        GammR=1i*(SigR-SigR'); %Gamma of right lead, used later.

        G10=g11;G11=G10;


        hx=-(d0*s0+d0c*sx);Bx=hx;%%%%%%%%%%%%%%%%%%%%% hopping to next lead
        for i=2:1:Nx %Computing trasfer matrix along x. 
            Vyx=0;
            g22=Yy/(en*Yy-Hx-Vyx);      %%%????????????????? Sigma? Vyx????Eq.4.56
            G22=Yy/(Yy-g22*Bx*G11*Bx')*g22;
            G21=G22*Bx*G10;
            G11=G22;G10=G21;
        end

        G22=Yy/(Yy-gNN*Bx*G11*Bx')*gNN;	% Eq.4.57???????
        G21=G22*Bx*G10;
        g=real(trace(G21'*GammR*G21*GammL)); %this is final conductance. % Eq.4.44
        gimag=imag(trace(G21));
        listc=g;
        
        
        listu(1,ie)=real(en);
        listumean(1,ie)=listc(:,1);
        listuimag(1,ie)=gimag;
    
end 

    

    % for system
    h0=em*sz-t0*sx;
    hx=-(d0*s0+d0c*sx);
    Hx=zeros(2*Nx,2*Nx);
    
for ix=1:Nx
        ax=2*(ix-1)+1:2*ix;
        ay=ax+2;
%         Hx(ax,ax)=h0+Vd*Vrd(ix)*s0;
 Hx(ax,ax)=h0+Vd*s0;
        if ay(end)>2*Nx
            continue; % when one acts, Open Boundary Condition
            ay=ay-2*Nx;
        end
        Hx(ax,ay)=hx;
        Hx(ay,ax)=hx';                                        
end
    Yy=eye(2*Nx);
    listden=zeros(iNlen,ielen);    
for ie=1:ielen
        ie;
        en=enmin+(enmax-enmin)*(ie-1)/(ielen-1)+1i*eta;%choose one energy for lead!
        
        Gf=Yy/(en*Yy-Hx);
        listu(1,ie)=real(en);
        gimag=imag(trace(Gf));        
        listden(1,ie)=gimag;         
end

figure;
subplot(3,1,1);
plot(listu,listumean,'.-');
xlabel('Energy'); ylabel('Conductance');
title([' Nx=' num2str(Nx)]);
axis([listu(1) listu(end) 0 1.1])
subplot(3,1,2);
plot(listu,log(abs(listden)),listu,log(abs(listuimag)),'.-');
xlim([listu(1) listu(end) ])

loop=loop+1;

end


if loop==2

% Plot band structure
Nkx=250;
count2=1;kxlist=zeros(Nkx,1);
count=1;Kx=zeros(1,1); Ky=zeros(1,1);En=zeros(1,1);
for ix=1:Nkx
    kx=-pi+2*pi*(ix-1)/Nkx;
    kxlist(ix)=kx;
    
    h0=-2*d0*cos(kx);
    hx=-t0-2*d0c*cos(kx);
    H=h0*s0+ep*s0+em*sz+hx*sx;
      
        dz=eig(H);        
        dzlist(:,ix)=real(dz);
        
end


%figure;
subplot(3,1,3);
plot(dzlist,kxlist,'.-');xlabel('E');ylabel('kx');title([' gam=' num2str(gam)]);
axis([listu(1) listu(end) -pi pi])

end







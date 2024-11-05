clc
clear
close all


%Input parameters for membrane mechanics===========================================================>>>.   
arc_domain =400; % dimensionless area 
mesh=(0:0.0005:1).^4;
lambda = 0.001;  % membrane tension (pN/nm)
R0 = 20; % dimensional length (nm)
gamma =20; % hyp_tan coeff for spon cur 
k0 =84; % bending rigidity 
dk =1;  % k_(p+l)/k_l ratio bending rigidity (do not change for fixed bending rigidity)
dkd=0;  % kd_(p+l)/kd_l-1 ratio of deviatoric bending rigidity(do not change for fixed bending rigidity)
P = 0; %  pressure
dkG = 0; %kg_(p+l)/kg_l-1 ratio of gaussian bending rigidity(do not change for fixed bending rigidity)
alpha=1.0
beta=1.0
%rIn=0;
%rOut=0;
%step=60;
%r1=linspace(1,5900/(2*pi*R0^2),step/3);
%r2=linspace(5900/(2*pi*R0^2),490000/(2*pi*R0^2),step/3);
%r3=492000/(2*pi*R0^2)*ones(1,step/3);
%rF=[r1 r2 r3];

%f0=0;
acoat=20;  %2*600000/(2*pi*R0^2); %%coat area of protein
%c1=zeros(1,2*step/3);
%c2=1*linspace(0,0.004,step/3);
%C0=0.02; %0.0265;  %0.65*0.0075;       %[c1 c2];
%C0=0.65*0.0075;
%D0=0.00;   %0.75*0.0075;
%C0=0.02-D0;
D1=0.0; %0.028;  ## collar deviatoric curvature (if any)
acn=(arc_domain-acoat);
%fc=0.25*7.95/(acn*2*pi*R0^2);
%z1=linspace(0,1200,step/3);
%z2=1200* ones(1,2*step/3);
%zpRng=[z1 z2];
%lamRng=lambda*(5:-1:1);

file_name=['phi.txt'];
fileID2=fopen(file_name,'w+');

file_name=['x.txt'];
fileID3=fopen(file_name,'w+');

file_name=['y.txt'];
fileID4=fopen(file_name,'w+');

file_name=['clv.txt'];
fileID6=fopen(file_name,'w+');

initSol =endoInit(arc_domain,mesh, lambda, k0,R0);

%for ki=1:5
%DRng=linspace(0,0.015,31);

 D0=0; % deviatoric curvature of protein coat (if any)
 ell=0.01431;
    C00=0.18*ell; % spontaneous curvature of protein coat (if any)
kr=0.1;
% phio=1/(1+kr);
    
% Input for kinetic model===========================================================>>>.    
phi=zeros(1,length(mesh)); %+0.005*rand(1,length(mesh));  % initial area fraction of linkers
g=20;

%phi=0.5*(tanh(g*(mesh*arc_domain+ acoat)) - tanh(g*(mesh*arc_domain - acoat)));
tt=0; % initial time
    dt=5e-3*8; % time step (ND)
    kl=0.5; % linkers stiffness (pN/nm)
    sig0=1e-4;  % sat density of linkers (#/nm^2)
    fc=kl*sig0;  % force applied by liners on membrane per unit area per unit deformation 
    kBT=4.14; %Boltzman_const X temp
    tm=0;  % 
    del=1; % standard bond length of dimers (nm)
    %C0=0.0023436;
    %k1=8*2.09;
    %k2=15;
    %k3=8*1.76;
    eta(1)=0.18;
    C0=eta(1)*ell;
    k1=1.75;
    k3=0.025;
    H=ones(1,length(mesh))
for i=2:501

   
  tm=tm+dt; % local time
  Sl=0.1*(k1*exp(alpha*H(i))*(1-eta(i-1)))-k3*exp(beta*H(i))*eta(i-1);
    eta(i)=eta(i-1)+Sl*(dt);
    C0=eta(i)*ell;
   %C0=C0+(k1*(0.02-C0)-k3*C0)*dt;
    %C0=C00*(tanh(10*tm)); % spontaneous curvature at local time
    
Nn=size(mesh);       % mesh size
Np=Nn(2); % number of mesh points
%phi=0.5*(1+ tanh(20*(arc_domain*mesh - acoat)));
    
try
[t, Sol]= bud_sage_microparticle_dyn(arc_domain,mesh, lambda, fc, phi, acoat,k0, dk,dkd,dkG, P, gamma, C0, D0, D1, R0, tt, initSol);
catch ME
        
        display(ME.message);
        display(sprintf('Error solution: \\lambda = %0.5f D= %0.3f', C0,D0))
        Sol=initSol;
        
%initSOl=Sol;

H=Sol(4,:);

end

initSol=Sol;


fprintf(fileID2,'%16.8f',phi(:)); 
fprintf(fileID2,'\n');

fprintf(fileID3,'%16.8f',Sol(1,:)); 
fprintf(fileID3,'\n');


fprintf(fileID4,'%16.8f',Sol(2,:)); 
fprintf(fileID4,'\n');

yy=Sol(2,:);
fprintf(fileID6,'%16.8f',C0,yy(1),max(Sol(3,:))); 
fprintf(fileID6,'\n');

%max(Sol(3,:))



    fs=-Sol(2,:)*R0; % dimensional deformation
    for j=1:length(mesh)
     phiS(j)=0.01*((0.9999-phi(j))/kr-exp(kl*del*fs(j)/kBT)*phi(j))/8;    % net binding rate k_on(1-phi)-k_off*exp(f*del/k_BT)*phi
    if(abs(Sol(2,j)*R0)>25)   %% checking for realistic value of density
    phi(j)=0;
    elseif (phi(j)+phiS(j)*dt<0)
    phi(j)=0;
    elseif (phi(j)+phiS(j)*dt>1)
    phi(j)=1;
    else
    phi(j)=phi(j)+phiS(j)*dt;  % linker density update
    end
    end

    tt=tt+dt;
  %  fc=fc*1.1
%Nn=size(Sol);
%N=Nn(2);
%xx(N)=Sol(1,1);
%yy(N)=Sol(2,1);

%for i=1:N-1
   %xx(N-i)=-Sol(1,i+1);
   %xx(N+i)=Sol(1,i+1);
  % yy(N-i)=Sol(2,i+1);
 %  yy(N+i)=Sol(2,i+1);
%end

%plot(xx,yy)

%pause(1)
%end


end

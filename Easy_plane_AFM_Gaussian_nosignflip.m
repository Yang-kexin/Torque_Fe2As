clear all
close all
clc

%syms xi
chipp=0.01; % dimensionless
chipr=0.0001;
psi=0;
%H0=15000; % [Oe]
%m=0.00128; % mass 
%Mmol=545; % molar mass
chix=[chipr,0,0;0,chipp,0;0,0,chipp]; %moment along x-axis
chiy=[chipp,0,0;0,chipr,0;0,0,chipp]; %moment along y-axis

K1=-2e5;
K2=-2e3;
K22=150;%[J/m^3]
K44=0;
%N=20;
%sigma=150;
%L=sigma*2;
%tau_nor=zeros(N+1,200);
%nor(1)=0;
% for j=1:1:N;
%     delta(j)=-L+2*L*j/N;
%     K22(j)=K4+delta(j);
for j=1:1:5;
    H0=j*1000+10000;
for i=1:1:200;
xi=i*2*pi/200;
H=H0*[cos(xi);sin(xi);0];
% Hc=sqrt(2*K22/(chipp-chipr));
% Hx=H(1);
% Hy=H(2);
% N_1(i)=0.5*exp(-2*(abs(Hx)/6568)^2);
% N_2(i)=0.5*exp(-2*(abs(Hy)/6568)^2);
% %N_parallel=0.5;
% N_x(i)=N_1(i)/(N_2(i)+N_1(i))*1;
% N_y(i)=1-N_x(i);
X0y=[0];%theta,omega
%Xsoly=fminsearch(@(X2)Etot2(X2,xi,H0,K1,K2,K22,K44,chiy,chipp),X0y);
Xsoly=fminbnd(@(X2)Etot2(X2,xi,H0,K22,K44,chix,chiy),-pi,pi);
thetay=Xsoly(1);
%phiy=Xsoly(2)
X2=Xsoly;
[Ey]=Etot2(X2,xi,H0,K22,K44,chix,chiy);
Rz=[cos(thetay),-sin(thetay),0;sin(thetay),cos(thetay),0;0,0,1];
R=Rz;
Ei(i)=Ey;
Hx=H0*cos(thetay-xi);
Hy=H0*sin(thetay-xi); 
N_1=0.5*exp(-2*(abs(Hx)/6568)^2);
N_2=0.5*exp(-2*(abs(Hy)/6568)^2);
%N_parallel=0.5;
N_x=N_1/(N_2+N_1)*1;
N_y=1-N_x;
chi=chix*N_x+chiy*N_y;
Ry=[cos(thetay),-sin(thetay),0;sin(thetay),cos(thetay),0;0,0,1];
thetax=thetay+pi/2;
Rx=Ry;
%Rx=[cos(thetax),-sin(thetax),0;sin(thetax),cos(thetax),0;0,0,1];
% R=Rz;
chiprimex=Rx*(chix*N_x)*inv(Rx);
chiprimey=Ry*(chiy*N_y)*inv(Ry);
Mx=chiprimex*H; % Oe->A/m
My=chiprimey*H;
taux=cross(Mx,H)*0.1/4/pi;
tauy=cross(My,H)*0.1/4/pi; %N/m or J/m^3; T*A/m=J/m^3
tau_tot(i)=taux(3)+tauy(3);
theta1(i)=thetay;
xaxis(i)=xi; % angle between easy axis and field
Nxi(i)=N_x;
Nyi(i)=N_y;
end
tau_j(j,:)=tau_tot(1,:);
figure(1)
plot(xaxis,tau_j(j,:))
hold on
end

% p(j)=exp(-(delta(j))^2/2/sigma^2)/sqrt(2*pi*sigma^2);
% tau_nor(j+1,:)=tau_nor(j,:)+tau_tot*p(j);
% tau_totj(j)=tau_tot(43);
% nor(j+1)=nor(j)+p(j);
% jaxis(j)=j;
%  end


%  tau_tot=diff(Ei);
%  x=xaxis(1:199);
% 
% figure(1)
% plot(H1,A)
% figure(3)
% plot(xaxis/pi*180,Ei)
% xlabel('angle of field')
% ylabel('energy')
% legend('total')


% figure(2)
% plot(delta,p)
% xlabel('x-u')
% ylabel('probablity')
% legend('domain x','domain y')
% 
% figure(4)
% plot(delta,tau_totj)
% xlabel('K-K0')
% ylabel('torque each iteration')
% legend('domain x','domain y')

% figure(1)
% plot(xaxis/pi*180,tau_nor(N+1,:)/nor(N+1)) 
% xlabel('angle of field')
% ylabel('Torque (N/m^2)')
% legend('domain y')

% figure(4)
% plot(xaxis,tau_tot)
% xlabel('angle of field')
% ylabel('Torque (N/m^2)')
% 
% figure(3)
% plot(xaxis,theta1)
% xlabel('angle of field')
% ylabel('theta')
% legend('domain x')
% % 
% figure(2)
% plot(xaxis,Nxi,xaxis,Nyi)
% xlabel('angle of field')
% ylabel('N')
% legend('domain x','domain y')
% 
% figure(4)
% plot(xaxis,Eyi)
% xlabel('angle of field')
% ylabel('Energy')
% legend('y')
%figure(5)
% plot(xaxis,moment2)
% xlabel('angle of field')
% ylabel('theta')
% legend('domain y')

% figure(4)
% plot(xaxis,omegai)
% xlabel('angle of field')
% ylabel('Angle between two sublattice')
% 1
% Ry=[cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
% Rz=[cos(phi),-sin(phi),0;sin(phi),cos(phi),0;0,0,1];
% R=Rz*Ry;
% 
% chi2=R*chi*inv(R);
% H=H0*[cos(psi)*sin(xi);sin(psi)*sin(xi);cos(xi)];
% Fz=-1/2*m/Mmol*transpose(H)*chi2*H; % Zeemann energy
% Fa=m/Mmol*K1*(sin(theta))^2; % anisotropy energy
% Etot=Fz+Fa;
% Etot



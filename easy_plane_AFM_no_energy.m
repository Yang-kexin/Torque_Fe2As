close all
clear all
clc

%H0=5000;
chipp=0.015; % chi_perpendicular
chipr=0.000;  % chi_parallel
chiprime = chipp*1.05;  % chi_perpendicular_prime
for j=1:1:4;
H0=20000*j+10000;
for i=1:1:100;
xi=i*2*pi/100;
H=H0*[sin(xi);0;cos(xi)];
Hx=H(1);
N_parallel(i)=0.5*exp(-2*(abs(Hx)/7000)^2);
N_perpendicular(i)=1-N_parallel(i);
% torque of domain along x axis
chi1=[chipr,0,0;0,chiprime,0;0,0,chipp];
m1=chi1*H;
tau1=-cross(m1,H); % Oe^2
tau1y(i)=tau1(2)*0.1/4/pi*N_parallel(i);
% torque of domain along y axis
chi2=[chiprime,0,0;0,chipr,0;0,0,chipp];
m2=chi2*H;
tau2=-cross(m2,H);
tau2y(i)=tau2(2)*0.1/4/pi*N_perpendicular(i);
tau_tot(j,i)=tau1y(i)+tau2y(i);
xaxis(i)=xi;
angle(i)=acos((m2(1)*H(1)+m2(2)*H(2)+m2(3)*H(3))/norm(m2)/norm(H))*180/pi;
end
tau1yj(j,:)=tau1y(1,:);
tau2yj(j,:)=tau2y(1,:);
figure(1)
plot(xaxis*180/pi,tau2y,'LineWidth',1.5)
xlabel('Sample position ')
ylabel('Torque of domains along y-direction')
legend('0.5 T','1 T','1.5 T','2 T','2.5 T','3 T')
hold on
figure(2)
plot(xaxis*180/pi,tau_tot(j,:),'LineWidth',1.5)
xlabel('Sample position ')
ylabel('Sum of torque of two domains')
legend('0.5 T','1 T','1.5 T','2 T','2.5 T','3 T')
hold on
% figure(3)
% plot(xaxis*180/pi,tau1y,'LineWidth',1.5)
% xlabel('Sample position ')
% ylabel('Torque of domains along x-direction')
% legend('0.5 T','1 T','1.5 T','2 T','2.5 T','3 T')
% hold on
end
% figure(3)
% plot(xaxis,N_parallel,xaxis,N_perpendicular)
% figure(1)
% plot(xaxis,tau1y,xaxis,tau2y,xaxis,tau_tot)
% legend('domain along x-axis','domain along y-axis','sum of two domains')
% figure(2)
% plot(xaxis,tau_tot)


function [E]=Etot2(X,xi,H0,K22,K44,chix,chiy)
psi=0;
mu0=1.257E-6; %{T*m/A=N/A^2}
sbM=4e5; %[A/m] sublattice magnetization
omega=0;
theta=X(1);
Ry=[cos(theta),-sin(theta),0;sin(theta),cos(theta),0;0,0,1];
Rx=Ry;
%Rx=[-sin(theta),-cos(theta),0;cos(theta),-sin(theta),0;0,0,1];
H=H0*[cos(xi);sin(xi);0];
Hx=H0*cos(theta-xi);
Hy=H0*sin(theta-xi); 
N_1=0.5*exp(-2*(abs(Hx)/6568)^2);
N_2=0.5*exp(-2*(abs(Hy)/6568)^2);
%N_parallel=0.5;
N_x=N_1/(N_2+N_1)*1;
N_y=1-N_x;
chi2=Rx*(chix*N_x)*inv(Rx)+Ry*(chiy*N_y)*inv(Ry);
E=K22*cos(4*theta)+K44*cos(8*theta)-0.5*transpose(H)*chi2*H*0.1/4/pi;%*(cos(theta-xi))^2;
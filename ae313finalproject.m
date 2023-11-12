% -------------------------------------------
% Program Name: final project option 2
% Program Purpose:^^
%
% Created By: Lee Weidow
%
% Created On: nov 25 2021
% Last Modified On:
%
% Comments: code used is taken from codes used for hw 4 and 5
% hardcoded for the values given for the project
% unused inputs for variables are for my own use in the future
% -------------------------------------------

clc
clear
close all

mu=398600.44;
Re=6378.1366; %km
K=[0 0 1];

%finding initial values

%phi=input('Phi in deg: ');
phi=32.248814; 
phi=phi*pi/180;
%lambda=input('Lambda in Deg: ');
lambda=36-74.99; 
lambda=lambda*pi/180; %deg to rad
%rho=input('Rho in km : '); %km
rho=822; 
%B=input('Beta in Deg: ');
B=18.0912; 
B=B*pi/180;
%sigma=input('Sigma in Deg: ');
sigma=61.7066; 
sigma=sigma*pi/180;
%rhodot=input('Rho dot in km/s: '); %km/s
rhodot=3.48499169;
%Bdot=input('Beta dot in deg/s: ');
Bdot= .269604966; 
Bdot=Bdot*pi/180; %rad/s
%sigmadot=input('Sigma dot in deg: ');
sigmadot=-.4321605433; 
sigmadot=sigmadot*pi/180;
%TOF=input('Time of flight in seconds: ');
%TOF=32*60; %case 1
TOF=10*3600; %case 2
H=[0 0 0]';
H(3)=0; %input('Ground station elevation in km: ') assuming the ground station is at sea level
we=[0 0 7.2921*10^-5]'; %rad/s rotational speed of earth ECI
rsitesez=[0 0 6378]'+H; %km

%finding ri and vi

rhosez=[rho*-cos(sigma)*cos(B) rho*cos(sigma)*sin(B) rho*sin(sigma)]';
rot=[cos(lambda) -sin(lambda) 0; sin(lambda) cos(lambda) 0; 0 0 1]*[sin(phi) 0 cos(phi); 0 1 0; -cos(phi) 0 sin(phi)];
rhoeci=rot*rhosez;
rsiteeci=rot*rsitesez;
rvececi=rhoeci+rsiteeci; %this and vvececi are initial values
rhodots=-rhodot*cos(sigma)*cos(B)+rho*sigmadot*sin(sigma)*cos(B)+rho*Bdot*cos(sigma)*sin(B);
rhodote=rhodot*cos(sigma)*sin(B)-rho*sigmadot*sin(sigma)*sin(B)+rho*Bdot*cos(sigma)*cos(B);
rhodotz=rhodot*sin(sigma)+rho*sigmadot*cos(sigma);
rhodotsez=[rhodots rhodote rhodotz]';
rhodoteci=rot*rhodotsez;
vvececi=rhodoteci+cross(we,rvececi);

%finding initial OEs

r=norm(rvececi);
v=norm(vvececi);
hvec=cross(rvececi,vvececi);
h=norm(hvec);
energy=.5*v^2-mu/r;
a=-mu/(2*energy);
%p=h^2/mu; %not needed here but im keeping it just incase
hhat=hvec./h;
i=acos(dot(hhat,K));
rhat=rvececi./r;
evec=cross(vvececi,hvec)/mu-rhat;
e=norm(evec);
ehat=evec./e;
if dot(rvececi,vvececi)>=0
    f=acos(dot(ehat,rhat));
elseif dot(rvececi,vvececi)<0
    f=2*pi-acos(dot(ehat,rhat));
end
nhat=cross(K,hvec)./norm(cross(K,hvec));
if dot(ehat,K)>=0
    w=acos(dot(nhat,ehat));
elseif dot(ehat,K)<0
    w=2*pi-acos(dot(nhat,ehat));
end
omega=atan2(nhat(2),nhat(1));

%finding final OEs w/ tof

J2=1.08264*10^-3;
T=2*pi*sqrt(a^3/mu); %time/orbit
n=2*pi/T;
omegadot=-(3*n*J2)/(2*(1-e^2)^2)*(Re/a)^2*cos(i);
wdot=-(3*n*J2)/(4*(1-e^2)^2)*(Re/a)^2*(5*(cos(i))^2-1);
omegaf=omega+omegadot*TOF;
wf=w+wdot*TOF;

%using OEtoRV
p=a*(1-e^2);
if p==0
    p=a;
end
h2=sqrt(p*mu);

%newtons method to find f2

M=n*TOF;
E1=M;
x=1; %place holder equal to E2-E1
while x>.05 %small enough value for this 
    E2=E1-(E1-e*sin(E1)-M)/(1-e*cos(E1));
    x=E2-E1;
    E1=E2;
end
f2=2*atan(sqrt((1+e)/(1-e))*tan(E2/2));

%using f2 to find vf and rf

r2=p/(1+e*cos(f2));
rpqw=[r2*cos(f2) r2*sin(f) 0]';
vpqw=(mu/h).*[-sin(f2) e+cos(f2) 0];
vpqw=vpqw';
v=norm(vpqw);
rot=[cos(omegaf) -sin(omegaf) 0; sin(omegaf) cos(omegaf) 0; 0 0 1]*[1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)]*[cos(wf) -sin(wf) 0; sin(wf) cos(wf) 0;0 0 1];
%for reference:
%rot1=[cos(omega) -sin(omega) 0; sin(omega) cos(omega) 0; 0 0 1];
%rot2=[1 0 0; 0 cos(i) -sin(i); 0 sin(i) cos(i)];
%rot3=[cos(w) -sin(w) 0; sin(w) cos(w) 0;0 0 1];
vvecf=rot*vpqw;
rvecf=rot*rpqw;

%outputs

f=f*180/pi; %rad to deg
f2=f2*180/pi;
i=i*180/pi;
w=w*180/pi;
wf=wf*180/pi;
omega=omega*180/pi;
omegaf=omegaf*180/pi;
fprintf('Position and Velocity Vectors:\n')
fprintf('r0=%.2fI + %.2fJ + %.2fK (km)\n',rvececi(1),rvececi(2),rvececi(3))
fprintf('v0=%.2fI + %.2fJ + %.4fK (km/s)\n',vvececi(1),vvececi(2),vvececi(3))
fprintf('rf=%.2fI + %.2fJ + %.2fK (km)\n',rvecf(1),rvecf(2),rvecf(3))
fprintf('vf=%.2fI + %.2fJ + %.4fK (km/s)\n\n',vvecf(1),vvecf(2),vvecf(3))
fprintf('Initial OEs:\na=%f (km)\nf=%f (deg)\ne=%f\ni=%f (deg)\nw=%f (deg)\nomega=%f (deg)\n\n',a,f,e,i,w,omega)
fprintf('Final OEs:\na=%f (km)\nf=%f (deg)\ne=%f\ni=%f (deg)\nw=%f (deg)\nomega=%f (deg)\n\n',a,f2,e,i,wf,omegaf)

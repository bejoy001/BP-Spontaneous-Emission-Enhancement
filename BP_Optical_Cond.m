%%% Optical conductivity calculation for few-layer black phosphorus
%%% William S. Whitney (wwhitney@caltech.edu)
%%% Based on calculation by Tony Low, et al., in doi:10.1103/PhysRevB.90.075434
%%% Please read that paper for details on the physics, including
%%% approximations and constants, and normalization of the output
% Physical constants
clc;
clear all;
close all;
tic
delete(gcp('nocreate'));
h = 4.135668e-15; % planck's constant in ev
hb = h/(2*pi); % reduced planck's constant in ev
hbj = 1.054571e-34; % reduced planck's constant in joules
c = 2.9979e8; % speed of light
e = 1.6021766e-19; % elementary charge, in coulombs
m0 = 9.10938e-31; % electron mass, in kg
eps0 = 8.854187e-12; % permittivity of free space
% Parameters to change. User must set these, and data / plotting file names at end
ml = 5.3*1e-10;
n = 38; % for mono layer
tz = n*ml; % flake thickness
% bandgap calculation
A = 1.5180;
B = 0.3;
b = -0.6861;
Eg = B+A*n^b;
axis = 'xx'; % crystal / polarization axis ('xx' or 'yy')
egap = Eg; % band gap
T = 300; % temperature, in k
%omega = logspace(13, 17, 100);
%aryhbwj = hbj*omega; % frequency band for calculation, in ev
aryhbwj = linspace(0.4*e, 1*e, 50);
kpartitions = 1000; % points to sample in k space (numerical instabilities below 1000 or so)
intsubbands = 10; % number of subbands to sum over (too many = slow, too few = large inaccuracies)
pp = parpool(8); % set up parallel computing resource (machine dependent, read appropriate manpages)
% Setup
a = 4.47e-10/2; % using low / lin convention
b = 3.34e-10/2; % using low / lin convention
Ev = -0.3*e; % bulk band energies, set Ec to 0
Ec = 0*e;
etac = (hbj^2)/(0.4*m0);
etav = (hbj^2)/(0.4*m0);
nuc = (hbj^2)/(1.4*m0);
nuv = (hbj^2)/(2.0*m0);
gamma = e*4*a/pi;
beta = e*2*(a^2)/(pi^2);
etad = 5e-3*e;
mvz = 0.46*m0; % valence eff mass, from Lin, NL 2016
mcz = 0.25*m0; % cond eff mass, from Lin, NL 2016
dgap = egap - ((Ec + hbj*hbj*pi*pi/(2*mcz*tz*tz))-(Ev - hbj*hbj*pi*pi/(2*mvz*tz*tz)))/e;
dsc = 1*e*dgap/2;
dsv = 1*e*dgap/2;
deltaEv1 = (hbj^2)*(pi^2)/(2*mvz*(tz^2)) + dsv;
deltaEc1 = (hbj^2)*(pi^2)/(2*mcz*(tz^2)) + dsc;
Ec1 = Ec + deltaEc1;
Ev1 = Ev - deltaEv1;
Ef = (Ec1 + Ev1)/2+0.1*e;
kbsi = 1.38064852e-23;
kT = kbsi*T;
aryw = aryhbwj./hbj;
arykx = linspace(-pi/a,pi/a,kpartitions);
aryky = linspace(-pi/b,pi/b,kpartitions);
arys = [1 -1];
arysp = [-1 1];
aryj = 1:intsubbands;
arysigxx = zeros(1,length(aryw));
% Start calculating
% Loop over frequencies
for idxw = 1:length(aryw)
fltSum = 0;
hbwj = aryhbwj(idxw);
fltPctDone = 100*idxw/length(aryw);
disp(sprintf('%0.1f%% Complete\n',fltPctDone));
% Loop over s quantum number
for idxs = 1:length(arys)
% Loop over j quantum number
for idxj = 1:length(aryj)
s = arys(idxs);
sp = arysp(idxs);
j = aryj(idxj);
deltaEvj = (j^2)*(hbj^2)*(pi^2)/(2*mvz*(tz^2)) + dsv;
deltaEcj = (j^2)*(hbj^2)*(pi^2)/(2*mcz*(tz^2)) + dsc;
Ecj = Ec + deltaEcj;
Evj = Ev - deltaEvj;
delj = Ecj - Evj;
aryIntValskx = zeros(1,length(arykx));
% Loop over kx
parfor idxkx = 1:length(arykx) % this loop is parallelized, make regular for to debug
aryIntValsky = zeros(1,length(aryky));
kx = arykx(idxkx);
% Loop over ky
for idxky = 1:length(aryky)
ky = aryky(idxky);
Hj = [Ecj + etac*kx^2 + nuc*ky^2, gamma*kx + beta*ky^2;
gamma*kx + beta*ky^2, Evj - etav*kx^2 - nuv*ky^2];
if axis == 'xx'
v = (1/hbj)*[2*etac.*kx gamma; gamma -2*etav.*kx];
elseif axis == 'yy'
v = (1/hbj)*[2*nuc.*ky 2*beta*ky; 2*beta*ky -2*nuv.*ky];
else
error('Bad axis specification!');
end
if (s == 1)
Esjk = (1/2)*((Ecj + Evj) + (etac-etav)*kx^2 + (nuc-nuv)*ky^2) ...
+ (1/2)*((delj^2 + delj*(2*(etac+etav)*kx^2 + 2*(nuc+nuv)*ky^2)) ...
+ ((etac+etav)*kx^2 + (nuc+nuv)*ky^2)^2 + 4*(gamma*kx + beta*ky^2)^2)^(1/2);
Phisjk = [ 1; (Esjk-Hj(1,1))/Hj(1,2)];
Phisjk = Phisjk/norm(Phisjk);
else
Esjk = (1/2)*((Ecj + Evj) + (etac-etav)*kx^2 + (nuc-nuv)*ky^2) ...
- (1/2)*((delj^2 + delj*(2*(etac+etav)*kx^2 + 2*(nuc+nuv)*ky^2)) ...
+ ((etac+etav)*kx^2 + (nuc+nuv)*ky^2)^2 + 4*(gamma*kx + beta*ky^2)^2)^(1/2);
Phisjk = [ 1; (Esjk-Hj(1,1))/Hj(1,2)];
Phisjk = Phisjk/norm(Phisjk);
end
if (sp == 1)
Espjk = (1/2)*((Ecj + Evj) + (etac-etav)*kx^2 + (nuc-nuv)*ky^2) ...
+ (1/2)*((delj^2 + delj*(2*(etac+etav)*kx^2 + 2*(nuc+nuv)*ky^2)) ...
+ ((etac+etav)*kx^2 + (nuc+nuv)*ky^2)^2 + 4*(gamma*kx + beta*ky^2)^2)^(1/2);
Phispjk = [ 1; (Espjk-Hj(1,1))/Hj(1,2)];
Phispjk = Phispjk/norm(Phispjk);
else
Espjk = (1/2)*((Ecj + Evj) + (etac-etav)*kx^2 + (nuc-nuv)*ky^2) ...
- (1/2)*((delj^2 + delj*(2*(etac+etav)*kx^2 + 2*(nuc+nuv)*ky^2)) ...
+ ((etac+etav)*kx^2 + (nuc+nuv)*ky^2)^2 + 4*(gamma*kx + beta*ky^2)^2)^(1/2);
Phispjk = [ 1; (Espjk-Hj(1,1))/Hj(1,2)];
Phispjk = Phispjk/norm(Phispjk);
end
fEsjk = 1/(exp((Esjk-Ef)/kT)+1);
fEspjkp = 1/(exp((Espjk-Ef)/kT)+1);
aryIntValsky(idxky) = ((fEsjk - fEspjkp)/(Esjk - Espjk)) ...
*(transpose(Phisjk)*v*Phispjk)*(transpose(Phispjk)*v*Phisjk)/(Esjk-Espjk+hbwj+etad*1i);
end
aryIntValskx(idxkx) = trapz(aryky,aryIntValsky);
end
fltIntVal = -1i*2*hbj*e*e*trapz(arykx,aryIntValskx)/(4*pi*pi);
fltSum = fltSum + fltIntVal;
end
end
arysigxx(idxw) = fltSum/(e*e/(4*hbj));
end
% Generate plots and save data. User must modify file names below
toc
clear all;
clc;
close all;

tic
h = 6.626*1e-34;  % Planck's constant
hcut = h/(2*pi); 
e = 1.6*1e-19;	% electron charge


m0 = 9.1*1e-31;	%electron mass	
mu0 = 4*pi*1e-7;	
epsilon0 = 8.854*1e-12;

eta0 = sqrt(mu0/epsilon0);
omega = logpsace(13, 17, 100);
c = 3*1e8;
f = omega/(2*pi);
lambda = c./f;		% operating wavelengths

K0=(2*pi)./lambda;

N = 2*1e3;		% number of points for intergation
kmin = -500;	
kmax = 500;


a = [1e-9, 1e-8, 1e-7]; 

mu = [1 0 0];
mux = mu(1);
muy = mu(2);
muz = mu(3);
mu0 = 4*pi*1e-7;
epsilon0 = 8.854*1e-12;

load cond_xx.mat;	% optical conductivity along x axis (AC) of BP
load cond_yy.mat;   % optical conductivity along y axis (ZZ) of BP
sigma_0=e^2/(4*hcut);	
SigmaX=sigma_0*arysigxx;
SigmaY=sigma_0*arysigyy;


PF=zeros(length(a),length(f));
for k=1:length(a)
    for j = 1:length(f)
        
        internal_int = zeros(1, N);
        sigmaX=SigmaX(j);
        sigmaY=SigmaY(j);
        z0=a(k);
        k0=K0(j);
        kx= linspace(kmin*k0, kmax*k0, N);
        ky = linspace(kmin*k0, kmax*k0, N);
        parfor i = 1:length(kx)
            krho = sqrt(kx(i).^2+ky.^2);
            kz = sqrt(k0.^2-krho.^2);
            
            sigmaxx = (kx(i).^2*sigmaX+ky.^2*sigmaY)./krho.^2;
            sigmaxy = (kx(i).*ky.*(sigmaY-sigmaX))./krho.^2; %
            sigmayx = sigmaxy;
            sigmayy= (ky.^2*sigmaX+kx(i).^2*sigmaY)./krho.^2;
            
            
            Zs = kz./k0;
            Zp = k0./kz;
            cp = kz./k0;
            den_r = (2*Zs+eta0.*sigmayy).*(2*Zp+eta0.*sigmaxx)-eta0^2.*sigmaxy.*sigmayx;
            Rss = (-eta0.*sigmayy.*(2*Zp+eta0.*sigmaxx)+eta0.^2.*sigmaxy.*sigmayx)./den_r;
            Rsp = -(2.*cp.*Zp.*eta0.*sigmaxy)./den_r;
            Rps = (2.*Zs.*eta0.*sigmayx)./(cp.*den_r);
            Rpp = (eta0.*sigmaxx.*(2.*Zs+eta0.*sigmayy)+eta0.^2.*sigmaxy.*sigmayx)./den_r;
            
            M11 = (Rss.*(ky.^2))./(kz.*krho.^2)+((-kx(i).*ky).*Rsp)./(k0.*krho.^2)+((kx(i).*ky).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i)^2.*kz))./(k0.^2.*krho.^2);
            M12 = (Rss.*(-kx(i).*ky))./(kz.*krho.^2)+((-ky.^2).*Rsp)./(k0.*krho.^2)+((-kx(i).^2).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i).*ky.*kz))./(k0^2.*krho.^2);
            M13 = (((-ky.*krho.^2)./kz).*Rsp)./(k0.*krho.^2)-(Rpp.*((kx(i).*krho.^2.*kz)./kz))./(k0^2.*krho.^2);
            M21 = (Rss.*(-kx(i).*ky))./(kz.*krho.^2)+((kx(i).^2).*Rsp)./(k0.*krho.^2)+((ky.^2).*Rps)./(k0.*krho.^2)-(Rpp.*(kx(i).*ky.*kz))./(k0^2.*krho.^2);
            M22 = (Rss.*(kx(i).^2))./(kz.*krho.^2)+((kx(i).*ky).*Rsp)./(k0.*krho.^2)-((kx(i).*ky).*Rps)./(k0.*krho.^2)-(Rpp.*(ky.^2.*kz))./(k0^2.*krho.^2);
            M23 = (((kx(i).*krho.^2)./kz).*Rsp)./(k0.*krho.^2)-(Rpp.*((ky.*krho.^2)))./(k0^2.*krho.^2);
            M31 = -(((ky.*krho.^2)./kz).*Rps)./(k0.*krho.^2)+(Rpp.*(kx(i).*krho.^2))./(k0^2.*krho.^2);
            M32 = (((kx(i).*krho.^2)./kz).*Rps)./(k0.*krho.^2)+(Rpp.*(ky.*krho.^2))./(k0^2.*krho.^2);
            M33 = (Rpp.*(krho.^2))./(k0^2.*kz);
            fun = (mux.*(mux.*M11+muy.*M21+muz.*M31)+muy.*(mux.*M12+muy.*M22+muz.*M32)+muz.*(mux.*M13+muy.*M23+muz.*M33));
            integrand = fun.*exp(1i*2.*kz*z0);
            internal_int(i) = trapz(ky, integrand);
            
        end
        
        gamma = 1i*trapz(kx, internal_int)/(8*pi^2);
        PF(k,j) = 1+(6*pi/K0(j))*imag(gamma);       
        
    end
end



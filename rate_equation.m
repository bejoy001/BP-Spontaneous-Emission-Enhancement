clc;
close all;
clear all;


load PF.mat; % Purcell Factors for different wavelength; can be obtained from PurcellFactor.m code
load a.mat;	 % distance values
omega_n=logspace(14.9,16.32,100); % desired frequencies
Pin=zeros(length(a),length(omega_n));

for i = 1:length(a)

    n = 8;
    G = zeros(8,1);
    G(1) = -1e4;
    G(8) = 1e4;
    A = zeros(n,n);
    A(1,3) = 6.2649*1e8*PF(i,4); % multiply purcell factor of 122 nm with this
    A(1,5) = 1.6725*1e8*PF(i,5);%  multiply purcell factor of 103 nm with this
    A(2,5) = 2.2448*1e7*PF(i,2);%  multiply purcell factor of 656 nm with this
    A(3,3) = -6.2649*1e8*PF(i,4);% multiply purcell factor of 122 nm with this
    A(3,8) = 2.0625*1e7*PF(i,3); % multiply purcell factor of 486 with this
    A(5,5) = -(2.2448*1e7*PF(i,2)+1.6725*1e8*PF(i,5));% multiply purcell factor of 656 nm and 103 respectively
    A(5,8) = 7.0376*1e6*PF(i,1); % multiply purcell factor of 1875 nm with this
    A(8,8) = -(7.0376*1e6*PF(i,1)+2.0625*1e7*PF(i,3));% multiply purcell factor of 1875 nm and 486 nm respectively
    A = -A;
    N = pinv(A)*G;

    lambda = [1875 656 486 122 103]*1e-9;
    c = 3*1e8;
    omega = (2*pi*c)./lambda;
    hcut = (6.636*1e-34)/(2*pi);

    P1875 = hcut*omega(1)*N(8)*7.0376*1e6*PF(i,1);% multiply purcell factor of 1875 nm with this
    P656 = hcut*omega(2)*N(5)*2.2448*1e7*PF(i,2);% multiply purcell factor of 656 nm with this
    P486 = hcut*omega(3)*N(8)*2.0625*1e7*PF(i,3);% multiply purcell factor of 486 nm with this
    P122 = hcut*omega(4)*N(3)*6.2649*1e8*PF(i,4);% multiply purcell factor of 122 nm with this
    P103 = hcut*omega(5)*N(5)*1.6725*1e8*PF(i,5);% multiply purcell factor of 103 nm with this
    P = [P1875 P656 P486 P122 P103];

    Pin_n=zeros(1,length(omega_n));

    for k=1:length(omega)
        [c ind]=min(abs(omega_n-omega(k)));
        omega_n(ind)=omega(k);
        Pin_n(ind)=P(k);
        dx=omega_n(ind)-omega_n(ind-1);
        omega_n(ind+1)=dx+omega_n(ind);
    end


     % normailzation
    Pin(i,:)=Pin_n/sum(Pin_n);
   
end

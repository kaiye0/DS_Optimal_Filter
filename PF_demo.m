%========================================================================%
%                          PF_demo.m.
%
%   This is an implementation of a particle filter on a modified
%   form of Example 15.1 in Simon.
%
%   Demoz Gebre-Egziabher
%   Last Modified 12/04/2019
%
%=========================================================================%

clear all;
close all;
clc;

%   Defien variables and place holders

drl = 101;                  %   Data record length
t = zeros(drl,1);           %   Time vector
x = zeros(drl,1);   
y = zeros(drl,1);           %   Measurement history (true)
ym = zeros(drl,1);          %   Measurement history (from sensors)
dt = 1;                     %   Time step
sigma_v = 0.02;            %   Standard deviation of noise v
sigma_w = 0.2;
a = 0.1;
b = 2;

%   Generate truth data

y(1) = b*x(1);
ym(1) = y(1) + (sigma_v*randn(1))^2;

for k = 2:drl
    
    t(k) = t(k-1) + dt;
    x(k) = a*x(k-1) + sigma_w*randn(1);
    y(k) = b*x(k);
    ym(k) = y(k) + (sigma_v*randn(1))^2;
    
end

figure(gcf)
subplot(211)
plot(t,y,'b-',t,ym,'ro');grid on;ylabel('y');xlabel('Time (sec)');
title('Truth and Measured Values of y')
subplot(212)
plot(t,ym-y,'b-');grid on;ylabel('\delta y');xlabel('Time (sec)');

    
%   Particle Filter (PF) Implementation

N = 10000;                               %   Number of particles
x_p =  zeros(N,drl);
x_pf = zeros(drl,1);
x_p(:,1) = 2*rand(N,1) - 1;
P_pf = zeros(drl,1);
x_pf(1) = mean(x_p(:,1));
P_pf(1) = (std(x_p(:,1)))^2;
figure
hist(x_p(:,1),50);       %   Plot histogram of initial sample of particles

wB = waitbar(0,'Particle Filter Simulations in Progress....');

for k = 2:drl
    
    waitbar(k/drl,wB);
    
    x_p(:,k) = a*x_p(:,k-1) + sigma_w*randn(1);
    y_pf = b*x_p(:,k);
    q_tilde = my_chisqr_pdf(sigma_v^2*(ym(k)*ones(N,1) - y_pf),1);
    q = q_tilde/(sum(q_tilde));             %   Normalized liklihoods
    
    for m = 1:N
        r = rand(1,1);
        s_q = 0;
        n = 0;
        while (s_q < r)
            n = n + 1;
            s_q = s_q + q(n);
        end
        x_p(m,k) = x_p(n,k);
    end
    P_pf(k) = var(x_p(:,k));
    x_pf(k) = mean(x_p(:,k));
end   

close(wB);

%   Plot Partilce Progression
%
% figure(gcf+1)
% subplot(131);hist(x_p(:,1),50);subplot(132);hist(x_p(:,5),50);subplot(133);hist(x_p(:,10),50);
figure
subplot(211)
plot(t,x,'b-',t,x_pf,'r--');grid on;
title('PF Performance');ylabel('x');xlabel('Time (sec)');legend('x','x_{PF}')
subplot(212)
plot(t,x_pf-x,'b-');hold on;grid on;
plot(t,sqrt(P_pf),'m-',t,-sqrt(P_pf),'m-');
ylabel('\delta x');xlabel('Time (sec)');


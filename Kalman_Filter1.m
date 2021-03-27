%=========================================================================%
%                                HW4_P3.m
%                           (AEM 5451/EE 5251)
%                   
%   This m-script simulates Problem 3 in HW4.
%
%   Kai Ye
%   October 17, 2019
%=========================================================================

%   (1)     Clear workspace

clear all;
close all;
clc;

%   (2)     Define constants

dt = 1;
t_end = 10;
t = [0:dt:t_end]';
drl = length(t);            %   Data record length

%   (3)     Generate true state, measurement history and system model

%   (i)     True states
r = 10^0.5;
x_2 = 250*ones(drl,1) + r*randn(drl,1);      %   True food supply
x_1 = 650*ones(drl,1);  %   Initialize population

% population true model
for iter = 2:drl
    x_1(iter) = x_1(iter-1)*0.5 + 2*x_2(iter-1);   
end

%   (ii)    Measurement

y = x_1 + r*randn(drl,1);     %   Measurement history with noise~N(0,10)

%   (iii)   System model

%   Assume state vector x = [phase error;frequncy error]

F = [0.5 2;0 1];
q = 10^0.5;
Q = r^2.*eye(2);  %%% Process noise
H = [1 0];
R = [r^2];
I = eye(2);


%   (3)     Begin Kalman Filter

x_hat = zeros(2,drl);           %   Place holder for state history
P = zeros(2,2,drl);             %   Place holder for covariance history
%K = zeros(2,drl);

%   (i)     Initialize state and covraince matrix using all available
%   information

x_hat(:,1) = [600;200];             %   Initial estimation
P(:,:,1) = [500 0;0 200];              %   Initial covariance matrix

%   (ii)    Begin KF loop

wb = waitbar(0,'Running KF simulation ...');

for k = 2:m
    
    waitbar(k/drl,wb);
    
    
    %   Forward filter
    %   Time update
    
    x_hat(:,k) = F*x_hat(:,k-1);
    P(:,:,k) = F*P(:,:,k-1)*F' + Q;
    
    %   Measurement update
    
    K(:,k) = P(:,:,k)*H'*inv(H*P(:,:,k)*H' + R);
    x_hat(:,k) = x_hat(:,k) + K(:,k)*(y(k) - H*x_hat(:,k));
    P(:,:,k) = (I - K(:,k)*H)*P(:,:,k);
    
    
end

close(wb)

%   Plot results

%   (i) State estimation history

figure              
subplot(221)
h1 = plot(t,x_hat(1,:),'r-');hold on;
h2 = plot(t,x_1,'b-');grid on;
h3 = ylabel('$\hat{p}$','Interpreter','latex');
h4 = title('Population');
h5 = legend('Estimated','True');
h6 = xlabel('Time (step)');

subplot(222)
h1 = plot(t,x_hat(2,:),'r-');hold on;
h2 = plot(t,x_2,'b-');grid on;
h3 = ylabel('$\hat{f}$','Interpreter','latex');
h4 = xlabel('Time (step)');
h5 = legend('Estimated','True');
h6 = title('Food supply');

%   (ii)    Estimation error history            
subplot(223)
h1 = plot(t,sqrt(squeeze(P(1,1,:))),'g-');hold on;
h2 = plot(t,sqrt(squeeze(P(2,2,:))),'r-');grid on;
h3 = ylabel('$\sigma$','Interpreter','latex');
h4 = title('Estimation Error History');
h5 = legend('Population','Food');
h6 = xlabel('Time (step)');


subplot(224)
h1 = plot(t(2:drl),K(1,2:drl),'r-');hold on;
h2 = plot(t(2:drl),K(2,2:drl),'g-');
%h3 = ylabel('Kalman Gain Matrix','Interpreter','latex');
h4 = xlabel('Time (step)');
h5 = legend('K1','K2');
h6 = title('Kalman Gain Matrix');

stdd = [std(x_1-x_hat(1,:)') std(x_2-x_hat(2,:)')];


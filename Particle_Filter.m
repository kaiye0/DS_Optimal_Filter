%=========================================================================%
%                                HW5_P5.m
%                           (AEM 5451/EE 5251)
%                   
%   This m-script simulates Problem 5 in HW5.
%
%   Kai Ye
%   November 29, 2019
%=========================================================================

%   (1)     Clear workspace
clearvars -except x y
% clear all;
% close all;
% clc;

%   (2)     Define constants

dt = 1;
t_end = 100;
t = [0:dt:t_end]';
drl = length(t);            %   Data record length
%Q = 0.1;
%Q = 1;
Q = 10;
R = 1;

%   (3)     Generate true state, measurement history and system model
% 
% %   (i)     True states
% 
% x = zeros(drl,1);%   Initialization
% x(1)= 0.1;
% for iter = 2:drl
%     x(iter) = 0.5*x(iter-1) + 25*x(iter-1)/(1 + x(iter-1)^2) + 8*cos(1.2*(iter-2)) + sqrt(Q)*randn(1,1);   
% end
% 
% %   (ii)    Measurement
% 
% y = 1/20*x.^2 + sqrt(R)*randn(drl,1);     %   Measurement history with noise~N(0,1)

%   (3)     Begin Filters, default variables for EKF

x_hat = zeros(1,drl);           %   Place holder for state history
x_hatukf = zeros(1,drl);
x_hatpf = zeros(1,drl);
P = zeros(1,drl);             %   Place holder for covariance history

%   (i)     Initialize state and covraince matrix using all available
%   information

x_hat(1) = 0.1;             %   Initial estimation
x_hatukf(1) = 0.1; 
x_hatpf(1) = 0.1; 
P(1) = 2;              %   Initial covariance matrix
P_ukf(1) = 2; 

%   (ii)   Paticle Filter initialization
N = 10;
%N = 100;
%N = 1000;
x_pf = ones(N,1)*x(1) + sqrt(P(1))*randn(N,1);

%   (iii)    Begin KF loop

wb = waitbar(0,'Running KF simulation ...');

for k = 2:drl
    
    waitbar(k/drl,wb);
    
    
    %   Extended Kalman Filter
    %   Time update
    
    x_hat(k) = 0.5*x_hat(k-1) + 25*x_hat(k-1)/(1+x_hat(k-1)^2) + 8*cos(1.2*(k-2));
    F = 0.5 + 25*(1- x_hat(k-1)^2)/(1+x_hat(k-1)^2)^2;
    P(k) = F*P(k-1)*F' + Q;
    
    %   Measurement update
    
    H = 0.1*x_hat(k);
    K(k) = P(k)*H'*(H*P(k)*H' + R)^(-1);
    x_hat(k) = x_hat(k) + K(k)*(y(k) - x_hat(k)^2/20);
    P(k) = (1 - K(k)*H)*P(k);
    
     
    %   Unscented Kalman Filter
    %   Time update
    
    x1 = x_hatukf(k-1) + sqrt(P_ukf(k-1));
    x2 = x_hatukf(k-1) - sqrt(P_ukf(k-1));
    x_hat1 =  0.5*x1 + 25*x1/(1+x1^2) + 8*cos(1.2*(k-2));
    x_hat2 =  0.5*x2 + 25*x2/(1+x2^2) + 8*cos(1.2*(k-2));
    x_hatukf(k) = (x_hat1 +x_hat2)/2;
    P_ukf(k) = ((x_hat1 - x_hatukf(k))^2 + (x_hat2 - x_hatukf(k))^2)/2 + Q;
    
    %   Measurement update
    x1 = x_hatukf(k) + sqrt(P_ukf(k));
    x2 = x_hatukf(k) - sqrt(P_ukf(k));
    y1 = 1/20*x1^2;
    y2 = 1/20*x2^2;
    y_hat = (y1 + y2)/2;
    P_y = ((y1 - y_hat)^2 + (y2- y_hat)^2)/2 + R;
    P_xy = ((y1 - y_hat)*(x1 - x_hatukf(k)) + (y2- y_hat)*(x2 - x_hatukf(k)))/2;
    K_ukf(k) = P_xy*P_y^(-1);
    x_hatukf(k) = x_hatukf(k) + K_ukf(k)*(y(k) - y_hat);
    P_ukf(k) = P_ukf(k) - K_ukf(k)*P_y*K_ukf(k)';

    %Particle filter
    x_pfpri = 0.5.*x_pf + 25.*x_pf./(1+x_pf.^2) + 8*cos(1.2*(k-2)) + sqrt(Q)*randn(N,1);
    y_pf = x_pfpri.^2/20;
    v = y(k)*ones(N,1) - y_pf;
    q = 1/sqrt(2*pi)/sqrt(R)*exp(v.^2/R/2);
    qsum = sum(q);
    q = q/qsum;
    
    %resampling
    qcumsum = cumsum([q(:).']);
    for i =1:N
        r = rand;   
        for j = 1:N
            if qcumsum(j)>= r
                x_pf(i) = x_pfpri(j);
                break;
            end
        end
    end
    x_hatpf(k) = mean(x_pf);
    
end
close(wb)

%Calculate average RMS state estimation error
%RMSE = sqrt((norm(x - x_hatpf'))^2/drl);
RMSE(1) = sqrt(mean((x-x_hat').^2));
RMSE(2) = sqrt(mean((x-x_hatukf').^2));
RMSE(3) = sqrt(mean((x-x_hatpf').^2));


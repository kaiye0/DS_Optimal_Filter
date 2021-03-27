%=========================================================================%
%                                HW4_P3.m
%                           (AEM 5451/EE 5251)
%                   
%   This m-script simulates Problem 3 in HW4.
%
%   Kai Ye
%   October 25, 2019
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
n = 1;
P_his = zeros(2,n);

for test = 1:n
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
Q = [1e-6 0; 0 10];  %%% Process noise
H = [1 0];
R = [r^2];


%   (3)     Begin Kalman Filter with smoother

%   (i)     Initialize state and covraince matrix after using smoother
x_hatsm = zeros(2,drl);           %   Place holder for state history
P_sm = zeros(2,2,drl);             %   Place holder for covariance history
K_sm = zeros(2,2,drl);           %   Place holder for state history

%   (ii)    Begin KF loop

%wb = waitbar(0,'Running KF simulation ...');

for m = 1:drl
    
%    waitbar(m/drl,wb);
    % Initialize state and covraince matrix every loop
    x_hat = zeros(2,drl);           %   Place holder for state history
    P = zeros(2,2,drl);             %   Place holder for covariance history
    K = zeros(2,drl);           %   Place holder for state history
    s = zeros(2,drl);           %   Place holder for state history
    I = zeros(2,2,drl);             %   Place holder for covariance history
    x_hat(:,1) = [600;200];             %   Initial estimation
    P(:,:,1) = [500 0;0 200];           %   Initial covariance matrix
    I(:,:,drl) = (1e-6)*eye(2);

if  m>1
for k = 2:m
    %Forward Filter
    x_hat(:,k) = F*x_hat(:,k-1);
    P(:,:,k) = F*P(:,:,k-1)*F' + Q;
    K(:,k) = P(:,:,k)*H'*inv(H*P(:,:,k)*H' + R);
    x_hat(:,k) = x_hat(:,k) + K(:,k)*(y(k) - H*x_hat(:,k));
    P(:,:,k) = (eye(2) - K(:,k)*H)*P(:,:,k);
end
end

if m<drl 
for j = drl:-1:m+1
    %Backward Filter
    I(:,:,j) = I(:,:,j) + H'*inv(R)*H;
    s(:,j) = s(:,j) + H'*inv(R)*y(j);
    %I(:,:,j-1) = F'*inv(inv(I(:,:,j))+Q)*F;
    I(:,:,j-1) = F'*(inv(Q)-inv(Q)*inv(I(:,:,j)+inv(Q))*inv(Q))*F;
    s(:,j-1) = I(:,:,j-1)*inv(F)*inv(I(:,:,j))*s(:,j);
end
end
    
    K_sm(:,:,m) = inv(I(:,:,m))*inv(P(:,:,m)+inv(I(:,:,m)));
    x_hatsm(:,m) = K_sm(:,:,m)*x_hat(:,m) + (eye(2)-K_sm(:,:,m))*inv(I(:,:,m))*s(:,m);
    P_sm(:,:,m) = inv(inv(P(:,:,m))+I(:,:,m));

end
    P_his(1,test) = x_1(1)-x_hatsm(1,1);
    P_his(2,test) = x_2(1)-x_hatsm(2,1);
    result = var(P_his');
%     P_his(1,test) = P_sm(1,1,1);
%     P_his(2,test) = P_sm(2,2,1);
%     result = mean(P_his,2);
end

%close(wb)

%   Plot results

figure
%      Estimation error history            
subplot(211)
h1 = plot(t,squeeze(P(1,1,:)),'g-');hold on;
h2 = plot(t,squeeze(P_sm(1,1,:)),'r-');grid on;
h3 = ylabel('$\sigma^2$','Interpreter','latex');
h4 = title('Estimation error variance of population');
h5 = legend('Posteriori of forward estimate','Smoothed state estimate');
h6 = xlabel('Time (step)');


subplot(212)
h1 = plot(t,squeeze(P(2,2,:)),'g-');hold on;
h2 = plot(t,squeeze(P_sm(2,2,:)),'r-');grid on;
h3 = ylabel('$\sigma^2$','Interpreter','latex');
h4 = title('Estimation error variance of food supply');
h5 = legend('Posteriori of forward estimate','Smoothed state estimate');
h6 = xlabel('Time (step)');




%=========================================================================%
%                                HW3_P3.m
%                           (AEM 5451/EE 5251)
%                   
%   This m-script simulates Problem 3 in HW3.
%
%   Kai Ye
%   October 17, 2019
%=========================================================================

%   (1)     Clear workspace

clear all;
close all;
clc;

%   (2)     Define constants
dt = 0.1;  %Time step at 10Hz
t_end = 2;
t = [0:dt:t_end]';
drl = length(t);            %   Data record length

%   (3)     Generate true state, measurement history and system model
% initial continuous model dx/dt = Ax +Bu
R = 3; L = 1; C = 0.5;
A = [0 1/C;-1/L -R/L];
B = [0 ; 1/L];
%   continuous to discrete
F = expm(A*dt);
G = (F - eye(2)) / A * B;

%   (i)     True states
u = randn(drl,1);      %   input voltage is zero-mean, unity variance white noise
x = zeros(2,drl);      %   Initialize x(t) = [Vc  IL]

for iter = 2:drl
    x(:,iter) = F*x(:,iter-1) + G*u(iter-1);  %   True state x_k = F*x_(k-1) + G*u_(k-1)
end        

%   (ii)    Measurement

r = 1;                     %   Capacitor measurement noise 
y = x(1,:) + r*randn(1,drl);     %   Measurement history with noise~N(0,1)

%   (iii)   System model

%   Assume state vector x = [capacitor voltage;inductor current]
q = 1;
Q = dt*[q^2 0;0  q^2];
H = [1 0];
R = [r^2];
I = eye(2);

%   (3)     Begin Kalman Filter

x_hat_pri = zeros(2,drl);           %   Place holder for state history
x_hat_post = zeros(2,drl);          
P_pri = zeros(2,2,drl);             %   Place holder for covariance history
P_post = zeros(2,2,drl); 

%   (i)     Initialize state and covraince matrix using all available
%   information

x_hat_pri(:,1) = [0;0];             %   Guessing as we have no knowledge
x_hat_post(:,1) = [0;0];             
P_pri(:,:,1) = 0*eye(2);        %   Initial P0 = 0
P_post(:,:,1) = 0*eye(2);

%   (ii)    Begin KF loop

wb = waitbar(0,'Running KF simulation ...');

for k = 2:drl
    
    waitbar(k/drl,wb);
    
    %   Time update
    
    x_hat_pri(:,k) = F*x_hat_post(:,k-1) + G*u(k-1);
    P_pri(:,:,k) = F*P_post(:,:,k-1)*F' + Q;
    
    %   Measurement update
    
    K = P_pri(:,:,k)*H'*inv(H*P_pri(:,:,k)*H' + R);
    x_hat_post(:,k) = x_hat_pri(:,k) + K*(y(k) - H*x_hat_pri(:,k));
    P_post(:,:,k) = (I - K*H)*P_pri(:,:,k);
    
    
end

close(wb)

var_Pri = zeros(drl,1);
var_Post = zeros(drl,1);
for counter = 1:drl
    var_Pri(counter) = var(x(2,1:counter)-x_hat_pri(2,1:counter));
    %var_Pri(counter) = var(x_hat_pri(2,1:counter));
    var_Post(counter) = var(x(2,1:counter)-x_hat_post(2,1:counter));
end  

%   Plot results
figure              
h1 = plot(t,var_Post,'b-');hold on;
h2 = plot(t,var_Pri,'r-');
h3 = plot(t,squeeze(P_pri(2,2,:)),'g-');grid on;
h4 = ylabel('Variance of $\hat{I}_L$ error','Interpreter','latex');
h5 = xlabel('Time (sec)');
h6 = legend('priori','posteriori');
h7 = title('Estimation Error History');


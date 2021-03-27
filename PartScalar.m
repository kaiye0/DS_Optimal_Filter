function [errRMSKalman, errRMSParticle] = PartScalar(N,Q)
if ~exist('N', 'var')
    N = 100;
end
if ~exist('Q', 'var')
    Q = 1;
end

x = 0.1;
R = 1;
tf = 50;

xhat = x;
P = 2;
xhatPart = x;

for i = 1:N
    xpart(i) = x+sqrt(P) * randn;
end

xArr = [x];
yArr = [x^2/20 + sqrt(R)*randn];
xhatArr = [xhat];
xhatPartArr = [xhatPart];

close all;

for k = 1:tf
    x = 0.5*x + 25*x/(1+x^2) + 8*cos(1.2*(k-1)) + sqrt(Q)*randn;
    y = x^2/20 + sqrt(R)*randn;
    F = 0.5 + 25*(1-xhat^2)/(1+xhat^2)^2;
    P = F*P*F' + Q;
    H = xhat/10;
    K = P * H' * (H * P * H' + R)^(-1);
    xhat = 0.5 * xhat + 25 * xhat / (1 + xhat ^ 2) + 8 * cos(1.2*(k-1));
    xhat = xhat + K * (y - xhat^2/20);
    P = (1 - K * H) * P;
    
    for i = 1 : N
        xpartminus(i) = 0.5 * xpart(i) + 25 * xpart(i) / (1 + xpart(i) ^ 2) + 8 *cos(1.2*(k-1)) + sqrt(Q) * randn;
        ypart = xpartminus(i)^2 / 20;
        vhat = y - ypart;
        q(i) = (1 / sqrt(R) / sqrt(2*pi)) * exp(-vhat^2 / 2 / R);
    end
    
    qsum = sum(q);
    if qsum < eps
        q = ones(size(q)) / N;
    else
        for i = 1 : N
            q(i) = q(i) / qsum;
        end
    end
    
    for i = 1:N
        u = rand;
        qtempsum = 0;
        for j = 1 : N
            qtempsum = qtempsum + q(j);
            if qtempsum >= u
                xpart(i) = xpartminus(j);
                break;
            end
        end
    end
    
    xhatPart = mean(xpart);
    
    xArr = [xArr x];
    yArr = [yArr y];
    xhatArr = [xhatArr xhat];
    xhatPartArr = [xhatPartArr xhatPart];
end

errRMSKalman = sqrt((norm(xArr - xhatArr))^2 / tf)
errRMSParticle = sqrt((norm(xArr - xhatPartArr))^2 / tf)

    
    
    
    
    
    
    
    
    
    
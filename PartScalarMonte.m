NMonte = 100;
% N = 100;
% Q = 0.1;

% N = 100;
% Q = 1;
% 
% N = 100;
% Q = 10;
% 
% N = 10;
% Q = 1;
% 
N = 1000;
Q = 1;

for i =1:NMonte
    [errKalman(i),errParticle(i)] = PartScalar(N,Q);
    fprintf('.');
end
fprintf('\n');
disp(['N = ', num2str(N),', Q = ',num2str(Q),', Ave RMS Est Err = ',num2str(mean(errKalman)), ' (Kalman), ', num2str(mean(errParticle)), '(Particle)'])
% disp(['N = ', num2str(N),', Ave RMS Est Err = ',num2str(mean(errKalman)), ' (Kalman), ', num2str(mean(errParticle)), '(Particle)'])

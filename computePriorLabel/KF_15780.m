%function [KF_xy]=KF(measurements, model)

%
% Load the T={0,..,100} data from each of the two observation towers
%
x = load('obs.dat')';
%
% Invariant definitions for the problem
%
% Time step is 1 second
delta_t = 1;
% Random acceleration, noise for w_t's Gaussian
Gamma = zeros(6,6);
Gamma(1,1) = 0.0001;
Gamma(2,2) = 0.0001;
Gamma(3,3) = 0.0001;
Gamma(4,4) = 0.0001;
Gamma(5,5) = 0.1;
Gamma(6,6) = 0.1;
% Position noise in observation data for v_t's Gaussian
Sigma = zeros(4,4);
Sigma(1,1) = 100;
Sigma(1,2) = 10;
Sigma(2,1) = 10;
Sigma(2,2) = 100;
Sigma(3,3) = 100;
Sigma(3,4) = 10;
Sigma(4,3) = 10;
Sigma(4,4) = 100;
% A is the transition matrix
A = zeros(6,6);
% Update distance based on d_t + v_t*t + 1/2*a_t*t^2
A(1,:) = [1 0 delta_t 0 (0.5 * (delta_t)^2) 0];
A(2,:) = [0 1 0 delta_t 0 (0.5 * (delta_t)^2)];
% Update velocity based on v_t + a_t*t
A(3,:) = [0 0 1 0 delta_t 0];
A(4,:) = [0 0 0 1 0 delta_t];
% Update acceleration based on a_t
A(5,:) = [0 0 0 0 1 0];
A(6,:) = [0 0 0 0 0 1];
% C is the observation matrix
C = zeros(4,6);
% Each observation is derived solely from the position of the UAV at a
% certain time t, plus some Gaussian noise that we'll add in later

C(1,:) = [1 0 0 0 0 0];
C(2,:) = [0 1 0 0 0 0];
C(3,:) = [1 0 0 0 0 0];
C(4,:) = [0 1 0 0 0 0];
%
% Kalman filter
% Start with an initial state (subject to some noise), then work our way
% through the next t_max timesteps using Kalman equations
%
% t_max is total number of time steps
t_max = 100;
z = zeros(6,t_max);
% Initial state: (0,0)m at velocity (6,6)m/s and acceleration (0,0)m/s^2
% subject to noise V_0 = 6x6 identity matrix
z(:,1) = [0; 0; 6; 6; 0; 0];
% Add on the z_0 = N(z_0 | u_0, V_0) Gaussian noise
z(:,1) = z(:,1) + mvnpdf(z(:,1), 0, eye(6));
% Error covariance matrix of our system
cov_V = eye(6);
% z_{t+1} = A z_t + w_t
% x_t = C z_t + v_t
for i=2:(t_max)
    % Kalman Filter equations from:
    % http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
    % Predict Step
    temp_z = A * z(:,i-1);
    temp_cov_V = A * cov_V * A' + Gamma;
    % Correct Step
    gain = temp_cov_V * C' / (C * temp_cov_V * C' + Sigma);
    z(:,i) = temp_z + gain*(x(:,i) - C * temp_z);
    cov_V = (eye(6) - gain*C) * temp_cov_V;
end
%
% Plot the sensor data and sequence of posterior distribution means for the
% (x, y)-positions of the UAV.
%
figure
plot(x(1,:),x(2,:),'r:',x(3,:),x(4,:),'r:',z(1,:),z(2,:),'b-');
xlabel('Horizontal Distance (m)');
ylabel('Vertical Distance (m)');
title('Sensor Data and Estimated Position of a Rogue UAV');
legend('Sensor #1', 'Sensor #2', 'Posterior');
%
% Report UAV's posterior position distribution at t = 100. This is a
% multivariate Gaussian, and can be reported as its mean and 6x6 covariance
% matrix.
%
z(:,t_max)
cov_V


%



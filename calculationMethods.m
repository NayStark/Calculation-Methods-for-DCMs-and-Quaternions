%% AERSP 450 - Homework #4
% Sunay Neelimathara - FA 24 - Professor Eapen

%% Problem 1
clc;
clear;
close all;

T = readtable('SensorData.csv'); 

wx = T.wx;
wy = T.wy;
wz = T.wz;

timeData = datetime(T.time, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS''Z''',...
'TimeZone', 'UTC');

timeDifferences =  timeData - timeData(1);

time_diff = seconds(diff(timeDifferences));
average_interval = mean(time_diff);
frequency = 1 / average_interval;
fprintf('Sensor Frequency: %.2f Hz\n\n', frequency);

t = seconds(timeDifferences);

yaw0 = 30;
pitch0 = 70;
roll0 = 20;
theta1 = yaw0;
theta2 = pitch0;
theta3 = roll0;

% Find Euler angle sequence
% Yaw-Pitch-Roll is 3-2-1 Euler Angle Sequence


% Yaw-Pitch-Roll into DCM
% From our Euler Angle Relations:
C_ypr = [cosd(theta1)*cosd(theta2)  cosd(theta2)*sind(theta1)  -sind(theta2)
         sind(theta3)*sind(theta2)*cosd(theta1)-cosd(theta3)*sind(theta1)  sind(theta3)*sind(theta2)*sind(theta1)+cosd(theta3)*cosd(theta1)  sind(theta3)*cosd(theta2)
         cosd(theta3)*sind(theta2)*cosd(theta1)+sind(theta3)*sind(theta1)  cosd(theta3)*sind(theta2)*sind(theta1)-sind(theta3)*cosd(theta1)  cosd(theta3)*cosd(theta2)];

% Printing the calculated DCM
fprintf('DISPLAYING DCM: \n');
fprintf('\n The DCM for the yaw-pitch-roll angles (3-2-1 euler angle sequence) is as follows: \n \n');
disp(C_ypr);
fprintf('\n');

% Verifying that DCM is valid, by seeing if its determinant is equal to 1.
fprintf('DCM VERIFICATION: \n \n');
verification_DCM = det(C_ypr);
fprintf('||DCM|| = ');
disp(verification_DCM);
fprintf('DCM is valid. \n \n \n');

% Yaw-Pitch-Roll into quaternion
% Using the equation from the roadmap, we can convert to quaternions:
beta0 = (1/2) * sqrt(C_ypr(1, 1) + C_ypr(2, 2) + C_ypr(3, 3) + 1);
beta1 = (C_ypr(2, 3) - C_ypr(3, 2)) / (4*beta0);
beta2 = (C_ypr(3, 1) - C_ypr(1, 3)) / (4*beta0);
beta3 = (C_ypr(1, 2) - C_ypr(2, 1)) / (4*beta0);

quat1 = [beta0, beta1, beta2, beta3]';

% Displaying the quaternions:
fprintf('DISPLAYING QUATERNIONS: \n');
fprintf('\n The quaternions for the yaw-pitch-roll angles (3-2-1 euler angle sequence) are as follows: \n \n \t Beta0 = ');
disp(beta0);
fprintf('\t Beta1 = ');
disp(beta1);
fprintf('\t Beta2 = ');
disp(beta2);
fprintf('\t Beta3 = ');
disp(beta3);
fprintf('\n');


% Plot angular velocities as a function of time
% We are given the omegas, and t.
figure;
hold on;
title('\omega vs. time');
plot(t, wx);
plot(t, wy);
plot(t, wz);
hold off;
xlabel('time (seconds)');
ylabel('\omega (m/s)');
legend('\omega_{x}', '\omega_{y}', '\omega_{z}');
grid on;

% Comment on the nature of the data: (a) what is the frequency of the sensor output. (b) are the sensor outputs at equal time intervals?
total_time = t(1544, :);
time_interval = size(t);
frequency = 1/total_time;

fprintf('\n TIME DIFFERENCES: \n');
%disp(t);

%% Problem 2

% Part 1
% Cdot_BN = -[w_tilde]*C_BN

omega = [wx, wy, wz];

w_tilde = @(w) [0 -w(3) w(2);
                w(3) 0 -w(1);
                -w(2) w(1) 0];

DCM_prop = @(time, C_BN) reshape(-w_tilde(interp1(t, omega, time, 'linear')') * reshape(C_BN, 3, 3), [], 1);

C0 = C_ypr(:);

[t_out, C_ypr_ode] = ode45(@(time, C_BN) DCM_prop(time, C_BN), t, C0);

C_ypr_matrices = zeros(3, 3, length(t_out));
DCM_Valid = zeros(length(t_out));

for i = 1:length(t_out)
    C_ypr_matrices(:, :, i) = reshape(C_ypr_ode(i, :), 3, 3);
end

C_ypr_matrix = C_ypr_matrices(:, :, end);
DCM_Valid = det(C_ypr_matrix);


% Part 2
% C_BN(tk+1) = expm(-[w_tilde]*deltatk) C_BN(tk)

C_ypr_matrices_3 = zeros(3, 3, length(t));
C_ypr_matrices_3(:, :, 1) = C_ypr;

for k = 1:length(t)-1

    dt = t(k+1) - t(k);
    omega_k = omega(k, :);
    w_skew = w_tilde(omega_k);
    exp_w = expm(-w_skew * dt);
    C_ypr_matrices_3(:, :, k+1) = exp_w * C_ypr_matrices_3(:, :, k);

end

disp('Final Numerical DCM:')
disp(C_ypr_matrix);


disp('Final Analytical DCM:');
disp(C_ypr_matrices_3(:, :, end));
C_ypr_matrix_3_f = C_ypr_matrices_3(:, :, end);

% Part 3

identity_3 = eye(3);
error_tk = zeros(3, 3, length(t));
for k = 1:length(t)
    C_ypr_num = C_ypr_matrices(:, :, k);
    C_ypr_analytic = C_ypr_matrices_3(:, :, k);
    error_tk(:, :, k) = C_ypr_matrices(:, :, k) * C_ypr_matrices_3(:, :, k)' - identity_3;
end

error_norm = zeros(length(t), 1);
for k = 1:length(t)
    error_norm(k) = norm(error_tk(:, :, k));
end

figure;
plot(t, error_norm);
xlabel('Time (s)');
ylabel('Error Norm');
title('Time History of DCM Error');
grid on;

% Part 4
yaw = zeros(1, length(t));
pitch = zeros(1, length(t));
roll = zeros(1, length(t));

for k = 1:length(t)
    C_ypr = C_ypr_matrices(:, :, k);
    roll(k) = atan2d(C_ypr(2, 3), C_ypr(3, 3));
    pitch(k) = atan2d(-C_ypr(1, 3), sqrt(C_ypr(1, 1)^2 + C_ypr(1, 2)^2));
    yaw(k) = atan2d(C_ypr(1, 2), C_ypr(1, 1));
end

figure;
subplot(2, 1, 1);
hold on;
plot(t, yaw);
plot(t, pitch);
plot(t, roll);
xlabel('Time (s)');
ylabel('Yaw (rad)');
title('Yaw-Pitch-Roll vs Time - Numerical');
legend('yaw', 'pitch', 'roll')
grid on;
hold off;


for k = 1:length(t)
    C_ypr = C_ypr_matrices_3(:, :, k);
    roll(k) = atan2d(C_ypr(2, 3), C_ypr(3, 3));
    pitch(k) = atan2d(-C_ypr(1, 3), sqrt(C_ypr(1, 1)^2 + C_ypr(1, 2)^2));
    yaw(k) = atan2d(C_ypr(1, 2), C_ypr(1, 1));
end

subplot(2, 1, 2);
hold on;
plot(t, yaw);
plot(t, pitch);
plot(t, roll);
xlabel('Time (s)');
ylabel('Yaw (rad)');
title('Yaw-Pitch-Roll vs Time - Analytical');
legend('yaw', 'pitch', 'roll')
grid on;
hold off;

%% Problem 3
% Part 1

quat0 = [beta0, beta1, beta2, beta3];
beta_calc = zeros(length(t), length(quat0));
beta_calc(1, :) = quat0;

B_omega = @(w) [0 -w(1) -w(2) -w(3);
              w(1) 0 w(3) -w(2);
              w(2) -w(3) 0 w(1);
              w(3) w(2) -w(1) 0];

omega_transpose = omega';

for k = 1:length(t)-1

    dt = t(k+1) - t(k);

    w = omega_transpose(:, k);
    B_w = B_omega(w);
    phi_quat = expm(0.5 * B_w * dt);
    beta_calc(k+1, :) = (phi_quat * beta_calc(k, :)')';

end

% Part 2
quat_const = ones(length(t));
figure;
subplot(2, 1, 1);
hold on;
plot(t, beta_calc(:, 1));
plot(t, beta_calc(:, 2));
plot(t, beta_calc(:, 3));
plot(t, beta_calc(:, 4));
plot(t, quat_const, 'g');
legend('\beta_{0}', '\beta_{1}', '\beta_{2}', '\beta_{3}', 'Quaternion Constant: \beta = 1');
xlabel('time (seconds)');
title('Quaternion Propagation from Equations of Motion');
hold off;
grid on;

% Part 3
% yaw-pitch-roll from previous section analytical to quaternions
beta0_f = zeros(1, length(t));
beta1_f = zeros(1, length(t));
beta2_f = zeros(1, length(t));
beta3_f = zeros(1, length(t));
% test = zeros(4, length(t));

quaternion_array = zeros(4, length(t));
quaternion_array(:, 1) = quat1;
norm_quat = zeros(1, length(t));
norm_quat(1) = norm(quaternion_array);

for k= 1:length(t)-1
    C_ypr = C_ypr_matrices_3(:, :, k);
    beta0_f = (1/2) * sqrt(C_ypr(1, 1) + C_ypr(2, 2) + C_ypr(3, 3) + 1);
    beta1_f = (C_ypr(2, 3) - C_ypr(3, 2)) / (4*beta0_f);
    beta2_f = (C_ypr(3, 1) - C_ypr(1, 3)) / (4*beta0_f);
    beta3_f = (C_ypr(1, 2) - C_ypr(2, 1)) / (4*beta0_f);
    quaternion_array(:, k+1) = [beta0_f; beta1_f; beta2_f; beta3_f];
    if dot(quaternion_array(:, k), quaternion_array(:, k+1)) < 0
        quaternion_array(:, k+1) = -quaternion_array(:, k+1);
    end
    quaternion_array(:, k+1) = quaternion_array(:, k+1) / norm(quaternion_array(:, k+1));
    norm_quat(k+1) = norm(quaternion_array(:, k+1));
end

subplot(2, 1, 2);
hold on;
plot(t, quaternion_array(1, :));
plot(t, quaternion_array(2, :));
plot(t, quaternion_array(3, :));
plot(t, quaternion_array(4, :));
plot(t, norm_quat, 'g');
legend('\beta_{0}', '\beta_{1}', '\beta_{2}', '\beta_{3}', 'Quaternion Constant: \beta = 1');
xlabel('time (seconds)');
title('Quaternions from Analytical Propagation of Yaw-Pitch-Roll Angles');
hold off;
grid on;

% Part 4
% Error
error_quaternion = zeros(4, length(t));
error_norm_quat = zeros(1, length(t));

for k = 1:length(t)
    % Analytical quaternion
    beta1f = beta_calc(k, :); % Row vector
    beta2f = quaternion_array(:, k)'; % Convert column to row vector
    % Quaternion inverse (conjugate)
    beta2_inv = [beta2f(1), -beta2f(2:4)];
    % Quaternion multiplication
    delta_beta = quatmultiply(beta1f, beta2_inv); % Output is a row vector
    delta_beta = delta_beta'; % Convert to column vector
    % Quaternion error
    error_quaternion(:, k) = delta_beta - [1; 0; 0; 0]; % Ensure column alignment
end

% Compute norm of the quaternion error
for k = 1:length(t)
    error_norm_quat(k) = norm(error_quaternion(:, k));
end

% Plot the error norm
figure;
plot(t, error_norm_quat);
xlabel('Time (seconds)');
ylabel('Error Norm');
title('Time History of Quaternion Propagation Error');
grid on;

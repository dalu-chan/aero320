clc
clear
close all
J = [17, -3, 2; -3, 20, -4; 2, -4, 15];
angular_velocity = [.01;-.1;.05];
tspan = [0 20];
initial_state = [0; 0; 0; 1; 0; 0; 0];

% Create a function to represent the differential equations
ode_fun = @(t, y) omega2eulerquat(t, y, angular_velocity);

% Solve the ODE using ode45 with the corrected initial condition
[t, y] = ode45(ode_fun, tspan, initial_state);
% Extract Euler Angles and Quaternion components
euler_angles = y(:, 1:3); % Extract Euler Angles
quaternions = y(:, 4:7); % Extract Quaternion components

% Plot Euler Angles versus time
figure;
subplot(3, 1, 1);
plot(t, euler_angles(:, 1), 'r', t, euler_angles(:, 2), 'g', t, euler_angles(:, 3), 'b');
xlabel('Time (s)');
ylabel('Euler Angles (radians)');
legend('Roll', 'Pitch', 'Yaw');
title('Euler Angles vs. Time');

% Plot Quaternion components versus time
subplot(3, 1, 2);
plot(t, quaternions(:, 1), 'r', t, quaternions(:, 2), 'g', t, quaternions(:, 3), 'b', t, quaternions(:, 4), 'k');
xlabel('Time (s)');
ylabel('Quaternion Components');
legend('q1', 'q2', 'q3', 'q4');
title('Quaternion Components vs. Time');

% Plot angular velocity versus time
subplot(3, 1, 3);
plot(t, ones(61,1)*angular_velocity(1), t, ones(61,1)*angular_velocity(2), t, ones(61,1)*angular_velocity(3))
xlabel('Time (s)');
ylabel('Angular Velocity Components');
legend('x', 'y', 'z');
title('Angular Velocity vs. Time');
function dy = omega2eulerquat(t, y, angular_velocity)
%OMEGA2EULERQUAT This function uses initial state matrix y: roll, pitch, yaw, q1,
%q2, q3, q4, and then angular velocity
% Extract Euler Angles and Quaternion components from the state vector
roll = y(1);
pitch = y(2);
yaw = y(3);
q1 = y(4);
q2 = y(5);
q3 = y(6);
q4 = y(7);

% Calculate the rotation matrix from Euler Angles (ZYX order)
clcR = euler2rotmat(roll, pitch, yaw);

% Calculate the time derivatives of Euler Angles (angular rates)
roll_rate = angular_velocity(1);
pitch_rate = angular_velocity(2);
yaw_rate = angular_velocity(3);

% Calculate the time derivatives of the Quaternion components
q1_dot = 0.5 * (q4 * roll_rate - q3 * pitch_rate + q2 * yaw_rate);
q2_dot = 0.5 * (q3 * roll_rate + q4 * pitch_rate - q1 * yaw_rate);
q3_dot = 0.5 * (-q2 * roll_rate + q1 * pitch_rate + q4 * yaw_rate);
q4_dot = 0.5 * (-q1 * roll_rate - q2 * pitch_rate - q3 * yaw_rate);

% Calculate the time derivatives of Euler Angles
roll_dot = roll_rate + sin(roll) * tan(pitch) * pitch_rate + cos(roll) * tan(pitch) * yaw_rate;
pitch_dot = cos(roll) * pitch_rate - sin(roll) * yaw_rate;
yaw_dot = sin(roll) / cos(pitch) * pitch_rate + cos(roll) / cos(pitch) * yaw_rate;

% Store the derivatives in a column vector
dy = [roll_dot; pitch_dot; yaw_dot; q1_dot; q2_dot; q3_dot; q4_dot];
end
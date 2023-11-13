clc
clear
close all

J = [17, -3, 2; -3, 20, -4; 2, -4, 15];
angular_velocity = [.01;-.1;.05];
initial_state = [0; 0; 0; 1; 0; 0; 0];

%% A
[V, D] = eig(J);
pmoment = diag(D);
disp("The principal moments are")
disp(pmoment)

%% B
% Assuming a mass of one
xL = sqrt(6*pmoment(1));
yL = sqrt(6*pmoment(2));
zL = sqrt(6*pmoment(3));
disp("The equivalent cuboid would have dimensions of length "+xL +" width "+yL+" and height "+zL)

%% C
tspan = [0 20];
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
legend('eta', 'q2', 'q3', 'q4');
title('Quaternion Components vs. Time');
% Plot angular velocity versus time
subplot(3, 1, 3);
plot(t, ones(61,1)*angular_velocity(1), t, ones(61,1)*angular_velocity(2), t, ones(61,1)*angular_velocity(3))
xlabel('Time (s)');
ylabel('Angular Velocity Components');
legend('x', 'y', 'z');
title('Angular Velocity vs. Time');

%% D
% Initialize the initial z-axis vector
initial_z_axis = [0; 0; 1];

% Initialize arrays to store z-axis motion in both frames
z_axis_inertial_frame = zeros(size(euler_angles));
z_axis_inertial_frame = zeros(size(euler_angles));

% Calculate z-axis motion for each time step
for i = 1:61
    % Convert Euler angles or quaternions to rotation matrices
    % You can use functions like euler2dcm or quat2dcm for this
    rotation_matrix = euler2rotmat(euler_angles(i,1),euler_angles(i,2),euler_angles(i,3));
    % Apply the rotation matrix to the initial z-axis
    z_axis_inertial_frame(i, :) = rotation_matrix * initial_z_axis;  % Replace 'rotation_matrix' with your actual rotation matrix
    z_axis_body_frame(i, :) = initial_z_axis;  % In the inertial frame, z-axis remains constant
end

% Plot the motion of the principal z-axis
figure;
subplot(2, 1, 1);
plot(t, z_axis_body_frame);
title('Principal Z-Axis Motion in Body Frame');
xlabel('Time (s)');
ylabel('Z-Axis Component');
legend('Zx', 'Zy', 'Zz');

subplot(2, 1, 2);
plot(t, z_axis_inertial_frame);
title('Principal Z-Axis Motion in Inertial Frame');
xlabel('Time (s)');
ylabel('Z-Axis Component');
legend('Zx', 'Zy', 'Zz');

%% E

tau = [1;-1;0];
initial = [0;0;0;0.01;-0.1;0.05;0;0;0;1;1;0;0;0;1;0;0;0;1];
wb0 = [0.01;-0.1;0.05]; % rad/s

ops = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t1,y1] = ode45(@eulerfun, tspan, initial, ops, D, tau);

figure
subplot(3,1,1)
plot(t1,y1(:,4:6))
grid on
title("Angular velocity vs. Time")
xlabel("Time(s)")
ylabel("w(rad/s)")
legend("wx","wy","wz","Location","best")

subplot(3,1,2)
plot(t1,y1(:,1:3))
grid on
title("Euler Angles vs. Time")
xlabel("Time(s)")
ylabel("Euler angles (rad)")
legend("phi","theta","psi","Location","best")

subplot(3,1,3)
plot(t1,y1(:,7:10))
grid on
title("Quaternions vs. Time")
xlabel("Time(s)")
ylabel("Quaternion")
legend("eps1","eps2","eps3","eta","Location","best")
hold off
disp(" ")
%% Problem 2

v = [20; 105; -10];
omega = [0.5; -0.1; 0.1];
m = 17474;
J = [2.44, 0, -1.2; 0, 27.0, 0; -1.2, 0, 30.0] * 1e6;

Ttrans = 0.5*m*norm(v)^2;

Trot = 0.5*omega'*J*omega;

Ttot = Ttrans + Trot;

disp("Total kinetic energy is "+Ttot+" J")

%% Functions Used
function dydt = eulerfun(t,y,J,Tc)
dydt = zeros(10,1);
phi = y(1);
theta = y(2);
psi = y(3);
w = y(4:6);
eps = y(7:9);
eta = y(10);
CbI = y(11:19); % Part d 

%phi, theta, and psi dot
dydt(1:3,1) = (1/cos(theta))*[cos(theta) sin(phi)*sin(theta) cos(phi)*sin(theta); 0 cos(phi)*cos(theta) -sin(phi)*cos(theta); 0 sin(phi) cos(phi)]*w;
%w dot
dydt(4:6) = J\(Tc-cross(w,J*w));
%eps
dydt(7:9) = (1/2)*(eta*eye(3)+crossm(eps))*w;
%eta
dydt(10) = -(1/2)*eps'*w;
%colums into 3X3
dydt(11:19) = -1*crossm(w)*reshape(CbI,[3,3]);
end

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

function R = euler2rotmat(roll, pitch, yaw)
%EULER2ROTMAT This function returns the rotation matrix for an inputed
%roll, pitch, and yaw
Rz = [cos(yaw), -sin(yaw), 0; sin(yaw), cos(yaw), 0; 0, 0, 1];
Ry = [cos(pitch), 0, sin(pitch); 0, 1, 0; -sin(pitch), 0, cos(pitch)];
Rx = [1, 0, 0; 0, cos(roll), -sin(roll); 0, sin(roll), cos(roll)];

R = Rz * Ry * Rx;
end
function crossProductMatrix = crossm(vector)
    % Check if the input vector is of the correct size
    if numel(vector) ~= 3
        error('Input vector must have exactly 3 elements');
    end

    % Extract components of the vector
    x = vector(1);
    y = vector(2);
    z = vector(3);

    % Create the cross product matrix
    crossProductMatrix = [0, -z, y;
                          z, 0, -x;
                         -y, x, 0];
end
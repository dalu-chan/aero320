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
%% Template for AuE 881 HW5
% This solution to the HW assignment is organized by task.  In each task,
% we follow the same template: set up the problem, solve it, plot/interpret
% the results, and answer the interpretation questions if there are any.

%% Initialization
% In general, the following two lines will make your code slower, but I
% suggest leaving them in for final testing so that you do not accidentally
% rely on stale result still existing in your workspace
clear all
close all

%% Base workspace initialization
% SI units have been used all through out the simulation unless otherwise
% stated

% 4>Point-mass road load model
m = 2000;    % Mass of the vehicle [kg]
r = 0.335;   % Radius of the tires [m]
C_d = 0.3;   % Aerodynamic drag coefficient [unitless]
A_f = 2.35;  % Frontal area [m^2]
C_r =0.009;  % Roll coefficienct [unitless] 
rho = 1.3;   % Density of air [kg/m^3]

% 1>Driver model
% A simple PI controller is used to implement reference speed tracking control
% PI controller parameters [unitless]
Kp = 3;
Ki = 1;
% Accelerator and brake pedal limits [unitless]
AccPedMax = 1;
AccPedMin = 0;
BrkPedMax = 0;
BrkPedMin = -1;

% 2>Engine-torqueConverter (ETc) model
eng_idle = 700*pi/30;    % Engine idle angular velocity [rad/s]
eng_redline = 5200*pi/30; % Engine redline angular velocity [rad/s]
% Vector definitions for 1-D map of maximum torque output from the ETc model  
% Map vector for maximum torque output from the ETc system [Nm]
% {It is the torque as measured at the torque converter turbine}
ETc_max_T = [ 500, 450, 350, 250, 200, 235.6, 257.2, 267.9, 267.3, 276.9,...
              277, 277.9, 283.2, 281.2, 281.2, 276.1, 270.9, 264.7, 255.4 ]; 
% Map vector for torque converter turbine angular speed [rad/s] 
ETc_w = [ 0, 200, 400, 600, 700, 999.1, 1397.6, 1800.5, 2199.5, 2600.6, 3002.9,...
          3402.5, 3802.4, 4002.4, 4203.4, 4602.8, 4802.5, 5003.5, 5203.8 ] *pi/30;   
% Vector definitions for 2-D map of fuel mass flow rate into the the engine(mfDot)
fmap_T = 0:10:290;  % Map vector for engine torque [Nm]
fmap_w = [ 700:100:5200 ]*pi/30;  % Map vector for engine angular speed [rad/s]
Q_lhv = 44e6; % Lower heating value of the fuel used (conventional gasoline) [J/kg] 
rho_f = 0.75; % Density of conventional gasoline [kg/l] 

% 3>Driveline model
gr = [3.2, 2.1, 1.5, 1, 0.85, 0.7];  % Gear ratio array of the automatic transmission [unitless] 
dr = 3.5;  % Gear ratio of the differential [unitless] 
fr = dr*gr;  % Final drive ratio array [unitless] 
g = 9.81;  % Acceleration due to gravity [m/s^2]
mu = 0.9;  % Traction coefficient [unitless]  {On dry asphalt flat road}
% Engine speed at shift points (assumed constant for simplification) [rad/s]
up_eng_speed = 2500*pi/30;
down_eng_speed = 1500*pi/30;
% Vehicle speed at shift points [m/s]
lim_up = r*up_eng_speed./fr;
lim_up(end) = 99999; 
lim_down = r*down_eng_speed./fr;
lim_down(1) = 0;
eff_dl =0.9;  % Driveline efficiency [unitless] 
BrkGain = 1e4;  % Brake effort gain [unitless]
 
%% Set simulation parameters
sim_end = 60;
time_step = 0.1;

%% Task 2
HW5 = sim('EmranAtef_Task2.slx',100);
figure(1)
plot(HW5.velocity);
title('2. Vehicle Velocity vs. Time [constant force of 10kN]')
xlabel('Time (s)')
ylabel('Vehicle Velocity (m/s)')
grid on;


%% Task 3 (500 Nm)
HW5 = sim('EmranAtef_Task3.slx',60);

figure(2)
plot(HW5.velocity);
title('3. The ACTUAL SPEED of the vehicle [Constant 500 Nm of the turbine]')
xlabel('Time (s)')
ylabel('Actual speed of the vehicle (m/s)')
grid on;

figure(3)
plot(HW5.F_wheels);
title('3. The FORCE at the wheels [Constant 500 Nm of the turbine]')
xlabel('Time (s)')
ylabel('The FORCE at the wheels (N)')
grid on;

figure(4)
plot(HW5.ETc_angVel);
title('3. The torque converter TURBINE VELOCITY [Constant 500 Nm of the turbine]')
xlabel('Time (s)')
ylabel('Torque Converter Turbine Velocity (rad/s)')
grid on;

figure(5)
plot(HW5.velocity.Data,HW5.F_wheels.Data);
title('3. Wheels Forces (N) vs. Vehicle Speed (m/s) [Constant 500 Nm of the turbine]')
xlabel('Actual speed of the vehicle (m/s)')
ylabel('The FORCE at the wheels (N)')
grid on;

%% Task 4

HW5 = sim('EmranAtef_Task4',200);

figure(6)
plot(HW5.velocity);
title('4. The ACTUAL SPEED of the vehicle')
xlabel('Time (s)')
ylabel('Actual speed of the vehicle (m/s)')
grid on;

figure(7)
plot(HW5.F_wheels);
title('4. The FORCE at the wheels')
xlabel('Time (s)')
ylabel('The FORCE at the wheels (N)')
grid on;

figure(8)
plot(HW5.ETc_angVel);
title('4. The torque converter TURBINE VELOCITY')
xlabel('Time (s)')
ylabel('Torque Converter Turbine Velocity (rad/s)')
grid on;

figure(9)
plot(HW5.velocity.Data,HW5.F_wheels.Data);
title('4. Wheels Forces (N) vs. Vehicle Speed (m/s)')
xlabel('Actual speed of the vehicle (m/s)')
ylabel('The FORCE at the wheels (N)')
grid on;

figure(10)
plot(HW5.ETc_Tq);
title('4. Turbine Torque')
xlabel('Time (s)')
ylabel('Torque Converter Turbine Torque (N.m)')
grid on;


%% Task 5

HW5 = sim('HW5_BaseModel',200);

figure(11)
plot(HW5.velocity);
title('5. The ACTUAL SPEED of the vehicle')
xlabel('Time (s)')
ylabel('Actual speed of the vehicle (m/s)')
grid on;

figure(12)
plot(HW5.F_wheels);
title('5. The FORCE at the wheels')
xlabel('Time (s)')
ylabel('The FORCE at the wheels (N)')
grid on;

figure(13)
plot(HW5.ETc_angVel);
title('5. The torque converter TURBINE VELOCITY')
xlabel('Time (s)')
ylabel('Torque Converter Turbine Velocity (rad/s)')
grid on;

figure(14)
plot(HW5.ETc_Tq);
title('5. Turbine Torque')
xlabel('Time (s)')
ylabel('Torque Converter Turbine Torque (N.m)')
grid on;

figure(15)
plot(HW5.RefVel);
hold on;
plot(HW5.velocity);
title('Ref vs. Act Velocities')
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on;
legend('RefVel','ActVel')

figure(16)
plot(HW5.error);
title('error')
xlabel('error')
ylabel('time (s)')

%% AuE881.HW4.Emran.Atef
% This solution to the HW assignment is organized by task.  In each task,
% we follow the same template: set up the problem, solve it, plot/interpret
% the results, and answer the interpretation questions if there are any.

%% Initialization
% In general, the following two lines will make your code slower, but I
% suggest leaving them in for final testing so that you do not accidentally
% rely on stale result still existing in your workspace
clear all;
close all;
clc
%% Simulink parameter definitions
% Define the vehicle parameters used in the Simulink model

mass = 1734;            % total vehicle mass [kg]

WD = 0.5;               % weight distribution [%] -- proportion on front axle
L = 2.7;               % wheelbase [m]
Jz = 2380;              % Moment of inertia of the vehicle [kgm^2]
Caf = 1450;             % Tire cornering stiffness front [N/deg]
Car = 3000;             % Tire cornering stiffness rear [N/deg]
Cmzf = 30;              % Tire aligning stiffness front [Nm/deg]
Cmzr = 20;              % Tire aligning stiffness rear [Nm/deg]

V = 40/3.6;             % velocity [m/s] (40 kmph is 40/3.6 m/s)

% some derived quantities
a = (1-WD)*L;           % front axle distance from cg [m]
b = WD*L ;              % rear axle distance from cg [m]

% Define constants used in the Simulink model
g = 9.81;               % gravity [m/s2]

% A few variable initializations necessary to make the initial Simulink
% template work -- Note: these values will need to be modified later for
% each simulation scenaria
steeringScenario = 1;  % start by running the first scenario (constant input)
steerAngle = 5;        % steer angle in [deg] for step input
enablePlotting = false; % turn on the plotting features (off is faster!)
maxYawRate = 20;       % [rad/s] -- stop the simulation if yaw rate is >

% %% Task 1
% 
% %calling the Simulink model, and displaying it
% figure('Name','Vehicle Position');
% HW4out1 = sim('HW4_Task2.slx',4);  
% disp(HW4out1.steer);
% 
% % ploting the Steering angle 
% figure(2)
% plot(HW4out1.steer);
% title('Steering Angle with Time')
% xlabel('Time (s)')
% ylabel('Steering Angle (degrees)')
% 
% % ploting the Yaw Rate
% figure(3)
% plot(HW4out1.yawrate);
% title('Yaw Rate (degrees/s) with Time')
% xlabel('Time (s)')
% ylabel('Yaw Rate (degrees/s)')
% 
% % ploting the Vehicle Position
% figure(4)
% plot(HW4out1.X.Data,HW4out1.Y.Data);
% title('Vehicle Position Y vs X')
% xlabel('X')
% ylabel('Y')
% % 
% %% Task 3
% 

R=100;
steeringScenario = 1;          % constant steering
steerAngle = 6;   % Ackermann steering angle in [deg] 
enablePlotting = false;        % turn on the plotting features (off is faster!)
maxYawRate = 20;               % [rad/s] -- stop the simulation if yaw rate is >

V = 111/3.6;                     % velocity [m/s] (1 kmph is 40/3.6 m/s)

%calling the Simulink model, and displaying it
HW4out1 = sim('HW4_Task2.slx',2000);  

% ploting the Vehicle Position
figure(15)
plot(HW4out1.X.Data,HW4out1.Y.Data);
grid on
title('Task 4.1: Vehicle Position - Ackermann steering angle and 1 km/h with 1000s run time')
xlabel('X')
ylabel('Y')

Kus=((WD*mass)/(2*Caf))-(((1-WD)*mass)/(2*Car))


% 
% % Scenario 1: A constant steering input – steeringScenario=1, steerAngle=5
% steeringScenario = 1;   
% steerAngle = 5;         % steer angle in [deg] for step input
% enablePlotting = false; % turn on the plotting features (off is faster!)
% maxYawRate = 20;        % [rad/s] -- stop the simulation if yaw rate is >
% 
% %calling the Simulink model, and displaying it
% HW4out1 = sim('HW4_Task2.slx',4);  
% 
% % ploting the Steering angle 
% figure(5)
% plot(HW4out1.steer);
% title('Scenario 1: Steering Angle with Time')
% xlabel('Time (s)')
% ylabel('Steering Angle (degrees)')
% 
% % ploting the Yaw Rate
% figure(6)
% plot(HW4out1.yawrate);
% title('Scenario 1: Yaw Rate (degrees/s) with Time')
% xlabel('Time (s)')
% ylabel('Yaw Rate (degrees/s)')
% 
% % ploting the Vehicle Position
% figure(7)
% plot(HW4out1.X.Data,HW4out1.Y.Data);
% title('Scenario 1: Vehicle Position Y vs X')
% xlabel('X')
% ylabel('Y')



% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% % Scenario 2: A left-right step steering input – steeringScenario=2, steerAngle=15
% steeringScenario = 2;   
% steerAngle = 15;        % steer angle in [deg] for step input
% enablePlotting = false; % turn on the plotting features (off is faster!)
% maxYawRate = 20;        % [rad/s] -- stop the simulation if yaw rate is >
% 
% %calling the Simulink model, and displaying it
% HW4out1 = sim('HW4_Task2.slx',4);  
% 
% % ploting the Steering angle 
% figure(8)
% plot(HW4out1.steer);
% grid on
% title('Scenario 2: Steering Angle with Time')
% xlabel('Time (s)')
% ylabel('Steering Angle (degrees)')
% 
% % ploting the Yaw Rate
% figure(9)
% plot(HW4out1.yawrate);
% grid on
% title('Scenario 2: Yaw Rate (degrees/s) with Time')
% xlabel('Time (s)')
% ylabel('Yaw Rate (degrees/s)')
% 
% % ploting the Vehicle Position
% figure(10)
% plot(HW4out1.X.Data,HW4out1.Y.Data);
% grid on
% title('Scenario 2: Vehicle Position Y vs X')
% xlabel('X')
% ylabel('Y')
%  
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%  
% % Scenario 3: A left-right saw-tooth input – steeringScenario=3, steerAngle=15 
% steeringScenario = 3;   
% steerAngle = 15;        % steer angle in [deg] for step input
% enablePlotting = false; % turn on the plotting features (off is faster!)
% maxYawRate = 20;        % [rad/s] -- stop the simulation if yaw rate is >
% 
% %calling the Simulink model, and displaying it
% HW4out1 = sim('HW4_Task2.slx',4);  
% 
% % ploting the Steering angle 
% figure(11)
% plot(HW4out1.steer);
% grid on
% title('Scenario 3: Steering Angle with Time')
% xlabel('Time (s)')
% ylabel('Steering Angle (degrees)')
% 
% % ploting the Yaw Rate
% figure(12)
% plot(HW4out1.yawrate);
% grid on
% title('Scenario 3: Yaw Rate (degrees/s) with Time')
% xlabel('Time (s)')
% ylabel('Yaw Rate (degrees/s)')
% 
% % ploting the Vehicle Position
% figure(13)
% plot(HW4out1.X.Data,HW4out1.Y.Data);
% grid on
% title('Scenario 3: Vehicle Position Y vs X')
% xlabel('X')
% ylabel('Y')
% 
% %% Task 4
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Task 4.1



% % calculating the simulated turning radius
Turning_data=[HW4out1.X.Data HW4out1.Y.Data];
[Lc,R,K] = curvature(Turning_data); % external function
% mean(R(50:end-1,:))                % excluding the first 50 points before the steady speed
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Task 4.2
% 
% i=1;
% V = [1, 30, 60, 90, 120, 140, 150, 160, 170, 180]./3.6;
% Radiud_per_velocity=zeros(length(V),2);
% for V=V 
%     % calling the Simulink model, and displaying it
%     HW4out1 = sim('HW4_Task2.slx',15);  
% 
%     % ploting the Vehicle Position
%     figure(16)
%     plot(HW4out1.X.Data,HW4out1.Y.Data,'LineWidth',1);
%     grid on
%     title('Task 4.2: Vehicle Position - Ackermann steering angle at different speeds')
%     xlabel('X')
%     ylabel('Y')
%     hold on
%     
%     % calculating the simulated turning radius for each velocity
%     Turning_data=[HW4out1.X.Data HW4out1.Y.Data];
%     [Lc,R,K] = curvature(Turning_data); % external function
%     m=mean(R(50:end-1,:));          % excluding the first 50 points before the steady speed
%     Radiud_per_velocity(i,1)=V;
%     Radiud_per_velocity(i,2)=m;
%     
%     i=i+1;
% end
% clear i;
% legend('1','30', '60', '90', '120', '140', '150', '160', '170', '180')
% 
% % plotting the simulated turning radius for each velocity
% figure(17)
% axis equal
% plot(Radiud_per_velocity(:,1),Radiud_per_velocity(:,2),'o-');
% grid on
% title('Task 4.2: Average Turning Radius per Velocity')
% xlabel('The Longitudinal  Velocity')
% ylabel('Average Turning Radius')
% 
% % computing understeer gradient
% Kus=((WD*mass)/(2*Caf))-(((1-WD)*mass)/(2*Car))
% if Kus>0
%     Vehicle_condition='The vehicle Understeers'
% else
%     Vehicle_condition='The vehicle Oversteers'
% end
% 
%     figure(19)
%     plot(HW4out1.lateral_acceleration);
%     grid on
%     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % Task 4.3
% 
% Caf=3000;         % Changing front tire cornering stiffness front [N/deg]
% maxYawRate = 20;  % [rad/s] -- stop the simulation if yaw rate is >
% 
% i=1;
% V = [1, 30, 60, 90, 120, 140, 150, 160, 170, 180]./3.6;
% Radiud_per_velocity=zeros(length(V),2);
% for V=V 
%     % calling the Simulink model, and displaying it
%     HW4out1 = sim('HW4_Task2.slx',15);  
% 
%     % ploting the Vehicle Position
%     figure(18)
%     plot(HW4out1.X.Data,HW4out1.Y.Data,'LineWidth',1);
%     grid on
%     title('Task 4.3: Vehicle Position - Ackermann steering angle at different speeds - Caf=3000[N/deg]')
%     xlabel('X')
%     ylabel('Y')
%     hold on
%         
%     i=i+1;
% end
% clear i;
% legend('1','30', '60', '90', '120', '140', '150', '160', '170', '180')
% 
% % computing understeer gradient
% Kus=((WD*mass)/(2*Caf))-(((1-WD)*mass)/(2*Car))
% if Kus>0
%     Vehicle_condition='The vehicle Understeers'
% else
%     Vehicle_condition='The vehicle Oversteers'
% end
% 
% % computing the oversteering critical velocity [km/h]
% The_critical_velocity=sqrt(-L/(Kus*pi/180))*3.6
% 
% 
% %% End
% 

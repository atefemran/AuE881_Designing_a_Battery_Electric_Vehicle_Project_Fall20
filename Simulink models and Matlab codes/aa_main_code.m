%% Group A

%% Initialization
clear all
close all

%% Design Variables
battery_capacity=40000;             %battery capcity in Wh
motorinverter_select='A';           %A, B, or C
configuration='RWD2';               %FWD1, FWD2, RWD1, RWD2, AWD2, or AWD4
drive_ratio=2;                      %final drive ratio
w=0.080;                            %(m) frame rail width
h=0.160;                            %(m) frame rail height
t=0.001;                            %(m) frame rail thickness
frame_rail_area=(w*h)-((w-(2*t))*(h-(2*t)));

if ((configuration =='FWD1') | (configuration =='RWD1'))
    num_motor=1;
    num_gearbox=1;                  %two motors could be connected to the samce GB to double the torque
    num_inverter=1;
    num_differential=1;
elseif ((configuration =='FWD2') | (configuration =='RWD2'))
    num_motor=2;
    num_gearbox=1;                  %two motors could be connected to the samce GB to double the torque
    num_inverter=2;
    num_differential=1;
elseif (configuration =='AWD2')
    num_motor=2;
    num_gearbox=2;
    num_inverter=2;
    num_differential=2;
elseif configuration =='AWD4'
    num_motor=4;
    num_gearbox=2;
    num_inverter=4;
    num_differential=2;
end    


%% Inputs
slope=0;                            %the road slope in %

%% Parameters 
g=9.81;                             %Gravitational acceleration 𝑔 = 9.81 𝑚/𝑠2
Cd=0.27;                            %Drag coefficient Cd = 0.27
r=659.7/1000/2;                     %for Tire size = 225/45𝑅18, the diameter = 659.7 mm
rho=1.26;                           %Density of air 𝜌 = 1.26 𝑘𝑔/𝑚3
steel_density= 7800;                %kg/m^3

frf=1.2;                            %Front ride requency (Hz)
frr=1.44;                           %Front ride requency (Hz)
damping_ratio=0.3;                  %damping ratio

% Dimensions
battery_width=1.25;                 %max battery width  
W103=1820/1000;                     %vehicle width (m)
W101=1550/1000;                     %Tread_front
L101=2700/1000;                     %the wheel base length

if motorinverter_select=='A'        %%motor, inverter, and gearbox paratmers
    motor_mass=85;                  %kg
    inverter_mass=35;               %kg
    motorinverter_cost=7500;        %dollars
    motorinverter_power=175000;     %W; continuous
    gearbox_mass=45;                %high torque gearbox is needed based on the map
    gearbox_cost=7500;              %high torque gearbox is needed based on the map
    motor_cylinder_dia = 0.390;
    inverter_box_h = 0.422;
    gearbox_h = 0.3;   
    load MotorA_Data.mat            %to load motor A map
elseif motorinverter_select=='B'
    motor_mass=41;                  %kg
    inverter_mass=16;               %kg
    motorinverter_cost=5200;        %dollars
    motorinverter_power=70000;      %W; continuous
    gearbox_mass=28;                %low torque gearbox is needed based on the map
    gearbox_cost=4000;              %low torque gearbox is needed based on the map
    motor_cylinder_dia = 0.276;
    inverter_box_h = 0.393;
    gearbox_h = 0.25;
    load MotorB_Data.mat            %to load motor B map
elseif motorinverter_select=='C'
    motor_mass=55;                  %kg
    inverter_mass=12;               %kg
    motorinverter_cost=5500;        %dollars
    motorinverter_power=100000;     %W; continuous
    gearbox_mass=45;                %high torque gearbox is needed based on the map
    gearbox_cost=7500;              %high torque gearbox is needed based on the map
    motor_cylinder_dia = 0.28;
    inverter_box_h = 0.2;
    gearbox_h = 0.3;    
    load MotorC_Data.mat            %to load motor C map
end

%Masses
passenger_mass=101;             %kg/passenger; left per passenger to be used in the CG calculations later
cargo_mass=28;                  %kg (4*7)
info_mass=15;                   %kg
sunroof_mass=10;                %kg
unsprung_mass=60*4;             %kg
total_gearbox_mass=gearbox_mass*num_gearbox;            %kg (based on the selection)
differential_mass=10;                                    %kg
total_differential_mass= num_differential*differential_mass; 
axle_mass=20;                                           %5kg*4halfs
battery_mass=battery_capacity/250;                              %battery capacity (Wh)/250(Wh/kg) = battery mass (kg)
total_motorinverter_mass=num_motor*(motor_mass+inverter_mass); %kg (based on the selection)

%Volumes
battery_volume=(battery_capacity/600/1000);             %m^3
battery_crosssection_area=battery_volume/battery_width;

%Costs
misc_electronics_cost=3000;     
wheel_assembly_cost=450*4;
differntial_cost=600;
axles_cost=350;

%% Assumptions

E=2.15*(10^11);                     %Young’s Modulus Pa
G=78000;                            %Shear Modulus N/mm2

Cr=0.009;                           %Rolling coefficient Cr=0.009
mu=0.9;                             %Coefficient of friction (dry asphalt) = 0.9
BrkGain = 1e4;                      %Max Brake effort gain [N]
eff_dl =0.9;                        %Driveline efficiency [unitless] 

useful_battery_capacity=0.95;       %the useful capacity of the battery is 90% of the capacity

Jz = 2380;                          % Moment of inertia of the vehicle [kgm^2]
Caf = 1290;                         % Tire cornering stiffness front [N/deg]
Car = 3000;                         % Tire cornering stiffness rear [N/deg]
Cmzf = 30;                          % Tire aligning stiffness front [Nm/deg]
Cmzr = 20;                          % Tire aligning stiffness rear [Nm/deg]

MR=0.9;                               % Motion Ratio


%% Simulink models Paramters
time_step = 0.05;

%% 2D Model inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOADS IN FRAME RAIL SECTION NEEDS TO BE UPDATED MANUALLY BASED ON THE
% DESIGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%the below values to be entered from the 2D model and are required for
%calculations, they are just assumed now to run the model
L108=3944/1000;                         %The total length of the vehicle
H101=1420/1000;                         %The total vehicle height
battery_thickness=0.05;                 %battery thickness 
L104 = 0.622;                           %front overhang
H156 = 0.140;                           % ground clearance
H5_1 = 0.500;                           % front passenger H-Point
H5_2 = 0.550;                           % rear passenger H-Point

% X loads locations // needed for the frame rail calculations later
cargofront_X_dis = 0.28;                % Cargo         (front)
motorf_X_dis = 0.4;                     % Motor         (front)       
inverterf_X_dis = 0.4;                  % Inverter      (front)
gearboxf_X_dis = 0.4;                   % Gearbox       (front)
differentialf_X_dis = 0.4;              % Differential  (front)
frontaxle_X_dis = 0.622;                % Axle          (front)
batteryfront_X_dis = 0.807;             % Battery       (front)
infotainment_X_dis = 0.9;               % Infortainment 
frontpass_X_dis = 1.3+0.622;            % Passenger     (front)
sunroof_X_dis = 2.0;                    % Sunroof       
body_X_dis = L108/2;                    % Body          
rearpass_X_dis = 2.04+0.622;            % Passenger     (rear)
motorr_X_dis = 0.622+2.7-0.3;           % Motor         (rear)
inverterr_X_dis = motorr_X_dis;         % Inverter      (rear)
gearboxr_X_dis = motorr_X_dis+(0.25/2); % Gearbox       (rear)
differentialr_X_dis = gearboxr_X_dis;   % Differential  (rear)
rearaxle_X_dis = 0.622+2.7;             % Axle          (rear)
batteryrear_X_dis = 3.4;                % Battery       (rear)
cargorear_X_dis = 3.944-0.305;          % Cargo         (rear)
framerail_X_dis_start = 0;              % Frame Rail
framerail_X_dis_end = L108;             % Frame Rail

% Y loads locations // needed for the frame rail calculations later
cargofront_Y_dis = H156 + 0.436 + 0.225/2 ;         % Cargo         (front)
motorf_Y_dis = H156 + motor_cylinder_dia/2 ;        % Motor         (front) 
inverterf_Y_dis = H156 + inverter_box_h/2 ;         % Inverter      (front)
gearboxf_Y_dis = H156 + gearbox_h/2 ;               % Gearbox       (front)
differentialf_Y_dis = H156 + r ;                    % Differential  (front)
frontaxle_Y_dis = H156 + r ;                        % Axle          (front)
batteryfront_Y_dis = H156 + battery_thickness/2 ;   % Battery       (front)
infotainment_Y_dis = H5_1;                          % Infortainment 
frontpass_Y_dis = H5_1;                             % Passenger     (front)
sunroof_Y_dis = H101;                               % Sunroof      
body_Y_dis = H101/3;                                % Body  
rearpass_Y_dis = H5_2;                              % Passenger     (rear)
motorr_Y_dis =  H156 + motor_cylinder_dia/2 ;       % Motor         (rear)
inverterr_Y_dis = H156 + inverter_box_h/2 ;         % Inverter      (rear)
gearboxr_Y_dis = H156 + gearbox_h/2 ;               % Gearbox       (rear)
differentialr_Y_dis = H156 + r;                     % Differential  (rear)
rearaxle_Y_dis = H156 + r ;                         % Axle          (rear)
batteryrear_Y_dis = H156 + battery_thickness/2 ;    % Battery       (rear)
cargorear_Y_dis = H156 + 0.436 + 0.225/2;           % Cargo         (rear)    
framerail_Y_dis = H156 + h/2;                       % Frame Rail

%% Frame Rail optimization - intialization (comment out after finishing it)
% % before activating add ; to stress and the deflection, and deactivvate
% % plots in SFBM
% frame_rail_optimization=zeros(0,6); %w,aspect ratio,t,area,stress, deflection
% aspect_ratio=1;
% i=0;
% for w=0.04:0.005:0.15
%     for aspect_ratio=1:0.5:5
%         h=w*aspect_ratio;
%         for t=0.001:0.001:0.015
%             i=i+1;
%             frame_rail_area=(w*h)-((w-(2*t))*(h-(2*t)));
%             
%% Vehicle Mass Calculations
% SI units have been used all through out the simulation unless otherwise stated

frame_rail_mass=3*2*frame_rail_area*L108*7800;          %kg (based on the frame rail design and 2D inputs 
body_mass=L108*H101*50;                                 %kg
mass=(4*passenger_mass)+cargo_mass+info_mass+sunroof_mass+unsprung_mass+total_gearbox_mass+total_differential_mass+axle_mass+body_mass+frame_rail_mass+battery_mass+total_motorinverter_mass;

%% Frame Rail calculations
% as per the assumptions in the report the below analysis was performed

% calculating the different loads (all of them in N)
load_cargof= -cargo_mass*g/4;
    if ((configuration =='RWD1') | (configuration =='RWD2'))
     load_motorf=0;
     load_inverterf=0;
     load_gearboxf=0;
     load_differentialf=0;  
     elseif ((configuration =='FWD2') | (configuration =='AWD4'))
     load_motorf= -2*motor_mass*g;
     load_inverterf= -2*inverter_mass*g;
     load_gearboxf= -1*gearbox_mass*g;
     load_differentialf= -1*differential_mass*g; 
     elseif ((configuration =='FWD1') | (configuration =='AWD2'))
     load_motorf= -1*motor_mass*g;
     load_inverterf= -1*inverter_mass*g;
     load_gearboxf= -1*gearbox_mass*g;
     load_differentialf= -1*differential_mass*g;
     end
load_axlef= -axle_mass*g/2;
load_batteryf= -battery_mass*g;
load_infotainment= -info_mass*g;
load_passf= -2*passenger_mass*g;
load_sunroof= -sunroof_mass*g;
load_body=-body_mass*g;
load_passr= -2*passenger_mass*g;
     if ((configuration =='RWD1') | (configuration =='AWD2'))
     load_motorr= -1*motor_mass*g;
     load_inverterr= -1*inverter_mass*g;
     load_gearboxr= -1*gearbox_mass*g;
     load_differentialr= -1*differential_mass*g;
     elseif ((configuration =='FWD2') | (configuration =='FWD1'))
     load_motorr=0;
     load_inverterr=0;
     load_gearboxr=0;
     load_differentialr=0;
     elseif ((configuration =='RWD2') | (configuration =='AWD4'))
     load_motorr= -2*motor_mass*g;
     load_inverterr= -2*inverter_mass*g;
     load_gearboxr= -1*gearbox_mass*g;
     load_differentialr= -1*differential_mass*g; 
     end
load_axler= -axle_mass*g/2;
load_batteryr= -battery_mass*g/4;
load_cargor= -cargo_mass*g;
load_framerail= -frame_rail_mass*g/L108;


% The SFBM function
Name1 = 'Frame Rail';
LengthSupport1 = [L108,L104,(L101+L104)];       % Length and Supports
% Concetrated Loads

F1 = {'CF',load_cargof,cargofront_X_dis};  
F2 = {'CF',load_motorf,motorf_X_dis};
F3 = {'CF',load_inverterf,inverterf_X_dis};
F4 = {'CF',load_gearboxf,gearboxf_X_dis};
F5 = {'CF',load_differentialf,differentialf_X_dis};  
F6 = {'CF',load_axlef,frontaxle_X_dis};  
F7 = {'CF',load_batteryf,batteryfront_X_dis};  
F8 = {'CF',load_infotainment,infotainment_X_dis};  
F9 = {'CF',load_passf,frontpass_X_dis};  
F10 = {'CF',load_sunroof,sunroof_X_dis};  
F11 = {'CF',load_body,body_X_dis};  
F12 = {'CF',load_passr,rearpass_X_dis};  
F13 = {'CF',load_motorr,motorr_X_dis};  
F14 = {'CF',load_inverterr,inverterr_X_dis};
F15 = {'CF',load_gearboxr,gearboxr_X_dis}; 
F16 = {'CF',load_differentialr,differentialr_X_dis};
F17 = {'CF',load_axler,rearaxle_X_dis};  
F18 = {'CF',load_batteryr,batteryrear_X_dis};  
F19 = {'CF',load_cargor,cargorear_X_dis};  

% Distributed Loads
D1 = {'DF',load_framerail,[framerail_X_dis_start,framerail_X_dis_end]}; 

% Call the function - output in the max bending mometn - N.m
bending_moment=SFBM(Name1,LengthSupport1,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18,F19,D1);
bending_moment_max=max(abs(bending_moment));                            %max bending moment

% Calculating the max stress on the beam (the hieghest moment at the most
% far line of the section (the biggest stresses)
second_moment_area=((w*(h^3))/12)-(((w-2*t)*((h-2*t)^3))/12);
stress=(bending_moment_max*(h/2))/(second_moment_area*1000000)          %MPa

% Calculating the deflection 
% rough assumption of having the vehicle weight affecting at the center
% distance betweeen the two wheels (the most sever condition)
% two supporters and one concentrated load at the middle equals (m.g)

deflection_max=(mass-unsprung_mass)*g*(L101^3)/(48*E*second_moment_area)

%% Frame Rail optimization - end (comment out after finishing it)
% 
% %w,aspect ratio,t,area,stress, deflection
% frame_rail_optimization(i,1)=w;
% frame_rail_optimization(i,2)=aspect_ratio;
% frame_rail_optimization(i,3)=t;
% frame_rail_optimization(i,4)=frame_rail_area;
% frame_rail_optimization(i,5)=stress;
% frame_rail_optimization(i,6)=deflection_max;
%         end
%     end
% end

%% Center of Gravity

mass_x_dis = -1*(load_cargof*cargofront_X_dis ...
    + load_motorf*motorf_X_dis ...
    + load_inverterf*inverterf_X_dis ...
    + load_gearboxf*gearboxf_X_dis ...
    + load_differentialf*differentialf_X_dis ...
    + load_axlef*frontaxle_X_dis ...
    + load_batteryf*batteryfront_X_dis ...
    + load_infotainment*infotainment_X_dis ...
    + load_passf*frontpass_X_dis ...
    + load_sunroof*sunroof_X_dis ...
    + load_body*body_X_dis ...
    + load_passr*rearpass_X_dis ...
    + load_motorr*motorr_X_dis ...
    + load_inverterr*inverterr_X_dis ...
    + load_gearboxr*gearboxr_X_dis ...
    + load_differentialr*differentialr_X_dis ...
    + load_axler*rearaxle_X_dis ...
    + load_batteryr*batteryrear_X_dis ...
    + load_cargor*cargorear_X_dis ...
    +(load_framerail*L108)*(L108/2))/g;

mass_CG=mass-unsprung_mass;

CG_a = (mass_x_dis/mass_CG) - L104 ;  %L104: front overhang
CG_b = L101 - CG_a;
CG_location=CG_a/L101


mass_Y_dis = -1*(load_cargof*cargofront_Y_dis ...
    + load_motorf*motorf_Y_dis ...
    + load_inverterf*inverterf_Y_dis ... 
    + load_gearboxf*gearboxf_Y_dis ...
    + load_differentialf*differentialf_Y_dis ...
    + load_axlef*frontaxle_Y_dis ...
    + load_batteryf*batteryfront_Y_dis ...
    + load_infotainment*infotainment_Y_dis ...
    + load_passf*frontpass_Y_dis ...
    + load_sunroof*sunroof_Y_dis ...
    + load_body*body_Y_dis ...
    + load_passr*rearpass_Y_dis ...
    + load_motorr*motorr_Y_dis ...
    + load_inverterr*inverterr_Y_dis ...
    + load_gearboxr*gearboxr_Y_dis ...
    + load_differentialr*differentialr_Y_dis ...
    + load_axler*rearaxle_Y_dis  ...
    + load_batteryr*batteryrear_Y_dis ...
    + load_cargor*cargorear_Y_dis ...
    + (load_framerail*L108)*framerail_Y_dis)/g;

CG_h = mass_Y_dis/mass_CG;

%% Vehicle Costing Calculations
body_cost=1000*W103*L108*H101;
frame_rail_cost=frame_rail_mass*2;
battery_cost=145*battery_capacity/1000;

cost=misc_electronics_cost+wheel_assembly_cost+(gearbox_cost*num_gearbox) ...
    +differntial_cost+axles_cost+(num_motor*motorinverter_cost)+body_cost ...
    +frame_rail_cost+battery_cost

%% Vehicle Perfromance Model Setting and Parameters
%road load model
vehicle_frontal_area=H101*W103;

%Forces Distributer Model
front_on=0;       %0 means no forces to the front, 1 means forces are going to front
rear_on=0;        %0 means no forces to the rear, 1 means forces are going to rear        
if ((configuration =='FWD1') | (configuration =='FWD2'))
    front_on=1;       %0 means no forces to the front, 1 means forces are going to front
    rear_on=0;        %0 means no forces to the rear, 1 means forces are going to rear      
elseif ((configuration =='RWD2') | (configuration =='RWD1'))
    front_on=0;       %0 means no forces to the front, 1 means forces are going to front
    rear_on=1;        %0 means no forces to the rear, 1 means forces are going to rear      
elseif ((configuration =='AWD4') | (configuration =='AWD2'))
    front_on=1;       %0 means no forces to the front, 1 means forces are going to front
    rear_on=1;        %0 means no forces to the rear, 1 means forces are going to rear      
end

%Torque multiplication based on the number of engines per axle
torque_multi_factor=0;
if (configuration =='AWD2')
    torque_multi_factor=1;   
elseif (configuration =='FWD1')
    torque_multi_factor=1;   
elseif (configuration =='RWD1')
    torque_multi_factor=1; 
elseif (configuration =='RWD2')
    torque_multi_factor=2; 
elseif (configuration =='AWD4')
    torque_multi_factor=2; 
elseif (configuration =='FWD2')
    torque_multi_factor=2;         
end

%Electrical Motor Model
%Breakpoint values
omega_BP=max(size(EM_Omega_Map));
omega_max_BP=max(size(EM_Omega_Max));
torque_BP=max(size(EM_Torque_Map));
torque__max_BP=max(size(EM_Torque_Max));
%max Torqua and omega
omega_max=max(EM_Omega_Max);
torque__max=max(EM_Torque_Max);


%Drive Cycle Model
load US06_Drive_Cycle.mat
cycle(:,1)=t_cyc;
cycle(:,2)=v_cyc;
Kp = 0.5;
Ki = 0.03;

%% Max speed calculations (1 of 3 don't change the sequence)

sim_max=sim('Vehicle_Performance_Model_max.slx');
max_speed=max(sim_max.max_velocity_km_per_h)                  %(m/s)

%% 0-100 km/h time calculations (2 of 3 don't change the sequence)

acceleration_time=max(sim_max.acceleration_time_to_hundred_kmh) %(s)

%% Range Calculations  calculations (3 of 3 don't change the sequence)

sim_range=sim('Vehicle_Performance_Model_range.slx');
range_battery=max(sim_range.range_km)                 %(km)

%% Ascend a 12% slope speed (minimum at 30 km/h) 

slope=0.12;                            %the road slope in %
sim_slope=sim('Vehicle_Performance_Model_slope.slx');

slope_difference_variation=sim_slope.slope_error.Data(end);

if abs(slope_difference_variation)<0.001
    slope_result='[Pass]'
else
    slope_result='[Fail]'
end

%% Profit calculations

profit_input=[cost; range_battery; acceleration_time; max_speed];
[profit,demand,price] = profitPredict(profit_input)

%% Assesing the stability at 100 m radius with up to 0.8 g lateral accleration

Kus=(((1-CG_location)*mass)/(2*Caf))-(((CG_location)*mass)/(2*Car))

steeringScenario = 1;
steerAngle = 6;   % Ackermann steering angle in [deg] 
enablePlotting = false;        % turn on the plotting features (off is faster!)
maxYawRate = 20;               % [rad/s] -- stop the simulation if yaw rate is >

%Intial values
V = 30;                     % velocity [m/s] (1 kmph is 40/3.6 m/s)
turning_data = sim('vehicle_dynamics_turning.slx',200); 
Turning_data=[turning_data.X.Data turning_data.Y.Data]; 
[Lc,R,K] = curvature(Turning_data); % external function 
radius=mean(R(50:end-1,:)); % excluding the first 50 points before the steady speed

%Determining the speed for 100 m diamter (50 m radius)
while (radius<49.5) || (50.5<radius) 
    if radius<50
        V=V+0.5;
    elseif radius>50
        V=V-0.5;
    end
    turning_data = sim('vehicle_dynamics_turning.slx',200); 
    Turning_data=[turning_data.X.Data turning_data.Y.Data]; 
    [Lc,R,K] = curvature(Turning_data); % external function 
    radius=mean(R(50:end-1,:)); % excluding the first 50 points before the steady speed
end

% calculating the simulated turning radius 
turning_data = sim('vehicle_dynamics_turning.slx',1000);  
Turning_data=[turning_data.X.Data turning_data.Y.Data]; 
[Lc,R,K] = curvature(Turning_data); % external function 
radius=mean(R(50:end-1,:)) % excluding the first 50 points before the steady speed

lateral_acceleration_max=max(abs(turning_data.lateral_acceleration.data))

if lateral_acceleration_max>=(g*0.8)
    turnning_result='[Pass]'
elseif lateral_acceleration_max<(g*0.8)
    turning_result='[Adjust]'    %Increase the speed and steering angle
end

figure(1)
plot(turning_data.X.Data,turning_data.Y.Data);
grid on
title('Turning Radius')
xlabel('X')
ylabel('Y')

%% Spring Constant Calculations

front_corners_sprung_mass=(mass-unsprung_mass)*(1-CG_location)/2;
rear_corners_sprung_mass=(mass-unsprung_mass)*(CG_location)/2;

spring_rate_front=4*(pi^2)*(frf^2)*front_corners_sprung_mass*(MR^2)     % N
spring_rate_rear=4*(pi^2)*(frr^2)*rear_corners_sprung_mass*(MR^2)       % N


%% Damping Coefficient 

spring_damping_coeffiecient_front = 4*pi*damping_ratio*frf*front_corners_sprung_mass;   %N.s/m
spring_damping_coeffiecient_rear = 4*pi*damping_ratio*frr*rear_corners_sprung_mass;     %N.s/m

%%


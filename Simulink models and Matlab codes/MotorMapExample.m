% This script illustrates what the data is in each Motor_Data.mat file:
% MotorA_Data.mat and MotorB_Data.mat
% clear all
% close all

%%
load MotorA_Data.mat

% the efficiency values for which to create contour lines
eta = [0.97 0.965 0.96 0.955 0.95 0.945 0.94 0.93 0.92 0.91 0.9 0.88 0.86 0.84 0.82 0.8 0.75 0.7 0.65 0.6 0.55 0.5];

figure (1)

% Create a contour plot of the efficiency map:
%   EM_Efficiency_Map is a 2D lookup table with efficiency values for given
%      rotational speed and torques corresponding to the scales below
%   EM_Omega_Map represents the speed in rad/s
%   EM_Toreque_Map represents the torque in Nm
[cc,h] = contour(EM_Omega_Map,EM_Torque_Map,EM_Efficiency,eta);
clabel(cc,h);
hold on
grid on

% The map extends beyond the rated range for the motor.  We therefore plot
% the maximum torque values and white out the portion of the map that is
% invalid 
plot([EM_Omega_Max EM_Omega_Max(end)],[EM_Torque_Max 0],'k','linewidth',3)
colormap(cool)
hold off

large = 1e6;
xp = [EM_Omega_Max EM_Omega_Max(end) large large 0];
yp = [EM_Torque_Max 0 0 large large];
patch(xp,yp,'white','Edgecolor','white');

% make it look good...
axis([0 max(EM_Omega_Map)+50 0 max(EM_Torque_Map)+20])
xlabel('Motor Speed [rad/s]')
ylabel('Motor Torque [Nm]')
title('Motor A Efficiency Map')
legend('Efficiency Contours','Motor Peak Torque Line','location','northeast')

% %%
% load MotorB_Data.mat
% 
% % the efficiency values for which to create contour lines
% eta = [0.97 0.965 0.96 0.955 0.95 0.945 0.94 0.93 0.92 0.91 0.9 0.88 0.86 0.84 0.82 0.8 0.75 0.7 0.65 0.6 0.55 0.5];
% 
% figure (2)
% 
% % Create a contour plot of the efficiency map:
% %   EM_Efficiency_Map is a 2D lookup table with efficiency values for given
% %      rotational speed and torques corresponding to the scales below
% %   EM_Omega_Map represents the speed in rad/s
% %   EM_Toreque_Map represents the torque in Nm
% [cc,h] = contour(EM_Omega_Map,EM_Torque_Map,EM_Efficiency,eta);
% clabel(cc,h);
% hold on
% grid on
% 
% % The map extends beyond the rated range for the motor.  We therefore plot
% % the maximum torque values and white out the portion of the map that is
% % invalid 
% plot([EM_Omega_Max EM_Omega_Max(end)],[EM_Torque_Max 0],'k','linewidth',3)
% colormap(cool)
% 
% large = 1e6;
% xp = [EM_Omega_Max EM_Omega_Max(end) large large 0];
% yp = [EM_Torque_Max 0 0 large large];
% patch(xp,yp,'white','Edgecolor','white');
% 
% % make it look good...
% axis([0 max(EM_Omega_Map)+50 0 max(EM_Torque_Map)+20])
% xlabel('Motor Speed [rad/s]')
% ylabel('Motor Torque [Nm]')
% title('Motor B Efficiency Map')
% legend('Efficiency Contours','Motor Peak Torque Line','location','northeast')
% 
% %%
% load MotorC_Data.mat
% 
% % the efficiency values for which to create contour lines
% eta = [0.97 0.965 0.96 0.955 0.95 0.945 0.94 0.93 0.92 0.91 0.9 0.88 0.86 0.84 0.82 0.8 0.75 0.7 0.65 0.6 0.55 0.5];
% 
% figure (3)
% 
% % Create a contour plot of the efficiency map:
% %   EM_Efficiency_Map is a 2D lookup table with efficiency values for given
% %      rotational speed and torques corresponding to the scales below
% %   EM_Omega_Map represents the speed in rad/s
% %   EM_Toreque_Map represents the torque in Nm
% [cc,h] = contour(EM_Omega_Map,EM_Torque_Map,EM_Efficiency,eta);
% clabel(cc,h);
% hold on
% grid on
% 
% % The map extends beyond the rated range for the motor.  We therefore plot
% % the maximum torque values and white out the portion of the map that is
% % invalid 
% plot([EM_Omega_Max EM_Omega_Max(end)],[EM_Torque_Max 0],'k','linewidth',3)
% colormap(cool)
% 
% large = 1e6;
% xp = [EM_Omega_Max EM_Omega_Max(end) large large 0];
% yp = [EM_Torque_Max 0 0 large large];
% patch(xp,yp,'white','Edgecolor','white');
% 
% % make it look good...
% axis([0 max(EM_Omega_Max)+50 0 max(EM_Torque_Max)+20])
% xlabel('Motor Speed [rad/s]')
% ylabel('Motor Torque [Nm]')
% title('Motor C Efficiency Map')
% legend('Efficiency Contours','Motor Peak Torque Line','location','northeast')
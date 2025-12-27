% Stability Thing
% Authors: Byrn Balangauan, Alex Burghardt

% Description: 
% This software takes the geometric parameters of the various
% flying surfaces of the airplane in order to build up the stability analysis
% coefficients needed to analyze the vehicle as a full assembly.

% Comments: 
% It currently only works for normalized tails, meaning that it will not
% work for a plane that uses a V-tail. This is being developed further.

folder = fileparts(which("new_base_code_pt_2.m"));
addpath(genpath(folder));

clear

% Further Author Notes:
% 

%% Flight Conditions
h = 30000; % altitude, ft
g = 32.2; % gravity acceleration, ft/s^2
M = 0.459; % Mach number, N/A

% Derived From Flight Conditions
c = 994.8; % speed of sound, from ISA, ft/s
rho = 8.893e-4; % air density, sl/ft^3
W = 6360; % weight of aircraft, lbf
mgc = 5.47; % mean geometric chord, ft
b = 33.8; % wing span, ft
S_w = 182; % ft^2

% Other relevant values
u1 = c * M; % trim velocity
q1 = 0.5 * rho * u1 ^ 2; % trim dynamic pressure
m  = W / g; % mass

% Aircraft Stability Characteristics

%% Thrust Coefficients

% Variation of airplane X-axis thrust coefficient with dimensionless speed
CTxu = -0.07;

% Steady state thrust coefficient (usually equal to steady state drag)
CTx1 = 0.030;

%% Drag Coefficients

% Variation of airplane drag coefficient with dimensionless speed
CDu = 0;

% Steady state drag coefficient (Induced drag)
CD1 = 0.03;

% Variation of airplane drag coefficient with angle of attack
CDa = 0.250;

% Variation of airplane drag coefficient with elevator deflection angle
CDde = 0;

%% Lift Force Coefficients

% Variation of airplane lift coefficient with dimensionless speed
CLu = 0;

% Variation of airplane lift coefficient with angle of attack
CLa = 5.15;

% Steady state lift coefficient
CL1 = 0.378;

% Variation of airplane lift coefficient with elevator deflection
CLde = 0.5;

% Variation of airplane lift coefficient with dimensionless pitch rate 
CLq = 4.1;

%% Pitching Moment Coefficients

% Variation of pitching moment coefficient with dimensionless speed
CMu = 0;

% Steady state pitching moment coefficient
CM1 = 0;

% Variation of pitching moment coefficient with angle of attack
CMa = -0.7;

% Variation of pitching moment coefficient with dimensioness pitch rate
CMq = -14.9;

% Variation of pitching moment coefficient with elevator deflection
CMde = -1.12;

%% Side Force Coefficients

% Variation of side force coefficient with sideslip angle
CYb = -0.346;

% Variation of side force coefficient with dimensionless rate of change of
% roll rate
CYp = -0.0827;

% Yariation of side force coefficient with dimensionless rate of change of
% yaw rate
CYr = 0.3;

% Variation of side force coefficient with rudder deflection
CYdr = 0.2;

% Variation of side force coefficient with aileron deflection
CYda = 0;

%% Rolling Moment Coefficients

% Variation of rolling moment coefficient with sideslip angle
Clb = -0.0944;

% Variation of rolling moment coefficient with rate of change of roll rate
Clp = -0.442;

% Variation of rolling moment coefficient with rate of change of yaw rate
Clr = 0.0926;

% Variation of rolling moment coefficient with aileron deflection
Clda = 0.1810;

% Variation of rolling moment coefficient with rudder deflection
Cldr = 0.0150;

%% Yawing Moment Coefficients

% Variation of yawing moment coefficient with sideslip angle
CNb = 0.1106;

% Variation of yawing moment coefficient due to thrust with sideslip angle
CNTb = 0;

% Variation of yawing moment coefficient due to rate of change of roll rate
CNp = -0.0243;

% Variation of yawing moment coefficient due to rate of change of yaw rate
CNr = -0.1390;

% Variation of yawing moment coefficient due to aileron deflection
CNda = -0.0254;

% Variation of yawing moment coefficient due to rudder deflection
CNdr = -0.0365;

%% Important Inertia Values
Ixx = 7985
Iyy = 3326
Izz = 11183

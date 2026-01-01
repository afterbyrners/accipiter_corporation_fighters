%%% FORCE BUILD-UP

% Description: This program uses the geometric and aerodynamic parameters
% of the subject aircraft in order to compute the necessary forces and
% moments needed to properly analyze the aircraft.

% Afterwards, these stability derivatives can then be implemented in a
% separate but accompanying program to assemble a perturbation apprehensive
% state-space model of the aircraft.

% It currently only works for normalized tails, meaning that it will not
% work for a plane that uses a V-tail. This is being developed further.

folder = fileparts(which("new_base_code_pt_1.m"));
addpath(genpath(folder));

clear

%% FUNDAMENTAL AXIS NOTES (For a later README)

% This presumes the Stability Axis System (SAS) of the aircraft, which is
% fundamentally a modification of the Body Axis System (BAS) of the normal
% aircraft with a rotation of the X and Z axes by the trim angle of attack.

% In the BAS:
% The X-axis is +X toward the nose & -X toward the aft of the plane
% The Y-axis is +Y toward the right & -Y toward the left of the plane
% The Z-axis is +Z toward the bottom & -Z toward the top of the plane

% In SAS, the sideslip is still considered, but the perturbations in the Y
% axis remain so small that only the sideslip is measured and the velocity
% projected onto the stability axis Xs, U1s, is presumed to be the total
% airspeed of the vehicle, VP1. Also, the velocity measured on the Zs axis
% is considered to be so small, that it is not to be measured. 

%% NOTATIONS

% X  : Axis running lengthwise with the plane, +X at nose, -X at aft
% Y  : Axis running spanwise with the plane, +Y at right, -Y at left
% Z  : Axis running depthwise with the plane, +Z at bottom, -Z at top

% θ  : Pitch angle/attitude of the plane, called "theta" usually
% Φ  : Bank/roll angle/attitude of the plane, called "phi" usually
% Ψ  : Heading angle/attitude of the plane, called "psi" usually
% α  : Angle of attack of the plane, called "alpha" or "AOA"
% β  : Sideslip angle of the plane, called "beta"

% VP : True aircraft air speed 
% U  : Component of VP along X
% V  : Component of VP along Y
% W  : Component of VP along Z

% FaX: Aerodynamic force component about the X axis
% FaY: Aerodynamic force component about the Y axis
% FaZ: Aerodynamic force component about the Z axis

% l  : Rolling moment made about the X axis
% m  : Pitching moment made about the Y axis
% n  : Yawing moment made about the Z axis

% P  : Roll rate
% Q  : Pitch rate
% R  : Yaw rate

% D : Coefficient of drag
% L : Coefficient of lift
% T : Coefficient of thrust (not an aerodynamic force)

% iw : Incidence angle of wing planform
% ih : Incidence angle of horizontal flight surface
% iv : Incidence angle of vertical flight surface

% δe : Deflection of elevators 
% δa : Deflection of ailerons
% δr : Deflection of rudders
% δf : Deflection of flaps
% δv : Deflection of ruddervators
% δn : Deflection of thruster nozzles

%% Plane & Flight Characteristics

% Define Altitude
h_ft = 30000; % ft
h_m = h_ft * 0.3048; % Convert altitude from feet to meters

% Get ISA properties (SI units)
[T_K, a_mps, P_Pa, rho_kgm3] = atmosisa(h_m);

% Convert to imperial units
T_R = T_K * 9/5;               % Kelvin → Rankine
a_fps = a_mps * 3.28084;       % m/s → ft/s
P_psf = P_Pa * 0.0208854342;   % Pa → lb/ft²
rho_sl = rho_kgm3 * 0.062428;  % kg/m³ → slugs/ft³

M_1 = 0.459; % Mach number

%u1 = c * M; % trim velocity
u_1 = 456;

%q1 = 0.5 * rho * u1 ^ 2; % trim dynamic pressure
q_1 = 92.7;

S_w = 182; % ft^2

% Maximum Thrust Vector Angle in both directions
dt_max_deg = 20; % deg
dt_max_rad = deg2rad(dt_max_deg); % rad

% It is worth noting that small angle approximation can be used such that
% Fz = Tsin(dt_max) = T*dt_max, since sin(dt_max) rough equals dt_max if
% dt_max is around less than 30 deg. In this sense, Fx = Tcos(dt_max) = T
% since cos(dt_max) is roughly about 0.94. % Linearizations can be made.

%% D Forces (Drag)

% Calculate Trim Drag Coefficient
C_D_x_1 = 0.03; % This value makes sense for an f-22. 

% Calculate Variation of Drag Coefficient with Angle of Attack
C_D_a = 0.07;

%% T Forces (Thrust)

% Calculate Trim Thrust Coefficient
C_T_x_1 = C_D_x_1;

% Calculate Variation of Thrust Coefficient with Forward Speed Perturbation
C_T_x_u = M_1 / (q_1 * S_w) - (2 * C_T_x_1); % Roskam Equation 3.206

% Calculate Variation of X Thrust with Nozzle Deflection
C_T_x_d_n = C_T_x_1 * (cos(dt_max_rad)) / dt_max_rad;

% Calculate Variation of Z Thrust with Nozzle Deflection
C_T_z_d_n = C_T_x_1 * (sin(dt_max_rad)) / dt_max_rad;

%% L Forces (Lift)

%% X Forces (Drag & Thrust)

%% Calculate Y Forces (Sideforces)

%% Calculate Z Forces (Lift)

%% Calculate L Moments (Rolling)

%% Calculate M Moments (Pitching)

%% Calculate N Moments (Yawing)

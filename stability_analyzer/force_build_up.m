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

%% Flight Conditions
M = 0.85;
eta_h = 0.95; % Horizontal Tail Dynamic Pressure Ratio
eta_v = 0.95; % Vertical Tail Dynamic Pressure Ratio

%% Airfoil Parameters (NACA 64-206)
a_L_0   = -1.0; % Zero-Lift Angle of Attack, deg
a_L_0   = deg2rad(a_L_0); % Zero-Lift Angle of Attack, rad
C_m_ac  = -0.040; % Pitching Moment at Zero-Lift Coefficient

%% Key Parameters Related to General Vehicle
X_cg = 22; % Actual Location of CG, ft

%% Wing Design Parameters
S_w = 656; % Wing Reference Area, ft^2
AR_w = 2; % Wing Aspect Ratio, N/A
sweep_w = 40; % Wing Leading Edge Sweep, deg
taper_w = 0.08; % Wing Taper Ratio, 
cf_c_w = 0.5; % Wing Chord-Flap Ratio
x_w = 12; % Absolute Wing Front Point X Location
y_w = 0; % Absolute Wing Front Point Y Position
z_w = 1; % Absolute Wing Front Point Z Position

% Convert Wing Leading Edge Sweep to radians
sweep_w = deg2rad(sweep_w); % rad

% Calculate Wing Span
b_w = sqrt(S_w * AR_w); % ft

% Calculate Wing Root Chord
c_r_w = (2 * S_w) / (b_w * (1 + taper_w)); % ft

% Calculate Wing Tip Chord
c_t_w = c_r_w * taper_w; % ft

% Calculate Wing Mean Geometric Chord (MGC)
mgc_w = 2/3 * c_r_w * (1 + taper_w + taper_w ^ 2) / (1 + taper_w); % ft

% Calculate Y-position of Wing MGC
y_mgc_w = (b_w / 6) * (1 + 2 * taper_w) / (1 + taper_w); % ft

% Calculate X-position of Wing MGC
x_mgc_w = y_mgc_w * tan(sweep_w); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_w = TransonicAppendixB(AR_w,sweep_w,taper_w,M); % 1/rad

% Calculate Effectiveness of Wing Control Surfaces (Ailerons)
tau_w = AppendixD(cf_c_w); % N/A

% Calculate Non-Dimensional Wing Aerodynamic Center Position X-Wise
x_ac_w = ((x_w + x_mgc_w + 0.25 * mgc_w) - (x_w + x_mgc_w)) / mgc_w; % N/A

% Calculate Non-Dimensional Center of Gravity Position X-Wise
x_cg = ((X_cg) - (x_w + x_mgc_w)) / mgc_w; % N/A

%% Horizontal Tail (H.T.) Design Parameters
S_h = 77.5123; % H.T. Reference Area, ft^2
AR_h = 2.1803; % H.T. Aspect Ratio, N/A
sweep_h = 39; % H.T. Leading Edge Sweep, deg
taper_h = 0.08409; % H.T. Taper Ratio, 
cf_c_h = 1.0; % H.T. Chord-Flap Ratio
x_h = 39; % Absolute H.T. Front Point X Location
y_h = 0; % Absolute H.T. Front Point Y Position
z_h = 1; % Absolute H.T. Front Point Z Position

% Convert leading edge sweep to radians
sweep_h = deg2rad(sweep_h); % rad

% Calculate Wing Span
b_h = sqrt(S_h * AR_h); % ft

% Calculate Root Chord
c_r_h = (2 * S_h) / (b_h * (1 + taper_h)); % ft

% Calculate Tip Chord
c_t_h = c_r_h * taper_h; % ft

% Calculate Mean Geometric Chord (MGC)
mgc_h = 2/3 * c_r_h * (1 + taper_h + taper_h ^ 2) / ... 
    (1 + taper_h); % ft

% Calculate Y-position of MGC
y_mgc_h = (b_h / 6) * (1 + 2 * taper_h) / (1 + taper_h); % ft

% Calculate X-position of MGC
x_mgc_h = y_mgc_h * tan(sweep_h); % ft

% Calculate Variation of Panel Lift Coefficient with Angle of Attack (AOA)
C_L_a_h = AppendixB(AR_h,sweep_h,taper_h,M); % 1/rad

% Calculate Effectiveness of H.T. Control Surfaces (Elevator)
tau_h = AppendixD(cf_c_h); % N/A

% Calculate Non-Dimensional H.T. Aerodynamic Center Position X-Wise
x_ac_h = ((x_h + x_mgc_h + 0.25 * mgc_h) - (x_w + x_mgc_w)) / mgc_w; % N/A

% Calculate X Distance Between ACs of H.T. and Wing
x_w_h = x_h - (c_r_w / 4) + (c_r_h / 4);

% Calculate Non-Dimensional X Distance between H.T. and Wing
r_h = x_w_h / (0.5 * b_w);

% Calculate Z Distance Between ACs of Horizontal Tail and Wing 
z_w_h = (z_h + y_mgc_h * sqrt(2)) - z_w;

% Calculate Non-Dimensional Z Distance between H.T. and Wing
m_h = z_w_h / (0.5 * b_w);

% Calculate downwash gradient for horizontal tail equivalent
dw_h = TransonicAppendixC(AR_w,sweep_w,taper_w,M,r_h,m_h);

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

% Calculate Level Flight Coefficient of Lift
C_L_0 = -(C_L_a_w * a_L_0);  % N/A

% Calculate Coefficient of Lift caused by General Plane
C_L = C_L_a_w + eta_h * (S_h / S_w) * C_L_a_h * (1 - dw_h);  % 1/rad

% Calculate Coefficient of Lift caused by Horizontal Stabilizer Incidence
C_L_ih = eta_h * (S_h / S_w) * C_L_a_h;  % 1/rad

% Calculate Coefficient of Lift caused by Elevator Flaps
C_L_de = C_L_ih * tau_h;  % 1/rad

%% X Forces (Drag & Thrust)

%% Calculate Y Forces (Sideforces)

%% Calculate Z Forces (Lift)

%% Calculate l Moments (Rolling)

%% Calculate m Moments (Pitching)

% Calculate Level Flight Coefficient of Moment
C_m_0 = C_m_ac + C_L_0 * (x_cg - x_ac_w);  % N/A

% Calculate Coefficient of Moment caused by General Plane
C_m_a = C_L_a_w * (x_cg - x_ac_w) - ...
        C_L_a_h * eta_h * (S_h / S_w) * ...
        (1 - dw_h) * (x_ac_h - x_cg); % 1/rad

% Calculate Coefficient of Moment caused by Horizontal Stabilizer Incidence
C_m_ih = -C_L_a_h * eta_h * (S_h / S_w) * (x_ac_h - x_cg); % 1/rad

% Calculate Coefficient of Moment caused by Flap Angles
C_m_de = -C_L_a_h * eta_h * (S_h / S_w) * tau_h * (x_ac_h - x_cg); % 1/rad

%% Calculate n Moments (Yawing)

%% Calculate Flying Qualities

% Calculate Neutral Point
x_np = (x_ac_w + ((C_L_a_h / C_L_a_w) * eta_h * (S_h / S_w) * ... 
       (1 - dw_h) * x_ac_h)) / (1 + (C_L_a_h / C_L_a_w) * ... 
       eta_h * (S_h / S_w) * (1 - dw_h));

X_np = x_np * mgc_w + x_mgc_w + x_w;

% Calculate Static Margin Percentage
SM = (x_np - x_cg) * 100;
% Conceptual Stability Analysis 

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
% Currently for a Roskam plane

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
%u1 = c * M; % trim velocity
u1 = 456;
%q1 = 0.5 * rho * u1 ^ 2; % trim dynamic pressure
q1 = 92.7;
m  = W / g; % mass

% Aircraft Stability Characteristics

%% Thrust Coefficients

% Variation of airplane X-axis thrust coefficient with dimensionless speed
C_T_xu = -0.07;

% Steady state thrust coefficient (usually equal to steady state drag)
C_T_x1 = 0.030;

%% Drag Coefficients

% Variation of airplane drag coefficient with dimensionless speed
C_D_u = 0;

% Steady state drag coefficient (Induced drag)
C_D_1 = 0.03;

% Variation of airplane drag coefficient with angle of attack
C_D_a = 0.250;

% Variation of airplane drag coefficient with elevator deflection angle
C_D_de = 0;

%% Lift Force Coefficients

% Variation of airplane lift coefficient with dimensionless speed
C_L_u = 0;

% Variation of airplane lift coefficient with angle of attack
C_L_a = 5.15;

% Steady state lift coefficient
C_L_1 = 0.378;

% Variation of airplane lift coefficient with elevator deflection
C_L_de = 0.5;

% Variation of airplane lift coefficient with dimensionless pitch rate 
C_L_q = 4.1;

%% Pitching Moment Coefficients

% Variation of pitching moment coefficient with dimensionless speed
C_m_u = 0;

% Steady state pitching moment coefficient
C_m_1 = 0;

% Variation of pitching moment coefficient with angle of attack
C_m_a = -0.7;

% Variation of pitching moment coefficient with dimensioness pitch rate
C_m_q = -14.9;

% Variation of pitching moment coefficient with elevator deflection
C_m_de = -1.12;

%% Side Force Coefficients

% Variation of side force coefficient with sideslip angle
C_Y_b = -0.346;

% Variation of side force coefficient with dimensionless rate of change of
% roll rate
C_Y_p = -0.0827;

% Yariation of side force coefficient with dimensionless rate of change of
% yaw rate
C_Y_r = 0.3;

% Variation of side force coefficient with rudder deflection
C_Y_dr = 0.2;

% Variation of side force coefficient with aileron deflection
C_Y_da = 0;

%% Rolling Moment Coefficients

% Variation of rolling moment coefficient with sideslip angle
C_l_b = -0.0944;

% Variation of rolling moment coefficient with rate of change of roll rate
C_l_p = -0.442;

% Variation of rolling moment coefficient with rate of change of yaw rate
C_l_r = 0.0926;

% Variation of rolling moment coefficient with aileron deflection
C_l_da = 0.1810;

% Variation of rolling moment coefficient with rudder deflection
C_l_dr = 0.0150;

%% Yawing Moment Coefficients

% Variation of yawing moment coefficient with sideslip angle
C_n_b = 0.1106;

% Variation of yawing moment coefficient due to thrust with sideslip angle
C_n_T_b = 0;

% Variation of yawing moment coefficient due to rate of change of roll rate
C_n_p = -0.0243;

% Variation of yawing moment coefficient due to rate of change of yaw rate
C_n_r = -0.1390;

% Variation of yawing moment coefficient due to aileron deflection
C_n_da = -0.0254;

% Variation of yawing moment coefficient due to rudder deflection
C_n_dr = -0.0365;

%% Important Inertia Values
Ixx = 7985;
Iyy = 3326;
Izz = 11183;

%% State Space Coefficients: X-Axis Components (Forward-Backward)

% Change in X position due to changes in forward airspeed
Xu  = -(q1 * S_w) / (m * u1) * (C_D_u + 2 * C_D_1);

% Change in X position due to changes in thrust caused by speed changes
XTu = (q1 * S_w) / (m * u1) * (C_T_xu + 2 * C_T_x1);

% Change in X position due to changes in angle of attack
Xa  = (q1 * S_w) / (m) * (-C_D_a + C_L_1);

% Change in X position due to changes in elevator deflection
Xde = (q1 * S_w) / (m) * (C_D_de);

%% State Space Coefficients: Z-Axis Components (Downward-Upward)

% Change in Z position due to changes in forward airspeed
Zu  = -(q1 * S_w) / (m * u1) * (C_L_u + 2 * C_L_1);

% Change in Z position due to changes in angle of attack
Za  = -(q1 * S_w) / (m) * (C_L_a + C_D_1);

% Variation of Z position with pitch rate
Zq  = -(q1 * S_w * mgc) / (2 * m * u1) * (C_L_q);

% Variation of Z position with elevator deflection
Z_de = -(q1 * S_w) / (m) * (C_L_de);

%% State Space Coefficients: M Components (Pitching Moment About Y)

% Variation of pitching moment due to changes in forward airspeed
Mu  = (q1 * S_w * mgc) / (Iyy * u1) * (C_m_u + 2 * C_m_1);

% Variation of pitching moment due to changes in angle of attack
Ma  = (q1 * S_w * mgc) / (Iyy) * (C_m_a);

% Variation of pitching moment due to changes in pitch rate
Mq  = (q1 * S_w * mgc^2) / (2 * Iyy * u1) * (C_m_q);

% Variation of pitching moment due to changes in pitch rate
Mde  = (q1 * S_w * mgc) / (Iyy) * (C_m_de);

%% State Space Coefficients: Y Components (Right-Left)

% Variation of Y position with sideslip angle
Yb = (q1 * S_w) / (m) * (C_Y_b);

% Variation of Y position with roll rate
Yp = (q1 * S_w) / (2 * m * u1) * (C_Y_p);

% Variation of Y position with rudder deflection
Ydr = (q1 * S_w) / (m) * (C_Y_dr);

%% State Space Coefficients: L Components (Rolling Moment About X)

% Variation of rolling moment with sideslip angle
Lb  = (q1 * S_w * b) / (Ixx) * (C_l_b);

% Variation of rolling moment with roll rate
Lp  = (q1 * S_w * b^2) / (2 * Ixx * u1) * (C_l_p);

% Variation of rolling moment with yaw rate
Lr  = (q1 * S_w * b^2) / (2 * Ixx * u1) * (C_l_r);

% Variation of rolling moment with aileron deflection
Lda = (q1 * S_w * b) / (Ixx) * (C_l_da);

% Variation of rolling moment with rudder deflection
Ldr = (q1 * S_w * b) / (Ixx) * (C_l_dr);

%% State Space Coefficients: N Components (Yawing Moment About Z)

% Variation of yawing moment with sideslip angle
Nb  = (q1 * S_w * b) / (Izz) * (C_n_b);

% Variation of yawing moment with roll rate
Np  = (q1 * S_w * b^2) / (2 * Izz * u1) * (C_n_p);

% Variation of yawing moment with yaw rate
Nr  = (q1 * S_w * b^2) / (2 * Izz * u1) * (C_n_r);

% Variation of yawing moment with aileron deflection
Nda = (q1 * S_w * b) / (Izz) * (C_n_da);

% Variation of yawing moment with rudder deflection
Ndr = (q1 * S_w * b) / (Izz) * (C_n_dr);

%% Longitudinal State Space Assembly

% Recursive Variables:
% airspeed (U)
% angle of attack (AOA)
% pitch rate (q)
% and pitch angle (theta)

% Control Input(s): 
% elevator deflection (delta e)


% Plant Matrix Shell
A_Long      = zeros(4,4);

% First Column Inputs (U)
A_Long(1,1) = Xu + XTu; 
A_Long(2,1) = Zu / u1; 
A_Long(3,1) = Mu;
A_Long(4,1) = 0;

% Second Column Inputs (AOA)
A_Long(1,2) = Xa;
A_Long(2,2) = Za / u1;
A_Long(3,2) = Ma;
A_Long(4,2) = 0;

% Third Column Inputs (q)
A_Long(1,3) = 0;
A_Long(2,3) = 1 + (Zq / u1);
A_Long(3,3) = Mq;
A_Long(4,3) = 1;

% Fourth Column Inputs (theta)
A_Long(1,4) = -g;
A_Long(2,4) = 0;
A_Long(3,4) = 0;
A_Long(4,4) = 0;

% Control Matrix Shell
B_Long      = zeros(4,1);

% First Column Inputs (Elevator deflection)
B_Long(1,1)   = Xde; 
B_Long(2,1)   = Z_de / u1; 
B_Long(3,1)   = Mde; 
B_Long(4,1)   = 0;

%% Lateral State Space Assembly

% Recursive Variables:
% Change in sideslip (beta)
% Change in roll rate (p)
% Change in yaw rate (r)
% Change in roll angle (phi)

% Control Input(s): 
% Aileron deflection
% Rudder deflection

% Plant Matrix Shell
A_Late      = zeros(4,4);

% First Column Inputs (Beta)
A_Late(1,1) = Yb/u1;
A_Late(2,1) = Lb; 
A_Late(3,1) = Nb;
A_Late(4,1) = 0;

% Second Column Inputs (p)
A_Late(1,2) = Yp/u1;
A_Late(2,2) = Lp;
A_Late(3,2) = Np;
A_Late(4,2) = 1;

% Third Column Inputs (r)
A_Late(1,3) = -1;
A_Late(2,3) = Lr;
A_Late(3,3) = Nr;
A_Late(4,3) = 0;

% Third Column Inputs (r)
A_Late(1,4) = g/u1;
A_Late(2,4) = 0;
A_Late(3,4) = 0;
A_Late(4,4) = 0;

% Control Matrix Shell
B_Late      = zeros(4,2);

% First Column Inputs (Aileron Deflection)
B_Late(1,1) = 0; 
B_Late(2,1) = Lda; 
B_Late(3,1) = Nda; 
B_Late(4,1) = 0;

% Second Column Inputs (Rudder Deflection)
B_Late(1,2) = Ydr / u1; 
B_Late(2,2) = Ldr; 
B_Late(3,2) = Ndr; 
B_Late(4,2) = 0;

%% Directional Dynamics Analysis

% Longitudinal Dynamics
[V_Long,G_Long] = eig(A_Long);
eig_Long = diag(G_Long); % Gives Poles
tau_Long = -1 * 1 ./ (real(eig_Long)); % Time constants
w_n_Long = abs(eig_Long); % Natural frequencies
zeta_Long = -1 * real(eig_Long) ./ w_n_Long; % Damping ratios

% Short Period Calculations
[w_n_sh, idx_sh]  = max(w_n_Long);
zeta_sh = zeta_Long(idx_sh);
tau_sh  = tau_Long(idx_sh);

% Phugoid Motion Calculations
[w_n_ph, idx_ph]  = min(w_n_Long);
zeta_ph = zeta_Long(idx_ph);
tau_ph  = tau_Long(idx_ph);

% Lateral Dynamics
[V_Late,G_Late] = eig(A_Late);
eig_Late = diag(G_Late); % Gives Poles

tau_Late = -1 ./ (real(eig_Late)); % Time constants
w_n_Late = abs(eig_Late); % Natural frequencies
zeta_Late = -1 * real(eig_Late) ./ w_n_Late; % Damping ratios

% Spiral Motion Calculations
[tau_sp, idx_sp]  = max(tau_Late);
zeta_sp = zeta_Late(idx_sp);
w_n_sp  = w_n_Late(idx_sp);

% Rolling Motion Calculations
[tau_ro, idx_ro]  = min(tau_Late);
zeta_ro = zeta_Late(idx_ro);
w_n_ro  = w_n_Late(idx_ro);

% Dutch Roll Calculations
[zeta_dr, idx_dr]  = min(zeta_Late);
tau_dr = tau_Late(idx_dr);
w_n_dr  = w_n_Late(idx_dr);

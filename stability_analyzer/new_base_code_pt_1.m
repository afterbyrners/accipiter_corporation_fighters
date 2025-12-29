% Conceptual Stability Analysis 

% Description: This software takes the geometric parameters of the various
% flying surfaces of the airplane in order to build up the stability analysis
% coefficients needed to analyze the vehicle as a full assembly.

% It currently only works for normalized tails, meaning that it will not
% work for a plane that uses a V-tail. This is being developed further.

folder = fileparts(which("new_base_code_pt_1.m"));
addpath(genpath(folder));

clear


%% Problem Characteristics

% Given flight parameters
rho = 7.568e-04; % density, sl/ft^3
M = 0.5; % flight mach number, N/A
alpha0L = -3.0; % zero-lift angle of attack, deg
cMAC = -0.01; % zero-lift moment about aerodynamic center, N/A
etaH = 0.95; % dynamic pressure ratio, N/A
etaV = 0.95; % dynamic pressure ratio, N/A
xbarCG = 0.3; % center of gravity in terms of normalized length, N/A
d = 12; % maximum fuselage depth, ft

% Reference Wing Characteristics
ct_wing = 4; % wing tip chord length, ft
cr_wing = 23; % wing root chord length, ft
x_r_t_LE_wing = 28; % chordwise distance from wing tip to wing root LEs, ft
y_r_t_LE_wing = 50; % spanwise distance from wing tip to wing root LEs, ft
b_wing = 113; % wing span, ft
cf_c_wing = 0.2; % aileron chord ratio, N/A
z_r_t_LE_wing = 9; % depthwise distance from wing tip to wing root, ft
phi_wing = atan(z_r_t_LE_wing/y_r_t_LE_wing); % dihedral angle of wing, rad
x_wing = 0; % relative x position of wing 
y_wing = 0; % relative y position of wing
z_wing = 0; % distance of wing root to fuselage center line

% Horizontal Stabilizer Characteristics
ct_hori = 4; % horizontal stabilizer tip chord, ft
cr_hori = 12; % horizontal stabilizer root chord, ft
x_r_t_LE_hori = 13; % chordwise distance from h. tip to h. root LEs, ft
y_r_t_LE_hori = 23; % spanwise distance from h. tip to h. root LEs, ft
b_hori = 46; % horizontal stabilizer span, ft
cf_c_hori = 0.3; % elevator chord ratio, from tip station
x_hori = 55; % long. distance from wing to horizontal tail, ft
y_hori = 0; % 
z_hori = 6; % lat. distance from wing to horizontal tail, ft

% Vertical Stabilizer Characteristics
ct_vert = 6; % vertical stabilizer tip chord, ft
cr_vert = 26; % vertical stabilizer root chord, ft
x_r_t_LE_vert = 27; % chordwise distance from v. tip to v. root LEs, ft
b_vert = 2 * 31; % vertical tail span, ft (multipled by 2 because we assume its the measurement of a full wing)
cf_c_vert = 0.3; % vertical tail flap chord tip length, ft
x_vert = 39; % long. distance from wing to vertical tail, ft
y_vert = 0;
z_vert = 6; % lat. distance from wing to vertical tail, ft

%% Wing Calculations

% Calculate Wing Leading Edge Sweep Angle
sweep_wing = atan(x_r_t_LE_wing/y_r_t_LE_wing);

% Calculate Wing taper ratio, reference area, aspect ratio, MGC + location
[taper_wing, s_wing, AR_wing, mgc_wing, x_mgc_wing, y_mgc_wing] = AppendixA(b_wing,cr_wing,ct_wing,sweep_wing);

% Calculate Wing AC location
xacWing = 0.25; % always 0.25 for subsonic flight
yacWing = 0; % symmetric, no change
zacWing = 0; % assume so that wing is on fuselage centerline

% Calculate Polhamus Slope for the wing
cLAlphaWing = AppendixB(AR_wing,sweep_wing,taper_wing,M); % 1/rad

% Calculate Flap Effectiveness for the wing.
tauWing = AppendixD(cf_c_wing); % assumed? lowkey dont remember

%% Horizontal Stabilizer Calculations

% Calculate hori leading edge sweep angle
sweep_hori = atan(x_r_t_LE_hori / y_r_t_LE_hori); % deg

% Calculate hori taper ratio, reference area, aspect ratio, MGC + location
[taper_hori, s_hori, AR_hori, mgc_hori, x_mgc_hori, y_mgc_hori] = AppendixA(b_hori,cr_hori,ct_hori,sweep_hori);

% Calculate hori ac coordinates
xacHori = (-x_mgc_wing + x_hori + x_mgc_hori + 0.25 * mgc_hori)/mgc_wing; % ft
yacHori = 0; % symmetry, no lateral AC movement
zacHori = z_hori;

% Calculate r & m values for hor. downwash gradient
x_wh = x_hori - (cr_wing/4) + (cr_hori/4);
rHori = x_wh/(0.5*b_wing);
mHori = z_hori/(0.5*b_wing);

% Calculate polhamus slope for the hori.
cL_A_hori = AppendixB(AR_hori,sweep_hori,taper_hori,M); % 1/rad

% Calculate downwash gradient for the hori.
dw_hori = AppendixC(AR_wing,sweep_wing,taper_wing,M,rHori,mHori); % N/A
% Notice that the flying surface properties used here are for the wing.
% This is because the wing produces the downwash that the h. stab.
% experiences.

% Calculate flap effectiveness for the hori.
tau_hori = AppendixD(cf_c_hori);

%% Vertical Stabilizer Calculations

% Calculate Vert. Leading Edge Sweep Angle
sweep_vert = atan(x_r_t_LE_vert / (b_vert/2));  

% Calculate vert taper ratio, reference area, aspect ratio, MGC + location
[taper_vert, s_vert, AR_vert, mgc_vert, x_mgc_vert, z_mgc_vert] = AppendixA(b_vert,cr_vert,ct_vert,sweep_vert);

% Correct area
s_vert = s_vert / 2;

% Calculate Vert. AC Coordinates
xacVert = -x_mgc_wing + x_vert + x_mgc_vert + 0.25 * mgc_vert;  % always this value
zacVert = z_vert + z_mgc_vert;

% Calculate r & m values for Vert. Downwash Gradient
rVert = x_vert / (0.5 * b_wing);
mVert = z_vert / (0.5 * b_wing);

% Calculate Polhamus slope for the Vert.
cLAlphaVert = AppendixB(AR_vert, sweep_vert, taper_vert, M);  % 1/rad

% Calculate Flap Effectiveness for the Vert.
tauVert = AppendixD(cf_c_vert);

%% Calculate Important Parameters

% Calculate Level Flight Coefficient of Lift
cL0 = -(cLAlphaWing * deg2rad(alpha0L));  % N/A

% Calculate Coefficient of Lift caused by General Plane
cLAlphaTotal = cLAlphaWing + etaH * (s_hori / s_wing) * ...
               cL_A_hori * (1 - dw_hori);  % 1/rad

% Calculate Coefficient of Lift caused by Horizontal Stabilizer Incidence
cLih = etaH * (s_hori / s_wing) * cL_A_hori;  % 1/rad

% Calculate Coefficient of Lift caused by Elevator Flaps
cLde = cLih * tau_hori;  % 1/rad

% Calculate Level Flight Coefficient of Moment
cM0 = cMAC + cL0 * (xbarCG - xacWing);  % N/A

% Calculate Coefficient of Moment caused by General Plane
cMalphaTotal = cLAlphaWing * (xbarCG - xacWing) - ...
               cL_A_hori * etaH * (s_hori / s_wing) * ...
               (1 - dw_hori) * (xacHori - xbarCG);

% Calculate Coefficient of Moment caused by Horizontal Stabilizer Incidence
cMih = -cL_A_hori * etaH * (s_hori / s_wing) * (xacHori - xbarCG);

% Calculate Coefficient of Moment caused by Flap Angles
cMde = -cL_A_hori * etaH * (s_hori / s_wing) * ...
        tau_hori * (xacHori - xbarCG);

% Calculate Neutral Point
xbarNP = (xacWing + (cL_A_hori / cLAlphaWing) * ...
          etaH * (s_hori / s_wing) * ...
          (1 - dw_hori) * xacHori) / ...
         (1 + (cL_A_hori / cLAlphaWing) * ...
          etaH * (s_hori / s_wing) * (1 - dw_hori));

xNP = xbarNP * mgc_wing;

% Calculate Static Margin Percentage
SM = (xbarNP - xbarCG) * 100;

% Calculate Level Flight Coefficient of Side Force
cY0 = 0;  % Given/assumed

% Calculate Coefficient of Side Force based on Beta
cYBv = -(s_vert / s_wing) * cLAlphaVert * AppendixE(s_vert, s_wing, z_wing, d, AR_wing, sweep_wing, taper_wing);

cYBw = -0.0001 * abs(phi_wing) * 180 / pi;  % rad^-1, basically zero
cYBwv = cYBv +cYBw;
cYBf = 0.3 * cYBw;

cYB = cYBwv + cYBf;

% Calculate change in Coefficient of Side Force caused by aileron incidence
cYda = 0;  % we can assume so due to ailerons not affecting that DOF

% Calculate change in Coefficient of Side Force caused by rudder incidence
cYdr = etaV * (s_vert / s_wing) * cLAlphaVert * tauVert;

%% PROBLEM B

% Calculate Level Flight Coefficient of Rolling Moment
cl0 = 0;  % Assumed

% Calculate Coefficient of Rolling Moment change based on beta
clBv = cYBv * (z_mgc_vert / b_wing);
clBw = -2 * cLAlphaWing * phi_wing * (y_mgc_wing / b_wing);
clB  = clBv + clBw;

% Calculate Coefficient of Rolling Moment change caused by ailerons
aileron_end = 51;
aileron_start = 35;

clda = AppendixG(cLAlphaWing,tauWing,cr_wing,s_wing,b_wing,taper_wing,aileron_start,aileron_end);

% Calculate Coefficient of Rolling Moment change caused by rudder
cldr = cYdr * (z_mgc_vert / b_wing);

%% PROBLEM C

% Calculate Level Flight Coefficient of Yawing Moment
cn0 = 0;  % assumed

% Calculate Coefficient of Yawing Moment change caused by beta
cnB = -cYBv * (xacVert - xbarCG * mgc_wing) / b_wing;

% Calculate Coefficient of Yawing Moment change caused by ailerons
cnda = 0;  % assumed so (low dihedral = min effect)

% Calculate Coefficient of Yawing Moment change caused by rudder
cndr = -cYdr * (xacVert - xbarCG * mgc_wing) / b_wing;

% Calculate Coefficient of Yawing Moment caused by change in roll rate
cndpv =  
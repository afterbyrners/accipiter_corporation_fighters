% Raymer Subsonic Lift Curve

% Based on Raymer's conceptual design book, this function is a
% semi-empirical formula that calculates the wing lift-curve slope per
% radian, with validity up to the drag-divergent mach number and up to a
% Mach of unity (1) for a swept wing.

% It is worth noting that some estimations and assumptions were made.

% Since fighters are usually lifting bodies, an assumption of 0.99
% (Sexp/Sref) * Fref is made. 

% This is for the 64-206 airfoil, which has an a0 of -1.0 and a section
% lift derivative of 0.110 per degree, and has max thickness at 0.4 of c.

%% INPUTS & OUTPUTS

% AR       = aspect ratio of the wing
% LE_sweep = leading edge sweep angle, in degrees
% lambda   = taper ratio of the wing
% M        = flight mach number

% CL_a     = coefficient of lift derivative to AoA for a lifting surface

function CL_a = RaymerSubLift(AR,sweep_LE,lambda,M) % CL_alpha function

   % Initialize, this should change later
   sweep_LE = deg2rad(sweep_LE);
   Seff = 0.99; 
   eta = 0.95;
   sweep_t_max = atan( tan(sweep_LE) - (4 * 0.4 * (1-lambda) ) / (AR * (1+lambda) ) );

   % Korn's relation for drag divergence mach number
   K = 0.95;
   tc = 0.06;
   Cl = 0.277478; 
   % (based on CL found from wing loading versus dynamic pressure at trim)
   Mdd = K - tc - (0.1 * Cl) ;
   
   % Subsonic Curve Section
   if M < Mdd
   B = sqrt(1 - M^2);
   CL_a = (2 * pi * AR * Seff) / (2 + sqrt(4 + ((AR^2 * B^2) / eta^2) * (1 + (tan(sweep_t_max)^2 / B^2)))) ;

   % Supersonic Curve Section
   elseif M > (1/cos(sweep_LE))
   B = sqrt(M^2 - 1);
   CL_a = 4/B * (1 - (1 / (2 * AR * B)));

   % Transonic Curve Section
   else 
   


   % Subsonic endpoint
   M_1 = Mdd;
   B_1 = sqrt(1 - M_1^2);
   CL_a_sub = (2*pi*AR*Seff) / (2 + sqrt(4 + ((AR^2*B_1^2)/eta^2) * (1 + (tan(sweep_t_max)^2 / B_1^2))));

   % Supersonic endpoint
   M_2 = 1/cos(sweep_LE);
   B_2 = sqrt(M_2^2 - 1);
   CL_a_sup = 4/B_2 * (1 - (1 / (2 * AR * B_2)));

   % Slender Wing M = 1 approximation
   M_p = 1; 
   CL_a_p = pi * AR / 2;


   % Interpolate
   M_fair  = [M_1, M_p, M_2];
   CL_fair = [CL_a_sub, CL_a_p, CL_a_sup];
   CL_a = interp1(M_fair, CL_fair, M, 'makima');
   
   end

end

% B term flips when in supersonic

% Stored artifacts
   % Cl_a = 0.110;
   % Cl_a = Cl_a / deg2rad(1);
   % eta = Cl_a / (2 * pi / B); May not be necessary

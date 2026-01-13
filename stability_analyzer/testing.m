% Literally just buidling new functiosdnhadn

% Currently testing Raymer Wing Lift with the main basis of Roskam added as
% a side bonus and quest lol 

folder = fileparts(which("testing.m"));
addpath(genpath(folder));


clear

M = linspace(0,3,300);

AR = 5.5851;
sweep_LE = 3.63529;
lambda = 0.55;

r_h = 1.3560;
m_h = -0.0329;

CL_a1 = zeros(size(M));
dw_h = zeros(size(M));

for i = 1:length(AR)
    for j = 1:length(M)
        [CL_a1(i,j),Mdd,Mss] = RaymerWingLift(AR(i), sweep_LE, lambda, M(j));
        %dw_h(i,j) = AppendixC(AR(i),deg2rad(40),lambda,M(j),r_h,m_h);
    end
end

for i = 1:length(AR)
    for j = 1:length(M)
        [CL_a2(i,j)] = AppendixB(AR(i), sweep_LE, lambda, M(j));
     end
end

hold on
plot(M,CL_a1)
plot(M,CL_a2)
xline(0.85)
xlabel('Mach')
ylabel('Lift Curve')
ylim([0,10])
xlim([0,3])


% Conclusion: I believe that the Raymer one is much more accurate. 
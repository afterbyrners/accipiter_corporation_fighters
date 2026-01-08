% Literally just buidling new functiosdnhadn
folder = fileparts(which("testing.m"));
addpath(genpath(folder));


clear

M = linspace(0,3,300);

AR = [1,1.5,2,2.5,3,4];
sweep_LE = 40;
lambda = 0.08;

CL_a = zeros(size(M));
for i = 1:length(AR)
    for j = 1:length(M)
        CL_a(i,j) = RaymerSubLift(AR(i), sweep_LE, lambda, M(j));
    end
end
plot(M,CL_a)
xlabel('Mach')
ylabel('Lift Curve')
ylim([0,10])
xlim([0,3])
clear
clc

tic

d2r = pi/180.0;
mu = 3.986e5; % km^3/s^2, Earth's gravitational parameter
secondsPerDay = 86400;

% Inputs are from the 30 day propagation of the orbit from TLE
r0 = [ 1895874.38838166,  4497007.59724938, -4507886.60518927];
v0 = [-7120.17746759,   3051.88529362,    73.34216649];
r0 = r0/1e3;
v0 = v0/1e3;

oe0 = [6645.82965, 2.06999581e-3, d2r*42.7320639, ...
       d2r*156.211622, d2r*8.19853023, d2r*80.7678293];
oef = [6657.391879, 0.002595, d2r*42.748082, ...
       d2r*345.325883, d2r*124.412552, d2r*287.718564];

n = sqrt(mu/(oe0(1)^3));
T = 2*pi/n;

% First task is to generate the disk of positions and velocities
% associated with the original and final orbits

nStepsR = 40;
nStepsT = 500;

f1Start = 3;
f1End = 3.5;
f1Step = (f1End - f1Start)/nStepsR;

f2Start = 0.3;
f2End = 0.8;
f2Step = (f2End - f2Start)/nStepsR;

tofs = linspace(7000, 9000, nStepsT);
r0 = zeros(nStepsR, 3);
v0 = zeros(nStepsR, 3);
rf = zeros(nStepsR, 3);
vf = zeros(nStepsR, 3);

oefn = oe0;
oefn(4) = oef(4);
%oef = oefn;

for i=1:nStepsR
    f1 = i*f1Step + f1Start;
    f2 = i*f2Step + f2Start;
    oe0(6) = f1;
    oef(6) = f2;
    [r0(i, :), v0(i, :)] = coe2rv(oe0, mu, 'rad');
    [rf(i, :), vf(i, :)] = coe2rv(oef, mu, 'rad');
end

% Initializing our brute force search
costs = zeros(nStepsR, nStepsR, nStepsT);
dv1 = zeros(nStepsR, nStepsR, nStepsT);
dv2 = zeros(nStepsR, nStepsR, nStepsT);

for i=1:nStepsR
    for j=1:nStepsR
        for k=1:nStepsT
            % Choose our elements
            r = r0(i, :);
            v = v0(i, :);
            
            re = rf(j,:);
            ve = vf(j,:);
            
            tof = tofs(k)/secondsPerDay;
            m = floor(tof/T) + 1;
            [V1, V2, extremal_distances, exitflag] = lamberti(r, re, tof, m, ...
                                                  mu);
            
            if isnan(V2)
                errV = 10000.0;
                dv2(i, j, k) = 10000.0;
            else
                errV = norm(V2 - ve');
                dv2(i, j, k) = norm(V2 - ve');
            end
            
            if isnan(V1)
                errV = 10000.0;
                dv1(i, j, k) = 10000.0;
            else
                errV = errV + norm(V1 - v');
                dv1(i, j, k) = norm(V1 - v');
            end
            
            costs(i, j, k) = errV;
        end
    end
end

%min(costs, [], 'all')
[m, index] = min3d(costs);
burn1Pos = r0(index(1), :);
burn2Pos = rf(index(2), :);
f1burn = index(1)*f1Step + f1Start;
f2burn = index(2)*f2Step + f2Start;
tof = tofs(index(3));

fprintf('Minimum Burn Cost: %6.4f\n', m);
fprintf('Position for Burn One: %6.4f %6.4f %6.4f\n', burn1Pos(1), ...
        burn1Pos(2), burn1Pos(3));
fprintf('Delta Velocity, Burn One: %6.4f\n', dv1(index(1), index(2), ...
                                                 index(3)));
fprintf('Position for Burn Two: %6.4f %6.4f %6.4f\n', burn2Pos(1), ...
        burn2Pos(2), burn2Pos(3));
fprintf('Delta Velocity, Burn Two: %6.4f\n', dv2(index(1), index(2), ...
                                                 index(3)));
fprintf('Time of Flight: %6.4f\n', tof);
fprintf('True Anomaly of First Burn: %6.4f\n', f1burn);
fprintf('True Anomaly of Second Burn: %6.4f\n', f2burn);


toc

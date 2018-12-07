clear
clc

tic

d2r = pi/180.0;
mu = 3.986e5; % km^3/s^2, Earth's gravitational parameter
secondsPerDay = 86400;

% Inputs are from the 30 day propagation of the orbit from TLE
r0 = [6635837.6238441  -434376.78212234 -183676.93413547];
v0 = [ 514.4119716  5660.01991905 5250.16257648];
r0 = r0/1e3;
v0 = v0/1e3;

oe0 = [6646.80789, 8.86585360e-4, d2r*42.7745499, ...
       d2r*2.03443527, d2r*0, d2r*2.32385299];
oef = [6657.391879, 0.002595, d2r*42.748082, ...
       d2r*345.325883, d2r*124.412552, d2r*287.718564];

n = sqrt(mu/(oe0(1)^3));
T = 2*pi/n;

% First task is to generate the disk of positions and velocities
% associated with the original and final orbits

nStepsR = 50;
nStepsT = 50;

f1Start = 1;
f1End = 2;
f1Step = (f1End - f1Start)/nStepsR;

f2Start = 1.5;
f2End = 2.5;
f2Step = (f2End - f2Start)/nStepsR;

tofs = linspace(7620.0, 7640, nStepsT);
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
ms = zeros(nStepsR, nStepsR, nStepsT);
v1s = zeros(nStepsR, nStepsR, nStepsT, 3);
v2s = zeros(nStepsR, nStepsR, nStepsT, 3);

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
                errV = 10.0;
                dv2(i, j, k, :) = 10.0;
            else
                errVec = V2 - ve;
                errV = norm(errVec);
                dv2(i, j, k, :) = norm(errVec);
            end
            
            if isnan(V1)
                errV = 10.0;
                dv1(i, j, k, :) = 10.0;
            else
                errVec = V1 - v;
                errV = errV + norm(errVec);
                dv1(i, j, k, :) = norm(errVec);
            end
            
            costs(i, j, k) = errV;
            ms(i, j, k) = m;
            v1s(i, j, k, :) = V1;
            v2s(i, j, k, :) = V2;
        end
    end
end

%min(costs, [], 'all')
fileID = fopen('results_maneuvers.txt','w');
[minCost, index] = min3d(costs);
burn1Pos = r0(index(1), :);
burn2Pos = rf(index(2), :);
f1burn = index(1)*f1Step + f1Start;
f2burn = index(2)*f2Step + f2Start;
tof = tofs(index(3));

minCost = dv1(index(1), index(2), index(3)) + dv2(index(1), index(2), ...
                                                  index(3));

ddv1 = squeeze(v1s(index(1), index(2), index(3), :)) - v0(index(1), :)';
ddv2 = squeeze(v2s(index(1), index(2), index(3), :)) - vf(index(2), :)';

fprintf(fileID,'Minimum Burn Cost: %6.4f km/s\n', minCost);
fprintf(fileID,'Position for Burn One: %6.4f %6.4f %6.4f km, ECI\n', burn1Pos(1), ...
        burn1Pos(2), burn1Pos(3));
fprintf(fileID,'Delta Velocity, Burn One: %6.4f km/s\n', dv1(index(1), index(2), ...
                                                 index(3)));
fprintf(fileID,'Position for Burn Two: %6.4f %6.4f %6.4f km, ECI\n', burn2Pos(1), ...
        burn2Pos(2), burn2Pos(3));
fprintf(fileID,'Delta Velocity, Burn Two: %6.4f km/s\n', dv2(index(1), index(2), ...
                                                 index(3)));
fprintf(fileID,'Time of Flight: %6.4f seconds\n', tof);
fprintf(fileID,'True Anomaly of First Burn: %6.4f rad\n', f1burn);
fprintf(fileID,'True Anomaly of Second Burn: %6.4f rad\n', f2burn);

fprintf(fileID,'Burn 1: %6.5f %6.5f %6.5f km/s, ECI\n', ddv1(1), ddv1(2), ddv1(3));
fprintf(fileID,'Burn 2: %6.5f %6.5f %6.5f km/s, ECI\n', ddv2(1), ddv2(2), ddv2(3));


toc

% Plot our results to show that we have, qualitatively, reached an
% optimum.

total_delta = dv1+dv2;

figure;
plot(tofs(:), squeeze(total_delta(index(1), index(2), :)));
hold on;
xlabel('Time of Flight');
ylabel('Total Maneuver Impulse (m/s)');
title('Cost Function of Maneuver');
saveas(gcf, 'ToF_Cost','fig');
% close(gcf);

figure;
plot(linspace(f1Start, f1End, nStepsR), squeeze(total_delta(:, index(2), index(3))));
hold on;
xlabel('True Anomaly of Startpoint (rad)');
ylabel('Total Maneuver Impulse (m/s)');
title('Cost Function of Maneuver');
saveas(gcf, 'Start_Cost','fig');
% close(gcf);

figure;
plot(linspace(f2Start, f2End, nStepsR), squeeze(total_delta(index(1), :, index(3))));
hold on;
xlabel('True Anomaly of Endpoint (rad)');
ylabel('Total Maneuver Impulse (m/s)');
title('Cost Function of Maneuver');
saveas(gcf, 'End_Cost','fig');
% close(gcf);

clear
clc

d2r = pi/180.0;
mu = 3.986e5; % km^3/s^2, Earth's gravitational parameter
secondsPerDay = 86400;

% Inputs are from the 30 day propagation of the orbit from TLE
r0 = [ 1895874.38838166,  4497007.59724938, -4507886.60518927];
v0 = [-7120.17746759,   3051.88529362,    73.34216649];
oe0 = [6645.82965, 2.06999581e-3, d2r*42.7320639, ...
       d2r*156.211622, d2r*8.19853023, d2r*80.7678293];

r0 = r0/1e3;
v0 = v0/1e3;
% km, NA, deg, deg, deg, deg
oef_orig = [6657.391879, 0.002595, d2r*42.748082, ...
            d2r*345.325883, d2r*124.412552, d2r*287.718564];

oef = oe0;
oef(4) = oef_orig(4);

n = sqrt(mu/(oe0(1)^3));
T = 2*pi/n;

Mf = E2M(f2E(oef(6), oef(2)), oef(2));
M0 = E2M(f2E(oe0(6), oe0(2)), oe0(2));

hk0 = -mu/(2*oe0(1));
hkf = -mu/(2*oef(1));

fprintf('Energy difference between orbits: %7.6f\n', hk0 - hkf);

% Try setting the true anomaly of the original position closer to
% the final
tDiff = 1; %s
Mnew = Mf - tDiff*n;
fnew = E2f(M2E(Mnew, oe0(2)), oe0(2));
oe0(6) = fnew;

[r, v] = coe2rv(oef, mu, 'rad');
[r0, v0] = coe2rv(oe0, mu, 'rad');

m = 0;

tof = tDiff;
tof = tof/secondsPerDay;
T = T/secondsPerDay;

[V1, V2, extremal_distances, exitflag] = lamberti(r0', r', tof, m, ...
                                                  mu);
minErr = 1000000;
tof_Final = tof;
errV = norm(V2 - v');

for i=1:1000000
    tof = tof + (1/secondsPerDay);
    m = floor(tof/T) + 1;
    [V1, V2, extremal_distances, exitflag] = lamberti(r0', r', tof, m, ...
                                                 mu);
    if isnan(V2)
        errV = 100.0;
    else
        errV = norm(V2 - v');
    end
    
    if errV < minErr
        minErr = errV;
        V1_Final = V1;
        V2_Final = V2;
        tof_Final = tof;
    end
    
    if m > 500
        break;
    end
end

oef_lambert = rv2coe(r', V2_Final, mu, 'rad');

fprintf('Error in the inclination after burn: %4.4f\n', oef_lambert(3) ...
        - oef_orig(3));

fprintf('The final orbital velocity error is %7.6f\n', minErr);
fprintf('The time of flight is %3.2f days \n', tof_Final);
fprintf('The total cost of the maneuver is %7.6f\n', norm(V1_Final ...
                                                  - v0) + norm(V2_Final ...
                                                  - v));
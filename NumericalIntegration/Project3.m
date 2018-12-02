clear
clc

d2r = pi/180.0;
mu = 3.986e5; % km^3/s^2, Earth's gravitational parameter

% Inputs are from the 30 day propagation of the orbit from TLE
r0 = [ 1895874.38838166,  4497007.59724938, -4507886.60518927];
v0 = [-7120.17746759,   3051.88529362,    73.34216649];
oe0 = [6645.82965, 2.06999581e-3, d2r*42.7320639, ...
       d2r*156.211622, d2r*8.19853023, d2r*80.7678293];

r0 = r0/1e3;
v0 = v0/1e3;
% km, NA, deg, deg, deg, deg
oef = [6657.391879, 0.002595, d2r*42.748082, ...
       d2r*345.325883, d2r*124.412552, d2r*287.718564];

Mf = E2M(f2E(oef(6), oef(2)), oef(2));
M0 = E2M(f2E(oe0(6), oe0(2)), oe0(2));

oef(6);
E2f(M2E(Mf, oef(2)), oef(2));
err = oef(6) - E2f(M2E(Mf, oef(2)), oef(2)) - 2*pi;

hk0 = -mu/(2*oe0(1));
hkf = -mu/(2*oef(1));

fprintf('Energy difference between orbits: %7.6f\n', hk0 - hkf);

if Mf < 0
    Mf = Mf + 2*pi;
end

n = sqrt(mu/(oe0(1)^3));

T = 2*pi/n;

[r, v] = coe2rv(oef, mu, 'rad');
% Try setting the true anomaly of the original position closer to
% the final
tDiff = 1800; %s
Mnew = Mf - tDiff*n;
fnew = E2f(M2E(Mnew, oe0(2)), oe0(2));
oe0(6) = fnew;
[r0, v0] = coe2rv(oe0, mu, 'rad');

m = 0;

tof = (Mf - Mnew)/n + m*T

[V1, V2, extremal_distances, exitflag] = lamberti(r0', r', tof, m, ...
                                                  mu);
errV = norm(V2 - v')

%while (errV >= 1)
    %tof = tof + 0.1;
    %m = floor(tof/T);
    %[V1, V2, extremal_distances, exitflag] = lamberti(r0, r', tof, m, ...
    %                                             mu);

    %errV = norm(V2 - v');
    
    %if m > 5
        %break;
        %end
    %end                                     
    
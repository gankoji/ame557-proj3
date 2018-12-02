function E = f2E(f, e)

% Given true anomaly in radians, return eccentric anomaly in
% radians. e is eccentricity

b = sqrt((1 - e)/(1 + e));

tanHalf = b*tan(f/2);

E = 2*atan(tanHalf);


function M = E2M(E, e)

% Given eccentric anomaly in radians, return mean anomaly in
% radians. e is eccentricity

M = E - e*sin(E);
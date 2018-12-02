function [ coe ] = rv2coe( r, v, mu, angletype )
%{
FILENAME: rv2coe

DESCRIPTION: This function computes the classical (Keplerian) orbit
elements given the Cartesian position and velocity vectors, the
gravitational parameter of the attracting body, and the desired angletype 

INPUTS:
    r           = Cartesian position vector [ km ]
    v           = Cartesian velocity vector [ km/s ]
    mu          = gravitational parameter of attracting body [ km^3/s^2 ]
    angletype   = string indicating desired angletype for inputs
                    'deg' = degrees, 'rad' = radians

OUTPUTS:
    coe         = vector which contains the classical orbit elements
                    [ a; e; i; Om; w; f ]

FUNCTIONS/SCRIPTS CALLED:
    none
%}

%  Create empty vector for keplerian elements
coe = zeros(6, 1);

%  Conversions
r2d = 180/pi; 

%  Coordinate frame unit vectors
I = [1 0 0]; J = [0 1 0]; K = [0 0 1];

rhat = r/norm(r);                       % position unit vector [km]     
h    = cross(r,v);                      % angular momentum {h = r x v}
hhat = h/norm(h);                       % normalized angular momentum

nhat = cross(K,h)/norm(cross(K,h));     % normalized ascending node vector

%  Eccentricity
e      = (1/mu)*cross(v,h) - rhat;      % eccentricity vector
coe(2) = norm(e);

energy = (1/2)*dot(v,v) - mu/norm(r);   % energy, km^2/s^2
%  If energy < 0, the orbit is closed (periodic)

%  Semi-major axis (a) and parameter (p)
if coe(2) ~= 1
    coe(1) = -mu/(2*energy);            % {energy = -mu/(2*a)}
else
    coe(1) = inf;
end
   
%  Inclination (i) of orbit
coe(3) = acos(dot(K,hhat));     
%  If i < 90 deg, the elliptical orbit is a direct (prograde) orbit

%  Right ascension of the ascending node (Omega)
coe(4) = mod(atan2(dot(J,nhat),dot(I,nhat)), 2*pi);

%  Argument of periapsis (w)
coe(5) = mod(atan2(dot(hhat,cross(nhat,e)),dot(nhat,e)), 2*pi); 

%  True anomaly (f) at epoch [rad]
coe(6) = mod(atan2(dot(hhat,cross(e,r)),dot(e,r)), 2*pi);

if(strcmp(angletype, 'deg'))
    coe(3:end) = coe(3:end)*r2d;
end

end % ---- end function
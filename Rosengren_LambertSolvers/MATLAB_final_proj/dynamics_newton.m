function [ dy ] = dynamics_newton( ~, y, const )

%  Inertial position and velocity vectors
r = y(1:3);
v = y(4:6);

rmag  = norm(r);                    % orbital radius [ km ]
rhat  = r/rmag;                     % radius unit vector
znorm = rhat'*[0 0 1]';             % z-comp of r normalized by radius

% ------------------------------------------------------------------------
%  CENTRAL BODY ACCELERATION
% ------------------------------------------------------------------------

%  Gravitational acceleration [ km/s^2 ] due to a point-mass 
a_kep = -const.mu/rmag^3*r;

%  Gravitational acceleration [ km/s^2 ] due to Earth's oblateness 
a_J2 = -3/2*const.mu*const.J2*const.Re^2/rmag^4*((1 - 5*znorm^2)*rhat + 2*znorm*[0 0 1]');
     
% ------------------------------------------------------------------------
%  ATMOSPHERIC DRAG ACCELERATION 
% ------------------------------------------------------------------------

%  Ballistic coefficient [ kg/m^2 ]
Bstar = 0.13071e-3;
BC    = 1/12.741621/Bstar; 

%  Atmospheric density [ kg/m^3 ]
rho = atmosphere( norm(r) - const.Re );

%  Gravitational acceleration [ km/s^2 ] due to atmospheric drag
a_drag = -0.5/BC*rho*norm(v)*v*(1e3);

% -------------------- %
%  STATE DIFFERENTIAL  %
% -------------------- %

dy = [ v ; a_kep + a_J2 + a_drag ];

end % ----- End Function
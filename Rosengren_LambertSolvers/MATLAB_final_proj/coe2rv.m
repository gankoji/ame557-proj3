function [ r, v ] = coe2rv( coe, mu, angletype )
%{
FILENAME:  coe2rv

DESCRIPTION: This function computes the cartesian position and velocity
vectors given the classical elements, gravitational parameter of the 
attracting body, and angletype used for the inputs.

INPUTS:
    coe         = vector which contains the classical orbit elements
                    [ a; e; i; Om; w; f ]
    mu          = gravitational parameter of attracting body [ km^3/s^2 ]
    angletype   = string indicating desired angletype for inputs
                    'deg' = degrees, 'rad' = radians

OUTPUTS:
    r           = Cartesian position vector [ km ]
    v           = Cartesian velocity vector [ km/s ]

FUNCTIONS/SCRIPTS CALLED:
    none
%}

%  Conversions
d2r = pi/180;
    
%  Classical orbit elements
a  = coe(1);    
ec = coe(2);   
in = coe(3);   
Om = coe(4);   
w  = coe(5);  
f  = coe(6);   

if(strcmp(angletype, 'deg'))
    in = in*d2r;
    Om = Om*d2r;
    w  = w*d2r;
    f  = f*d2r;
end    

%  Argument of latitude [ rad ]
th = w + f;

%  Magnitude of position vector [ km ]
r = a*(1 - ec^2)/(1 + ec*cos(f));       % trajectory equation

%  Magnitude of velocity vector [ km/s ]
v = sqrt(mu*(2/r - 1/a));               % vis-viva equation

%  Vector kinematic quantities 
nhat = [ cos(Om) ; sin(Om) ; 0 ];   
rT   = [ -cos(in)*sin(Om) ; cos(in)*cos(Om) ; sin(in) ];

%  Flight-path angle [ rad ]
gamma = atan2(ec*sin(f), 1 + ec*cos(f));          

%  Normalized position and velocity vectors
rhat = cos(th)*nhat + sin(th)*rT;                   
vhat = sin(gamma - th)*nhat + cos(gamma - th)*rT;  

%  Position [ km ] and velocity [ km/s ]  vectors
r = r*rhat;                                         
v = v*vhat;                                          

end % ----- End Function
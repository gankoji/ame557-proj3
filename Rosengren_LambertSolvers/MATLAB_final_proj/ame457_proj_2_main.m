clearvars; close all; clc

tic

% ------------------------------------------------------------------------
%     CONSTANTS & PARAMETERS | EARTH, SATELLITE     
% ------------------------------------------------------------------------

%  Planetary constants
const.mu = 3.98600442e5;            % Earth gravitational parameter [ km^3/s^2 ]
const.Re = 6378.137;                % Earth radius [ km ]
const.J2 = 0.0010826267;            % Earth oblateness 

%  Conversions
jday = 86400;                       % seconds in a day

%  Initial classical orbit elements [ a; e; i; Om; w; f ]
coe0 = [ 6657.3918791 ;
         0.0025945 ;
         42.7480815 ;
         345.3258830 ;
         124.4125522 ;
         287.3948092 ];

%  Conversion cartesian position and velocity 
[ r0, v0 ] = coe2rv( coe0, const.mu, 'deg' );

% ------------------------------------------------------------------------
%     ORBIT PROPAGATION | NEWTONIAN (ODE45)     
% ------------------------------------------------------------------------

%  Initial state and propagation timespan
state0 = [ r0; v0 ];
tvec   = linspace(0, 175*jday, 1000);

%  Set integrator options
tol = 1e-9;
options = odeset( 'RelTol', tol, 'AbsTol', tol, 'Events', @reentryEvents );

%  State propagation                 
[time, statet] = ode45( @dynamics_newton, tvec, state0, options, const );      

%  Trim last entry (reentry condition is inside Earth)
idx = length(time) - 1;

%  Extrate position and velocity vectors from state
r = statet(1:idx,1:3)';   
v = statet(1:idx,4:6)';

%  Convert to the classical elements
coe = zeros(6,idx);
for k = 1:idx
   
    coe(:,k) = rv2coe( r(:,k), v(:,k), const.mu, 'deg' );  
    
end

% ------------------------------------------------------------------------ 
%     FIGURES     
% ------------------------------------------------------------------------
tdays = time(1:idx)/jday;

%  Reentry time [ days ]
disp(' ')
str1 = 'The total lifetime of the satellite';
fprintf('%s is %4.4f days.\n', str1, tdays(end))
disp(' ')

hFig = figure(1); 
    hold all
    plot( tdays , coe(1,:) , 'LineWidth' , 1.2 );
    hXLabel = xlabel( [ '$\textrm{Time [days]}', '$' ] , ...
                        'Interpreter', 'LaTex' );
    hYLabel = ylabel( [ '$a\, \textrm{ [km]}', '$' ] , ...
                        'Interpreter', 'LaTex' , 'Rotation', 90 );      
    set( gca , 'FontName' , 'Helvetica' );
    set( hXLabel , 'FontSize' , 12 )
    set( hYLabel , 'FontSize' , 14 )
    set( hFig , 'units' , 'inches' , ... 
         'NumberTitle' , 'off' , 'Name' , 'SMA' );
    set( hFig , 'position' , [0,6.75,6,4.75] );         
    
hFig = figure(2);     
    hold all
    plot( tdays , coe(2,:) , 'LineWidth' , 1.2 );
    hXLabel = xlabel( [ '$\textrm{Time [days]}', '$' ] , ...
                        'Interpreter', 'LaTex' );    
    hYLabel = ylabel( [ '$e\,', '$' ] , ...
                        'Interpreter', 'LaTex' , 'Rotation', 90 );           
    set( gca , 'FontName' , 'Helvetica' );
    set( hXLabel , 'FontSize' , 12 )
    set( hYLabel , 'FontSize' , 14 )
    set( hFig , 'units' , 'inches' , ... 
         'NumberTitle' , 'off' , 'Name' , 'ECC' );
    set( hFig , 'position' , [6.0,6.75,6,4.75] );         

hFig = figure(3);     
    hold all
    plot( tdays , coe(4,:) , 'LineWidth' , 1.2 );    
    hXLabel = xlabel( [ '$\textrm{Time [days]}', '$' ] , ...
                        'Interpreter', 'LaTex' );    
    hYLabel = ylabel( [ '$\Omega\, \textrm{ [deg]}', '$' ] , ...
                        'Interpreter', 'LaTex' , 'Rotation', 90 );
    set( gca , 'FontName' , 'Helvetica' );
    set( hXLabel , 'FontSize' , 12 )
    set( hYLabel , 'FontSize' , 14 )
    set( hFig , 'units' , 'inches' , ... 
         'NumberTitle' , 'off' , 'Name' , 'RAAN' );
    set( hFig , 'position' , [0,0.75,6,4.75] ); 

hFig = figure(4);     
    hold all
    plot( tdays , coe(5,:) , 'LineWidth' , 1.2 );    
    hXLabel = xlabel( [ '$\textrm{Time [days]}', '$' ] , ...
                        'Interpreter', 'LaTex' );    
    hYLabel = ylabel( [ '$\omega\, \textrm{ [deg]}', '$' ] , ...
                        'Interpreter', 'LaTex' , 'Rotation', 90 );
    set( gca , 'FontName' , 'Helvetica' );
    set( hXLabel , 'FontSize' , 12 )
    set( hYLabel , 'FontSize' , 14 )
    set( hFig , 'units' , 'inches' , ... 
         'NumberTitle' , 'off' , 'Name' , 'ARP' );
    set( hFig , 'position' , [6.0,0.75,6,4.75] ); 

toc
function [ value, isterminal, direction ] = reentryEvents( ~, y, const )

value = norm( y(1:3) ) - const.Re;      % Detect reentry = 0; that is, when r is less than the Earth's radius
isterminal = 1; 
direction = 0; 

end
function [condition, isterminal, direction] = EventFcnVR(t, x)

global R_earth 

h = sqrt(x(1)^2+x(2)^2+x(3)^2)-R_earth;                                     % Height [m]

% VR
condition(1) = h-100;                                                       % Event at 100 m  height
isterminal(1) = 1;                                                          % Ode stops
direction(1) = 0;



end
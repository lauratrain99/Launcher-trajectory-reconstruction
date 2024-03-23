function [condition, isterminal, direction] = EventFcnCP(t, x)

global R_earth w_earth r0Vect AoA
h = sqrt(x(1)^2 + x(2)^2 + x(3)^2) - R_earth;                               % Height [m]
r = x(1:3);
v = x(4:6);
v_r = v-(cross(w_earth,r));
kick = 0.2*pi/180;

rHat = r0Vect/norm(r0Vect);
eastVect = (cross([0 0 1], rHat))';                                         % East Vector
eastHat = eastVect/norm(eastVect);                                          % East Versor
uVect = cos(kick)*rHat'+sin(kick)*eastHat;                                  % Thrust Vector
uHat = uVect/norm(uVect);                                                   % Thrust Versor
flight_path = acos(dot(v_r,eastHat)/(norm(v_r)*norm(eastHat)));             % Flight path angle
pitch = acos(dot(uHat,eastHat)/(norm(uHat)*norm(eastHat)));                 % Pitch angle

AoA = flight_path-pitch;
condition(1) = flight_path-pitch;                                           % Event at reached desired kick angle
isterminal(1) = 1;                                                          % Ode stops
direction(1) = 0;

if condition(1) == h
    condition(1) = h;                                                       % Event at crash
    isterminal(1) = 1;                                                      % Ode stops
    direction(1) = 0;
    error('Launchers hit the ground!')
end

end
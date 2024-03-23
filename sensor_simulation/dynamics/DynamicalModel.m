function [dx] = DynamicalModel(t, x, phase_code)

global mu g0 Isp R_earth w_earth mdot P_e A_e A_a M0 cD0 m0 te_VR rHat ...
    pitch_deg uVect r0Vect

C = Isp*g0;
C_1 = C(1);

r = x(1:3);                                                                 % Position Vector
v = x(4:6);                                                                 % Velocity Vector

h = norm(r)-R_earth;                                                        % Altitude [m]
[rho, P_a, ~, a] = expEarthAtm(h);                                          % Density, pressure, temperature and speed of sound
v_r = v-cross(w_earth,r);                                                   % Relative velocity [m/s]
Th = C_1*mdot+(P_e-P_a)*A_e;                                                % Thrust [N]

m = m0-mdot*t;                                                              % Mass [kg]
M = norm(v_r)/a;                                                            % Mach number
cD = interp1(M0, cD0, M, 'linear', 'extrap');                               % Drag coefficient


% Dynamical eqs

rdot = v;

if phase_code == 1
    vdot = -mu/norm(r)^3*r+Th/m*r/norm(r)- ...
           0.5*rho*A_a/m*cD*v_r*norm(v_r);
    
elseif phase_code == 2
    w_pitch = 1*pi/180;                                                     % Pitch angular velocity [rad/s]
    pitch=w_pitch*(t-te_VR);                                                % Pitch angle [rad]
    
    eastVect = (cross([0 0 1], rHat))';                                     % East Vector
    uVect = cos(w_pitch*t)*r+sin(w_pitch*t)*eastVect;                       % Thrust Vector
    
    vdot = -mu/norm(r)^3*r+Th/m*uVect/norm(uVect)- ...
           0.5*rho*A_a/m*cD*v_r*norm(v_r);
    
    pitch_deg=pitch*180/pi;                                                 % Pitch angle [deg]
    
elseif phase_code == 3    
    kick = 0.2*pi/180;                                                      % Kick angle
    
    eastVect = (cross([0 0 1], r0Vect))';                                   % East Vector
    uVect=cos(kick)*r0Vect'+sin(kick)*eastVect;                             % Thrust Vector
    
    vdot = -mu/norm(r)^3*r+Th/m*uVect/norm(uVect)- ...
        0.5*rho*A_a/m*cD*v_r*norm(v_r);
    
elseif phase_code == 4    
    
    vdot = -mu/norm(r)^3*r+Th/m*v_r/norm(v_r)- ...
        0.5*rho*(A_a/m)*cD*v_r*norm(v_r);

end

dx = [rdot; vdot; -mdot];

end
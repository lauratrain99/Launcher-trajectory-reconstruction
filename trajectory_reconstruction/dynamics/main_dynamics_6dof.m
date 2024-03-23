%% Dynamic simulation

format longG
close all
clear 
clc

%% PARAMETERS

global mu g0 R_earth w_earth h_orbit lat m_pl N Isp eps P_e A_a A_e ...
    M0 cD0 m0 m_str a mdot mf1 mf2 mi2 te_VR te_CPR rHat pitch_deg uVect ...
    r0Vect AoA 

% Earth Constants
mu= 398602*10^9;                                                            % Gravitational parameter [m^3/s^2]
g0 = 9.80665;                                                               % Gravity constant [m/s^2]
w_earth = [0; 0; 2*pi/86164];                                               % Earth angular velocity [rad/s]
R_earth = 6378137;                                                          % Earth radius [m]

% Mission Requirements Data
h_orbit = 500000;                                                           % Orbit altitude [m]
lat = 5.2*pi/180;                                                           % Kourou latitude [rad]
a = h_orbit+R_earth;                                                        % Semi-major axis [m]
a_max1 = 7*g0;                                                              % First stage maximum acceleration [m/s^2]
a_max2 = 4*g0;                                                              % Second stage maximum acceleration [m/s^2]
m_pl = 300;                                                                 % Payload mass [kg]

% Launcher Data
N = 2;                                                                      % Number of stages 
eps = [0.11 0.14];                                                          % Structural Coefficient 
Isp = [310 330];                                                            % Vacuum Specific Impulse [s]

A_e = 0.3;                                                                  % Exhaust nozzle area [m^2] 
P_e = 40000;                                                                % Nozzle exit pressure [Pa]
A_a = 1;                                                                    % Reference aerodynamic area [m^2]
P_dyn_max = 45000;                                                          % Maximum dynamic pressure
M0 = [0.2 0.5 0.8 1.2 1.5 1.75 2 2.25 2.5 2.75 3 3.5 4 4.5 5 5.5 6 6.5];    % Mach number               
cD0 = [0.27 0.26 0.25 0.5 0.46 0.44 0.41 0.39 0.37 0.35 0.33...
      0.3 0.28 0.26 0.24 0.23 0.22 0.21];                                   % Drag coefficient
acc_VR = [0,0,0];
acc_CPR = [0,0,0];
acc_CP = [0,0,0];
acc_GT = [0,0,0];

t_VRacc = 0;
t_CPRacc = 0;
t_CPacc = 0;
t_GTacc = 0;

%% STAGING 

%The best distribution among the stages is furnished by @Staging that
%exploits the stage optimization by mean of the Lagrange multiplier method

[m0, m_subR, m_stage, m_str, m_prop, DV_req, DV_stage, lambda, ...
    MR] = Staging(N, Isp, eps, m_pl, h_orbit);                              % Staging

mf1 = m0-m_prop(1);                                                         % Mass at first stage burn out [kg]
mi2 = m0-m_stage(1);                                                        % Mass at second stage ignition [kg]
mf2 = mi2-m_prop(2);                                                        % Mass at second stage burn out [kg]

Th1 = mf1*a_max1;                                                           % First stage adapted thrust [N]
Th2 = mf2*a_max2;                                                           % Second stage thrust [N]

C = Isp*g0;                                                                 % Effective exhaust velocity [m/s]
mdot = Th1/C(1);                                                            % First stage Mass flow rate [kg/s]
t_bo1 = m_prop(1)/mdot;                                                     % First stage burning time [s]


%% INITIAL STATE VECTOR

r0Vect = [(R_earth)*cos(lat) 0 (R_earth)*sin(lat)];                         % Initial position vector [m]
rHat = r0Vect/norm(r0Vect);                                                 % Initial position versor 
eastVect = cross([0 0 1], rHat);                                            % East vector 
eastHat = eastVect/norm(eastVect);                                          % East versor 

v0Vect = cross(w_earth, r0Vect);                                            % Initial velocity vector [m/s]
vHat = v0Vect/norm(v0Vect);                                                 % Initial velocity versor [m/s]


%% VERTICAL RISING (VR)

phase_code = 1;                                                             % Phase identifier

tspan_VR = [0 100];                                                         % VR integration time [s]
x0_VR = [r0Vect v0Vect m0];                                                 % VR initial state vector

options_VR = odeset('Events', @(t,x) EventFcnVR(t, x), ...
    'AbsTol', 1E-8, 'RelTol', 1E-8);                                        % VR integration options
[t_VR, x_VR, te_VR, xe_VR, ie_VR] = ode45(@(t,x) DynamicalModel(t, x, phase_code),...
    tspan_VR, x0_VR,options_VR);                                            % VR state integration

H_VR = sqrt(xe_VR(1)^2+xe_VR(2)^2+xe_VR(3)^2)-R_earth;                      % VR altitude check [m]

for i = 1:length(x_VR)
    h_VR(i) = norm(x_VR(i,1:3))-R_earth;                                    % VR altitude [m]
end

r_VR = x_VR(:,1:3);                                                         % VR position vectors [m]
v_VR = x_VR(:,4:6);                                                         % VR velocity vectors [m/s]
m_VR = x_VR(:,7);                                                           % VR mass [kg]

rfVect_VR = xe_VR(1:3);                                                     % VR final position vector [m]
vfVect_VR = xe_VR(4:6);                                                     % VR final velocity vector [m/s]
mf_VR = xe_VR(7);                                                           % VR final mass [kg]

for i = 1:length(t_VR)
    pos_VR(i) = norm(r_VR(i,:));                                            % VR position [m]
    vel_VR(i) = norm(v_VR(i,:));                                            % VR velocity [m/s]
    vrelVector_VR(i,:) = v_VR(i,:)-cross(w_earth, r_VR(i,:));               % VR relative velocity vector [m/s]
    velrel_VR(i) = norm(vrelVector_VR(i,:));                                % VR relative velocity [m/s]
end
h_VR=h_VR';
pos_VR  =pos_VR';                                        
vel_VR =vel_VR';                                     
vrelVector_VR =vrelVector_VR';       
velrel_VR=velrel_VR';

% figure(1)
% subplot(1,3,1)
% xlabel('{t}  [s]')
% ylabel('{h}  [km]')
% title('Altitude')
% plot(t_VR, h_VR*1E-3,'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% subplot(1,3,2)
% plot(t_VR, vel_VR*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V}  [km/s]')
% title('Velocity')
% subplot(1,3,3)
% plot(t_VR, velrel_VR*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V_{rel}}  [km/s]')
% title('Relative Velocity')
% sgtitle('Vertical Rising')


%% CONSTANT PITCH RATE(CPR)

phase_code = 2;                                                             % Phase identifier

w_pitch = 1*pi/180;
kick = 0.2*pi/180;                                                          % Kick angle (desired pitch angle)
deltat_CPR = kick/w_pitch;

tspan_CPR = [te_VR te_VR+deltat_CPR];                                       % CPR integration time [s]
x0_CPR = [rfVect_VR vfVect_VR mf_VR];                                       % CPR initial state vector

options_CPR = odeset('AbsTol', 1E-8, 'RelTol', 1E-8);                       % CPR integration options
[t_CPR, x_CPR] = ode45(@(t,x) DynamicalModel(t, x, phase_code),...
    tspan_CPR, x0_CPR,options_CPR);                                         % CPR state integration

xe_CPR = x_CPR(end,:);                                                      % Final state vector
te_CPR = t_CPR(end);                                                        % Final time [s]

H_CPR = sqrt(xe_CPR(1)^2+xe_CPR(2)^2+xe_CPR(3)^2)-R_earth;                  % CPR altitude check [m]
aus=size(x_CPR);                                                            % Auxiliary variable
for i = 1:aus(1)
    h_CPR(i) = norm(x_CPR(i,1:3))-R_earth;                                  % CPR altitude [m]
end

r_CPR = x_CPR(:,1:3);                                                       % CPR position vectors [m]
v_CPR = x_CPR(:,4:6);                                                       % CPR velocity vectors [m/s]
m_CPR = x_CPR(:,7);                                                         % CPR mass [kg]

rfVect_CPR = xe_CPR(1:3);                                                   % CPR final position vector [m]
vfVect_CPR = xe_CPR(4:6);                                                   % CPR final velocity vector [m/s]
mf_CPR = xe_CPR(7);                                                         % CPR final mass [kg]

for i = 1:length(t_CPR)
    pos_CPR(i) = norm(r_CPR(i,:));                                          % CPR position [m]
    vel_CPR(i) = norm(v_CPR(i,:));                                          % CPR velocity [m/s]
    vrelVector_CPR(i,:) = v_CPR(i,:)-cross(w_earth, r_CPR(i,:));            % CPR relative velocity vector [m/s]
    velrel_CPR(i) = norm(vrelVector_CPR(i,:));                              % CPR relative velocity [m/s]
end

% figure(2)
% subplot(1,3,1)
% plot(t_CPR, h_CPR*1E-3,'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{h}  [km]')
% title('Altitude')
% subplot(1,3,2)
% plot(t_CPR, vel_CPR*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V}  [km/s]')
% title('Velocity')
% subplot(1,3,3)
% plot(t_CPR, velrel_CPR*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V_{rel}}  [km/s]')
% title('Relative Velocity')
% sgtitle('Constant Pitch Rate')


%% COSTANT PITCH (CP)

phase_code = 3;                                                             % Phase identifier

tspan_CP = [te_CPR 100];                                                    % CP integration time [s]
x0_CP = [rfVect_CPR vfVect_CPR mf_CPR];                                     % CP initial state vector

options_CP = odeset('Events', @(t,x) EventFcnCP(t, x), ...
    'AbsTol', 1E-8, 'RelTol', 1E-8);                                        % CR integration options

[t_CP, x_CP, te_CP, xe_CP, ie_CP] = ode45(@(t,x) DynamicalModel(t, x, phase_code),...
    tspan_CP, x0_CP,options_CP);                                            % CP state integration


H_CP = sqrt(xe_CP(1)^2+xe_CP(2)^2+xe_CP(3)^2)-R_earth;                      % CP altitude check [m]
aus=size(x_CP);
for i = 1:aus(1)
    h_CP(i) = norm(x_CP(i,1:3))-R_earth;                                    % CP altitude [m]
end

r_CP = x_CP(:,1:3);                                                         % CP position vectors [m]
v_CP = x_CP(:,4:6);                                                         % CP velocity vectors [m/s]
m_CP = x_CP(:,7);                                                           % CP mass [kg]

rfVect_CP = xe_CP(1:3);                                                     % CP final position vector [m]
vfVect_CP = xe_CP(4:6);                                                     % CP final velocity vector [m/s]
mf_CP = xe_CP(7);                                                           % CP final mass [kg]

for i = 1:length(t_CP)
    pos_CP(i) = norm(r_CP(i,:));                                            % CP position [m]
    vel_CP(i) = norm(v_CP(i,:));                                            % CP velocity [m/s]
    vrelVector_CP(i,:) = v_CP(i,:)-cross(w_earth, r_CP(i,:));               % CP relative velocity vector [m/s]
    velrel_CP(i) = norm(vrelVector_CP(i,:));                                % CP relative velocity [m/s]
end

% figure(3)
% subplot(1,3,1)
% plot(t_CP, h_CP*1E-3,'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{h}  [km]')
% title('Altitude')
% subplot(1,3,2)
% plot(t_CP, vel_CP*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V}  [km/s]')
% title('Velocity')
% subplot(1,3,3)
% plot(t_CP, velrel_CP*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V_{rel}}  [km/s]')
% title('Relative Velocity')
% sgtitle('Constant Pitch')


%% GRAVITY TURN (GT)

phase_code = 4;                                                             % Phase identifier

tspan_GT = [te_CP t_bo1];                                                   % GT integration time [s]
x0_GT = [rfVect_CP vfVect_CP mf_CP];                                        % GT initial state

options_GT = odeset('AbsTol', 1E-8, 'RelTol', 1E-8); 
[t_GT, x_GT] = ode45(@(t,x) DynamicalModel(t, x, ...
    phase_code), tspan_GT, x0_GT, options_GT);                              % GT state integration

xe_GT = x_GT(end,:);
te_GT = t_GT(end);

H_GT = norm(xe_GT(1:3))-R_earth;                                            % GT altitude check [m]
for i = 1:length(x_GT)
    h_GT(i) = norm(x_GT(i,1:3))-R_earth;                                    % GT altitude [m]
end

rf_GT = xe_GT(1:3);                                                         % GT final position [m]
vf_GT = xe_GT(4:6);                                                         % GT final velocity [m/s]
mf_GT = xe_GT(7);                                                           % GT final mass [kg]

r_GT = x_GT(:,1:3);                                                         % GT position vector [m]
v_GT = x_GT(:,4:6);                                                         % GT velocity vector [m/s]
m_GT = x_GT(:,7);                                                           % GT mass [kg]

for i = 1:length(t_GT)
    pos_GT(i) = norm(r_GT(i,:));                                            % GT position [m]
    vel_GT(i) = norm(v_GT(i,:));                                            % GT velocity [m/s]
    vrelVector_GT(i,:) = v_GT(i,:)-cross(w_earth, r_GT(i,:));               % GT relative velocity vector [m/s]
    velrel_GT(i) = norm(vrelVector_GT(i,:));                                % GT relative velocity [m/s]
end

% figure(4)
% subplot(1,3,1)
% plot(t_GT, h_GT*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{h}  [km]')
% title('Altitude')
% subplot(1,3,2)
% plot(t_GT, vel_GT*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V}  [km/s]')
% title('Velocity')
% subplot(1,3,3)
% plot(t_GT, velrel_GT*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{V_{rel}}  [km/s]')
% title('Relative Velocity')
% sgtitle('Gravity Turn')


%% STATE EVOLUTION

t = [t_VR', t_CPR', t_CP', t_GT'];                                          % Time vector [s]
r = [pos_VR', pos_CPR, pos_CP, pos_GT];                                     % Position vector [m]
v = [vel_VR', vel_CPR, vel_CP, vel_GT];                                     % Velocity vector [m/s]
m = [m_VR', m_CPR', m_CP', m_GT'];                                          % Mass vector [kg]
h = [h_VR', h_CPR, h_CP, h_GT]';                                            % Altitude vector [m]
vrel = [velrel_VR', velrel_CPR, velrel_CP, velrel_GT];                      % Relative velocity vector [m/s]

% figure(5)
% 
% subplot(1,5,1)
% plot(t, h*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{h}  [km]')
% title('Altitude')
% subplot(1,5,2)
% plot(t, r*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{r}  [km]')
% title('Range')
% subplot(1,5,3)
% plot(t, v*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{v}  [km/s]')
% title('Velocity')
% subplot(1,5,4)
% plot(t, m, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{m}  [kg]')
% title('Mass')
% subplot(1,5,5)
% plot(t, vrel*1E-3, 'Color', '[0.6350 0.0780 0.1840]')
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{v_{r}}  [km/s]')
% title('Relative Velocity')


%% DYNAMICAL PRESSURE CHECK

[rho_VR, Pa_VR, ~, a_VR] = expEarthAtm(h_VR);                               % VR atmospheric properties
P_dyn_VR_Vector = 0.5*rho_VR'.*(vrelVector_VR.^2);                          % VR dynamic pressure Vector [Pa]
for i = 1:length(t_VR)
   P_dyn_VR(i) = norm(P_dyn_VR_Vector(:,i));                                % VR dynamic pressure [Pa]
end

[rho_CPR, Pa_CPR, ~, a_CPR] = expEarthAtm(h_CPR);                           % CPR atmospheric properties
P_dyn_CPR_Vector = 0.5*rho_CPR'.*(vrelVector_CPR.^2);                       % CPR dynamic pressure Vector [Pa]
for i = 1:length(t_CPR)
   P_dyn_CPR(i) = norm(P_dyn_CPR_Vector(i,:));                              % CPR dynamic pressure [Pa]
end

[rho_CP, Pa_CP, ~, a_CP] = expEarthAtm(h_CP);                               % CP atmospheric properties
P_dyn_CP_Vector = 0.5*rho_CP'.*(vrelVector_CP.^2);                          % CP dynamic pressure Vector [Pa]
for i = 1:length(t_CP)
   P_dyn_CP(i) = norm(P_dyn_CP_Vector(i,:));                                % CP dynamic pressure [Pa]
end

[rho_GT, Pa_GT, ~, a_GT] = expEarthAtm(h_GT);                               % GT atmospheric properties
P_dyn_GT_Vector = 0.5*rho_GT'.*(vrelVector_GT.^2);                          % GT dynamic pressure Vector [Pa]
for i = 1:length(t_GT)
   P_dyn_GT(i) = norm(P_dyn_GT_Vector(i,:));                                % GT dynamic pressure [Pa]
end

P_dyn = [P_dyn_VR, P_dyn_CPR, P_dyn_CP, P_dyn_GT];                          % Dynamic pressure Vector [Pa]  

M_VR = vrelVector_VR./a_VR';                                                % VR Mach number [ ]
M_CPR = vrelVector_CPR./a_CPR';                                             % CPR Mach number [ ]
M_CP = vrelVector_CP./a_CP';                                                % CP Mach number [ ]
M_GT = vrelVector_GT./a_GT';                                                % GT Mach number [ ]

% figure(6)
% plot(t_VR, P_dyn_VR*10^-3, 'Color', '[0.4940 0.1840 0.5560]')
% hold on
% plot(t_CPR, P_dyn_CPR*10^-3, 'Color', '[0.6350 0.0780 0.1840]')
% hold on
% plot(t_CP, P_dyn_CP*10^-3, 'Color', '[0.8500 0.3250 0.0980]')
% hold on
% plot(t_GT, P_dyn_GT*10^-3, 'Color', '[0.9290 0.6940 0.1250]')
% hold on
% yline(P_dyn_max*10^-3, 'Color', '[0.4660 0.6740 0.1880]','LineWidth', 2);
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{P_{dyn}}  [kPa]')
% ylim([0 45])
% legend('Dynamic Pressure VR', 'Dynamic Pressure CPR', 'Dynamic Pressure CP',...
%        'Dynamic Pressure GT', 'Maximum Dynamic Pressure')
% title({'Dynamic Pressure'})


%% GRAVITY LOSSES

vrelVector = [vrelVector_VR, vrelVector_CPR', vrelVector_CP', vrelVector_GT'];

for i = 2:length(t)
    flight_path(i) = acos(dot(vrelVector(1:3,i),eastHat')/...               %Flight path angle Vector
                    (norm(vrelVector(1:3,i))*norm(eastHat')));
end

rVect = [r_VR', r_CPR', r_CP', r_GT'];
for i = 2:length(t)
    acc_G(i) = mu/(norm(rVect(1:3,i))^2)*sin(flight_path(i));               % Gravity acceleration [m/s^2]
end

DVLoss_G = trapz(t, acc_G);                                                 % Gravity losses [m/s]


%% AERODYNAMIC LOSSES

for i = 1:length(M_VR)
    cD_VR(i) = interp1(M0, cD0, M_VR(i), 'linear', 'extrap');               % VR drag coefficient [ ]
    acc_D_VR(i) = 0.5*A_a*rho_VR(i)*(velrel_VR(i)^2)*cD_VR(i)/m_VR(i);      % VR drag acceleration [m/s^2]
end

for i = 1:length(M_CPR)
    cD_CPR(i) = interp1(M0, cD0, M_CPR(i), 'linear', 'extrap');             % CPR drag coefficient [ ]
    acc_D_CPR(i) = 0.5*A_a*rho_CPR(i)*(velrel_CPR(i)^2)*cD_CPR(i)/m_CPR(i); % CPR drag acceleration [m/s^2]
end

for i = 1:length(M_CP)
    cD_CP(i) = interp1(M0, cD0, M_CP(i), 'linear', 'extrap');               % CP drag coefficient [ ]
    acc_D_CP(i) = 0.5*A_a*rho_CP(i)*(velrel_CP(i)^2)*cD_CP(i)/m_CP(i);      % CP drag acceleration [m/s^2]
end

for i = 1:(length(M_GT)-20)                                                 % Last 20 values of Mach vector are inf (expEarthAtm error)
    cD_GT(i) = interp1(M0, cD0, M_GT(i), 'linear', 'extrap');               % GT drag coefficient [ ]
    acc_D_GT(i) = 0.5*A_a*rho_GT(i)*(velrel_GT(i)^2)*cD_GT(i)/m_GT(i);      % GT drag acceleration [m/s^2]
end


DVLoss_D = trapz(acc_D_VR)+trapz(acc_D_CPR)+trapz(acc_D_CP)+...             % Drag losses [m/s]
           trapz(acc_D_GT);


%% PRESSURE LOSSES

for i = 1:length(Pa_VR)
    acc_P_VR(i) = Pa_VR(i)*A_e/m_VR(i);                                     % VR pressure acceleration [m/s^2]
end

for i = 1:length(Pa_CPR)
    acc_P_CPR(i) = Pa_CPR(i)*A_e/m_CPR(i);                                  % CPR pressure acceleration [m/s^2]
end

for i = 1:length(Pa_CP)
    acc_P_CP(i) = Pa_CP(i)*A_e/m_CP(i);                                     % CP pressure acceleration [m/s^2]
end

for i = 1:length(Pa_GT)
    acc_P_GT(i) = Pa_GT(i)*A_e/m_GT(i);                                     % GT pressure acceleration [m/s^2]
end

DVLoss_P = trapz(acc_P_VR)+trapz(acc_P_CPR)+trapz(acc_P_CP)+ ...            % Pressure losses [m/s]
           trapz(acc_P_GT);
DVLoss = DVLoss_D + DVLoss_G + DVLoss_P;                                    % Total losses [m/s]


% Results Display
% format short
% fprintf('Computed Losses during the trajectory: \n \n')
% fprintf('Gravity Losses: %5.4e [km/s] \n' , DVLoss_G*10^-3)
% fprintf('Aerodynamic Losses: %5.4e [km/s] \n' , DVLoss_D*10^-3)
% fprintf('Pressure Losses: %5.4e [km/s] \n' , DVLoss_P*10^-3)
% fprintf('DeltaV Total Losses: %5.4f [km/s] \n' , DVLoss*10^-3)




%% QUATERNIONS & ANGULAR VELOCITIES

rVector = [r_VR; r_CPR; r_CP; r_GT];                                        % Phases position Vector
uVect_VR = r_VR;                                                            % Thrust Vector VR
uVect_CPR = cos(w_pitch*t_CPR).*r_CPR+sin(w_pitch*t_CPR).*eastVect;         % Thrust Vector CPR
uVect_CP = (cos(kick)*r0Vect+sin(kick)*eastVect).*t_CP;                     % Thrust Vector CP
uVect_GT = vrelVector_GT;                                                   % Thrust Vector GT

uVector = [uVect_VR; uVect_CPR; uVect_CP; uVect_GT];                        % Phases thrust Vector 

t = [t_VR', t_CPR(2:end)', t_CP(2:end)', t_GT(2:end)'];

for i=1:length(t)
    eastVector = cross([0 0 1], rVector(i,:));                              % East Vector
    C_Matrix = CosineDirectorMatrix(uVector(i,:),eastVector);               % Cosine Director matrix
    q0(i) = (1/2)*sqrt(1+C_Matrix(1,1)+C_Matrix(2,2)+C_Matrix(3,3));
    q123Vect(i,:) = (1/(4*q0(i)))*[C_Matrix(2,3)-C_Matrix(3,2) C_Matrix(3,1)-C_Matrix(1,3) C_Matrix(1,2)-C_Matrix(2,1)];

end

q0 = q0';
qVect = [q0 q123Vect];                                                      % Quaternions Vector

% figure(1)
% plot(t,q123Vect)
% hold on
% plot(t,q0)
% grid on
% grid minor
% xlabel('{t}  [s]')
% ylabel('{Q}')
% legend('q1', 'q2', 'q3','q0')
% title({'Quaternions'})


% q0_dot = diff(q0')./diff(t);
for i=1:length(t)-1
      q0_dot_(i) = (q0(i+1)-q0(i))./(t(i+1)-t(i));
      q1_dot_(i) = (q123Vect(i+1,1)-q123Vect(i,1))./(t(i+1)-t(i));
      q2_dot_(i) = (q123Vect(i+1,2)-q123Vect(i,2))./(t(i+1)-t(i));
      q3_dot_(i) = (q123Vect(i+1,3)-q123Vect(i,3))./(t(i+1)-t(i));
end

% delete spurious terms
q0_dot_(95) = q0_dot_(96);
q1_dot_(95) = q1_dot_(96);
q2_dot_(95) = q2_dot_(96);
q3_dot_(95) = q3_dot_(96);
q3_dot_(82) = q3_dot_(83);

q0_dot = movmean(q0_dot_,50);
q1_dot = movmean(q1_dot_,50);
q2_dot = movmean(q2_dot_,50);
q3_dot = movmean(q3_dot_,50);

q_dot = [q0_dot; q1_dot; q2_dot; q3_dot];     
q_dot_ = [q0_dot_; q1_dot_; q2_dot_; q3_dot_];   % Quaternion Derivatives Vector
% q_dot = q0_dot;


figure(1)
plot(t(2:end), q0_dot, 'k', t(2:end), q1_dot, 'r', t(2:end), q2_dot, 'g', t(2:end), q3_dot, 'b','LineWidth',2.5)
hold on
plot(t(2:end), q0_dot_, 'k', t(2:end), q1_dot_, 'r', t(2:end), q2_dot_, 'g', t(2:end), q3_dot_, 'b','LineWidth',0.5)
grid minor
xlabel('t [s]')
ylabel('Quaternion derivatives')
legend('$q0_{der,filt}$', '$q1_{der,filt}$', '$q2_{der,filt}$', '$q3_{der,filt}$', ...
       '$q0_{der,non-filt}$', '$q1_{der,non-filt}$', '$q2_{der,non-filt}$', '$q3_{der,non-filt}$', ...
       'NumColumns',2,'FontSize',12)
title('Quaternion derivatives')

%%


w_att = zeros(1,4);                                                         % Initialize Angular Velocities
% q_conj = quatconj(qVect);
% w_att = 2*quatmultiply( q_dot', q_conj(2:end,:) );
for i=1:length(t)-1
    Q_matrix = [qVect(i,1) qVect(i,2) qVect(i,3) qVect(i,4);                % Quaternion products
                -qVect(i,2) qVect(i,1) qVect(i,4) -qVect(i,3);
                -qVect(i,3) -qVect(i,4) qVect(i,1) qVect(i,2);
                -qVect(i,4) qVect(i,3) -qVect(i,2) qVect(i,1)];
    w_att(i,1) = 0;
    w_att(i,:) = 2*Q_matrix*q_dot(:,i);                                     % Angular Velocities
    w_att(i,1) = 0;
end

w1 = movmean(w_att(:,2),50);
w2 = movmean(w_att(:,3),50);
w3 = movmean(w_att(:,4),50);


w = [w1, w2, w3];
figure(2)
plot(t(2:end), w1, 'b', t(2:end), w2, 'r', t(2:end), w3, 'g', 'LineWidth', 2)
hold on
plot(t(2:end), w_att(:,2), 'b', t(2:end), w_att(:,3), 'r', t(2:end), w_att(:,4), 'g', 'LineWidth', 0.5)
grid on
grid minor
xlabel('t [s]')
ylabel('$\omega [rad/s]$')
legend('$\omega_{x,filt}$', '$\omega_{y,filt}$', '$\omega_{z,filt}$', ...
       '$\omega_{x,non-filt}$', '$\omega_{x,non-filt}$', '$\omega_{x,non-filt}$', ...
       'NumColumns',2,'FontSize',14)
title('Angular Velocities')

%% NON_GRAVITATIONAL ACCELERATIONS

vrelVector = vrelVector';                                                   % Relative velocity Vector
vrelVector(41,:) = [];                                                      % Deleting repetitive components
vrelVector(83,:) = [];
vrelVector(96,:) = [];
uVector(41,:) = [];                                                         % Deleting repetitive components
uVector(83,:) = [];
uVector(96,:) = [];
m(:,41) = [];                                                               % Deleting repetitive components
m(:,83) = [];
m(:,96) = [];
[rho, P_a, ~, a] = expEarthAtm(h);                                          % Atmospheric Model
rho(41,:) = [];                                                             % Deleting repetitive components
rho(83,:) = [];
rho(96,:) = [];
rho(470:end) = rho(469);
rvec = [r_VR; r_CPR; r_CP; r_GT];
vvec = [v_VR; v_CPR; v_CP; v_GT];

for i=1:length(t)
    M = norm(vrelVector(i,:))/a;                                            % Mach number
    cD = interp1(M0, cD0, M, 'linear', 'extrap');                           % Drag coefficient
    cD(:,41) = [];                                                          % Deleting repetitive components
    cD(:,83) = [];
    cD(:,96) = [];
    if t(i) < t_VR(end)
        accVect(i,:) = Th1/m(i)*rvec(i,:)/norm(rvec(i,:))-...             % Non-gravitational Accelerations
        0.5*rho(i)*A_a/m(i)*cD(i)*vrelVector(i,:)*norm(vrelVector(i,:));
    elseif t(i) < t_CPR(end)
        accVect(i,:) = Th1/m(i)*uVector(i,:)/norm(uVector(i,:))-...             % Non-gravitational Accelerations
        0.5*rho(i)*A_a/m(i)*cD(i)*vrelVector(i,:)*norm(vrelVector(i,:));
    elseif t(i) < t_CP(end)
        accVect(i,:) = Th1/m(i)*uVector(i,:)/norm(uVector(i,:))-...             % Non-gravitational Accelerations
        0.5*rho(i)*A_a/m(i)*cD(i)*vrelVector(i,:)*norm(vrelVector(i,:));
    else
        accVect(i,:) = Th1/m(i)*vrelVector(i,:)/norm(vrelVector(i,:))-...             % Non-gravitational Accelerations
        0.5*rho(i)*A_a/m(i)*cD(i)*vrelVector(i,:)*norm(vrelVector(i,:));
    end
end


figure(3)
plot(t, accVect)
grid on
grid minor
xlabel('t [s]')
ylabel('$acc_{non_G} [m/s^2]$')
legend('$acc_x$', '$acc_y$', '$acc_z$')
title('Non-gravitational accelerations: inertial referece frame')

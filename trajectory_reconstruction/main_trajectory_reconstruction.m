% This script is the main file for trajectory reconstruction

clear;clc;

addpath test/
addpath navigation/
addpath dynamics/

set(groot, 'defaultTextInterpreter',            'latex');
set(groot, 'defaultAxesTickLabelInterpreter',   'latex'); 
set(groot, 'defaultLegendInterpreter',          'latex');
set(groot, 'defaultLegendLocation',             'northeast');

%% Run IMU parameters and dynamic simulation

run('dynamics/main_dynamics_6dof.m')
run('navigation/imu_params.m')

%% perpare variables to send to simulink

% Eliminate eps variable not to confuse simulink
clear eps;

% hardcode deletion NaN in trajectory
tf = t(488);
t_ = t(1:488);
clear t;
quat = [q0, q123Vect];

% initiate variables to corrupt with IMU
accVect_nongrav = zeros(488,3);
accVect_nongrav(:,1) = accVect(1:488,1);
accVect_nongrav(:,2) = accVect(1:488,2);
accVect_nongrav(:,3) = accVect(1:488,3);

wVect_ = zeros(488,3);
wVect_(:,1) = t_;
wVect_(:,2) = w(:,1);
wVect_(:,3) = w(:,2);
wVect_(:,4) = w(:,3);

% rotate non gravitational acceleration to body reference system
acc_Bnongrav = zeros(length(accVect_nongrav),3);
for i = 1:length(accVect_nongrav)
      acc_Bq = quatmultiply(quatmultiply(quatinv(quat(i,:)),[0,accVect_nongrav(i,:)]),quat(i,:));
      acc_Bnongrav(i,1) = t_(i);
      acc_Bnongrav(i,2:4) = acc_Bq(2:4);
end

%% run simulink model
simOut = sim('navigation/trajectory_IMU','SaveOutput','on');

% rearrange vectors for plotting
acc_meas = zeros(length(simOut.acc_meas.Data(1,1,:)),3);
omega_meas = simOut.omega_meas.Data;

for i = 1:length(acc_meas)
    acc_meas(i,1) = simOut.acc_meas.Data(1,1,i);
    acc_meas(i,2) = simOut.acc_meas.Data(1,2,i);
    acc_meas(i,3) = simOut.acc_meas.Data(1,3,i);
end


figure(4)
plot(simOut.acc_meas.Time, acc_meas(:,1),'r', ...
     simOut.acc_meas.Time, acc_meas(:,2),'b', ...
     simOut.acc_meas.Time, acc_meas(:,3),'g','LineWidth',1)
title("Non-grav acceleration measurements",'FontSize',14)
legend("$a_{meas,Bx}$","$a_{meas,By}$","$a_{meas,Bz}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Acceleration[$m/s^2$]",'FontSize',14)
grid minor


figure(5)
plot(simOut.omega_meas.Time, omega_meas(:,1),'r', ...
     simOut.omega_meas.Time, omega_meas(:,2),'b', ...
     simOut.omega_meas.Time, omega_meas(:,3),'g','LineWidth',1)

title("Angular velocity measurements",'FontSize',14)
legend("$\omega_{meas,Bx}$","$\omega_{meas,By}$","$\omega_{meas,Bz}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Angular velocity [$^\circ/s$]",'FontSize',14)
grid minor

figure(6)
plot(t_, acc_Bnongrav(:,2),'r', ...
     t_, acc_Bnongrav(:,3),'b', ...
     t_, acc_Bnongrav(:,4),'g','LineWidth',1)

title("Non-grav acceleration from dynamics",'FontSize',14)
legend("$a_{dyn,Bx}$","$a_{dyn,By}$","$a_{dyn,Bz}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Acceleration [$m/s^2$]",'FontSize',14)
grid minor

figure(7)
plot(t_, wVect_(:,2)*180/pi,'r', ...
     t_, wVect_(:,3)*180/pi,'b', ...
     t_, wVect_(:,4)*180/pi,'g','LineWidth',1)

title("Angular velocity from dynamics",'FontSize',14)
legend("$\omega_{dyn,Bx}$","$\omega_{dyn,By}$","$\omega_{dyn,Bz}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Angular velocity [$^\circ/s$]",'FontSize',14)
grid minor


%% angular rates results

omega_meas_filt = simOut.omega_meas_filt.Data*pi/180;
t_of = simOut.omega_meas_filt.Time;

figure(10)
plot(t_, w1*180/pi, 'b', t_, w2*180/pi, 'r', t_, w3*180/pi, 'g', 'LineWidth', 2)
hold on
plot(t_of, omega_meas_filt(:,1)*180/pi, 'b', t_of, omega_meas_filt(:,2)*180/pi, 'r', t_of, omega_meas_filt(:,3)*180/pi, 'g', 'LineWidth', 0.5)
hold off
grid minor
xlabel('t [s]')
ylabel('$\omega  [^\circ/s] $')
legend('$\omega_{x}$', '$\omega_{y}$', '$\omega_{z}$', ...
       '$\omega_{x,rec}$', '$\omega_{y,rec}$', '$\omega_{z,rec}$', ...
       'NumColumns',2,'FontSize',14)
title('Reconstructed angular velocity')

w1_interp = interp1(t_,w1,t_of);
w2_interp = interp1(t_,w2,t_of);
w3_interp = interp1(t_,w3,t_of);

figure(8)
plot(t_of,abs(w1_interp - omega_meas_filt(:,1))*180/pi,'b', ...
     t_of,abs(w2_interp - omega_meas_filt(:,2))*180/pi,'r', ...
     t_of,abs(w3_interp - omega_meas_filt(:,3))*180/pi,'g','LineWidth',0.5)
title("Reconstructed angular velocity error")
legend('$\Delta \omega_x$','$\Delta \omega_y$','$\Delta \omega_z$','FontSize',14)
xlabel("t [s]")
ylabel("Angular velocity error [$ ^\circ /s$]")
grid minor
hold off


%% attitude results

% numerical integration to get quaternions from angular velocity
qinit = quat(1,:);
tspan = [t_(1),t_(end)];
[t_qmeas,q_meas] = ode45(@(t,q) qderivative(t,q,t_of,omega_meas_filt),tspan,qinit);

% initialize euler from dynamics and from actual measurements
eul_th = zeros(length(quat),3);
eul_meas = zeros(length(q_meas),3);

for i = 1:length(quat)
    eul_th(i,:) = quat2eul(quat(i,:)) * 180/pi;
end

for i = 1:length(q_meas)
    eul_meas(i,:) = quat2eul(q_meas(i,:)) * 180/pi;
end

figure(9)
plot(t_,eul_th(1:1:end-1,1),'b',t_,eul_th(1:1:end-1,2),'r',t_,eul_th(1:1:end-1,3),'k','LineWidth',2)
title("Reconstructed attitude")
hold on
plot(t_qmeas,eul_meas(:,1),'b',t_qmeas,eul_meas(:,2),'r',t_qmeas,eul_meas(:,3),'k','LineWidth',0.5)
legend('$\psi$','$\theta$','$\phi$','$\psi_{rec}$','$\theta_{rec}$','$\phi_{rec}$','NumColumns',2,'FontSize',14)
xlabel("t [s]")
ylabel("Euler angles [$ ^\circ $]")
grid minor
hold off

eul_th_interp(:,1) = interp1(t_,eul_th(1:1:end-1,1),t_qmeas);
eul_th_interp(:,2) = interp1(t_,eul_th(1:1:end-1,2),t_qmeas);
eul_th_interp(:,3) = interp1(t_,eul_th(1:1:end-1,3),t_qmeas);

figure(10)
plot(t_qmeas,abs(eul_meas(:,1) - eul_th_interp(:,1)),'b', ...
     t_qmeas,abs(eul_meas(:,2) - eul_th_interp(:,2)),'r', ...
     t_qmeas,abs(eul_meas(:,3) - eul_th_interp(:,3)),'k','LineWidth',2)
title("Reconstructed attitude error")
legend('$\Delta \psi$','$\Delta \theta$','$\Delta \phi$','FontSize',14)
xlabel("t [s]")
ylabel("Euler angles error [$ ^\circ $]")
grid minor
hold off

%% position and velocity reconstruction 

% get initial r,v values for integration
vx_ = vvec(1,1); vy_ = vvec(1,2); vz_ = vvec(1,3);
rx_ = rvec(1,1); ry_ = rvec(1,2); rz_ = rvec(1,3);

% initialize vectors
acc_Nnongrav_rec = zeros(length(acc_meas),3);
acc_N_rec = zeros(length(acc_meas),3);
t_measacc = simOut.omega_meas.Time;
v_rec(1,:) = [vx_,vy_,vz_];
r_rec(1,:) = [rx_,ry_,rz_];

% adapt q_meas to same size as acc_meas
q_meas_interp(:,1) = interp1(t_qmeas,q_meas(:,1),t_measacc);
q_meas_interp(:,2) = interp1(t_qmeas,q_meas(:,2),t_measacc);
q_meas_interp(:,3) = interp1(t_qmeas,q_meas(:,3),t_measacc);
q_meas_interp(:,4) = interp1(t_qmeas,q_meas(:,4),t_measacc);

for i = 1:length(acc_meas)
      
      % rotate acc-nongravB to acc-nongravN
      acc_Nq = quatmultiply(quatmultiply(q_meas_interp(i,:),[0,acc_meas(i,:)]),quatinv(q_meas_interp(i,:)));
      acc_Nnongrav_rec(i,:) = acc_Nq(2:4);
      
      % compute gravity contribution in N frame
      grav_contrib_n = -mu*r_rec(i,:)'/norm(r_rec(i,:)')^3;
      
      % add gravity contribution to linear acceleration
      acc_N_rec(i,:) = acc_Nnongrav_rec(i,:) + grav_contrib_n';
      
      % numerical integration of velocity and position
      if i < (length(acc_meas) - 1)
          vx_ = vx_ + (acc_N_rec(i,1) * (t_measacc(i+1) - t_measacc(i)));
          vy_ = vy_ + (acc_N_rec(i,2) * (t_measacc(i+1) - t_measacc(i)));
          vz_ = vz_ + (acc_N_rec(i,3) * (t_measacc(i+1) - t_measacc(i)));
          v_rec(i+1,:) = [vx_,vy_,vz_];

          rx_ = rx_ + (vx_ * (t_measacc(i+1) - t_measacc(i)));
          ry_ = ry_ + (vy_ * (t_measacc(i+1) - t_measacc(i)));
          rz_ = rz_ + (vz_ * (t_measacc(i+1) - t_measacc(i)));
          r_rec(i+1,:) = [rx_,ry_,rz_];
      else
          v_rec(i+1,:) = v_rec(i,:);
          r_rec(i+1,:) = r_rec(i,:);
      end
end

t = [t_VR', t_CPR', t_CP', t_GT'];
figure(11)
plot(t, vvec(:,1),'r', ...
     t, vvec(:,2),'b', ...
     t, vvec(:,3),'g','LineWidth',2)
hold on
plot(t_measacc, v_rec(1:end-1,1),'r', ...
     t_measacc, v_rec(1:end-1,2),'b', ...
     t_measacc, v_rec(1:end-1,3),'g','LineWidth',0.5)
hold off
xlim([t_measacc(1), t_measacc(end)])
title("Reconstructed velocity",'FontSize',14)
legend("$v_{x,true}$","$v_{y,true}$","$v_{z,true}$", ...
    "$v_{x,rec}$","$v_{y,rec}$","$v_{z,rec}$",'NumColumns',2,'FontSize',14) 
xlabel("Time [s]")
ylabel("Velocity [$m/s$]",'FontSize',14)
grid minor

figure(12)
sgtitle("Reconstructed position",'FontSize',14)
subplot(3,1,1)
plot(t, rvec(:,1),'r', 'LineWidth',2)
hold on
plot(t_measacc, r_rec(1:end-1,1),'r','LineWidth',0.5)
hold off
legend("$r_{x,true}$","$r_{x,rec}$",'FontSize',14)
xlim([t_measacc(1), t_measacc(end)])
xlabel("Time [s]")
ylabel("Position [$m$]",'FontSize',14)
grid minor

subplot(3,1,2)
plot(t, rvec(:,2),'b', 'LineWidth',2)
hold on
plot(t_measacc, r_rec(1:end-1,2),'b','LineWidth',0.5)
hold off
legend("$r_{y,true}$","$r_{y,rec}$",'FontSize',14)
xlim([t_measacc(1), t_measacc(end)])
xlabel("Time [s]")
ylabel("Position [$m$]",'FontSize',14)
grid minor

subplot(3,1,3)
plot(t, rvec(:,3),'g', 'LineWidth',2)
hold on
plot(t_measacc, r_rec(1:end-1,3),'g','LineWidth',0.5)
hold off
legend("$r_{z,true}$","$r_{z,rec}$",'FontSize',14)
xlim([t_measacc(1), t_measacc(end)])
xlabel("Time [s]")
ylabel("Position [$m$]",'FontSize',14)
grid minor

% correct spurious time positions
t(42)= t(42) + 0.001;
t(83)= t(83) + 0.001;
t(96)= t(96) + 0.001;

vx_interp = interp1(t,vvec(:,1),t_measacc);
vy_interp = interp1(t,vvec(:,2),t_measacc);
vz_interp = interp1(t,vvec(:,3),t_measacc);

rx_interp = interp1(t,rvec(:,1),t_measacc);
ry_interp = interp1(t,rvec(:,2),t_measacc);
rz_interp = interp1(t,rvec(:,3),t_measacc);

figure(13)
plot(t_measacc, abs(vx_interp - v_rec(1:end-1,1)),'r', ...
     t_measacc, abs(vy_interp - v_rec(1:end-1,2)),'b', ...
     t_measacc, abs(vz_interp - v_rec(1:end-1,3)),'g','LineWidth',1)
title("Reconstructed velocity error",'FontSize',14)
legend('$\Delta v_x$','$\Delta v_y$','$\Delta v_z$','FontSize',14) 
xlabel("Time [s]")
ylabel("Velocity error [$m/s$]",'FontSize',14)
grid minor

figure(14)
plot(t_measacc, abs(rx_interp - r_rec(1:end-1,1)),'r', ...
     t_measacc, abs(ry_interp - r_rec(1:end-1,2)),'b', ...
     t_measacc, abs(rz_interp - r_rec(1:end-1,3)),'g','LineWidth',1)
title("Reconstructed position error",'FontSize',14)
legend('$\Delta r_x$','$\Delta r_y$','$\Delta r_z$','FontSize',12) 
xlabel("Time [s]")
ylabel("Position error [$m$]",'FontSize',14)
grid minor


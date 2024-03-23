% Simulate navigation data using STIM300 IMU model

clear;clc;

addpath test/
addpath navigation/
addpath dynamics/

%% Run IMU parameters and dynamic simulation

run('dynamics/main_dynamic_6dof.m')
run('navigation/imu_params.m')

% Eliminate eps variable not to confuse simulink
clear eps;

% hardcode deletion NaN in trajectory
tf = t(469);
t_ = t(1:469);

accVect_ = zeros(469,3);
accVect_(:,1) = t_;
accVect_(:,2) = accVect(1:469,1);
accVect_(:,3) = accVect(1:469,2);
accVect_(:,4) = accVect(1:469,3);

wVect_ = zeros(469,3);
wVect_(:,1) = t_;
wVect_(:,2) = w_att(1:469,2);
wVect_(:,3) = w_att(1:469,3);
wVect_(:,4) = w_att(1:469,4);

% run simulink model
simOut = sim('navigation/trajectory_IMU','SaveOutput','on');

% rearrange vectors for plotting
acc_meas = zeros(length(simOut.acc_meas.Data(1,1,:)),3);
omega_meas = simOut.omega_meas.Data;

for i = 1:length(acc_meas)
    acc_meas(i,1) = simOut.acc_meas.Data(1,1,i);
    acc_meas(i,2) = simOut.acc_meas.Data(1,2,i);
    acc_meas(i,3) = simOut.acc_meas.Data(1,3,i);
end

% plot results
set(groot, 'defaultTextInterpreter',            'latex');
set(groot, 'defaultAxesTickLabelInterpreter',   'latex'); 
set(groot, 'defaultLegendInterpreter',          'latex');
set(groot, 'defaultLegendLocation',             'northeast');

figure(1)
plot(simOut.acc_meas.Time, acc_meas(:,1),'r', ...
     simOut.acc_meas.Time, acc_meas(:,2),'b', ...
     simOut.acc_meas.Time, acc_meas(:,3),'g','LineWidth',1)

title("Non-grav acceleration measurements",'FontSize',14)
legend("$a_{meas,x}$","$a_{meas,y}$","$a_{meas,z}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Acceleration[$m/s^2$]",'FontSize',14)
grid minor


figure(2)
plot(simOut.omega_meas.Time, omega_meas(:,1),'r', ...
     simOut.omega_meas.Time, omega_meas(:,2),'b', ...
     simOut.omega_meas.Time, omega_meas(:,3),'g','LineWidth',1)

title("Angular velocity measurements",'FontSize',14)
legend("$\omega_{meas,x}$","$\omega_{meas,y}$","$\omega_{meas,z}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Angular velocity [$^\circ/s$]",'FontSize',14)
grid minor

figure(3)
plot(accVect_(:,1), accVect_(:,2),'r', ...
     accVect_(:,1), accVect_(:,3),'b', ...
     accVect_(:,1), accVect_(:,4),'g','LineWidth',1)

title("Non-grav acceleration from dynamics",'FontSize',14)
legend("$a_{dyn,x}$","$a_{dyn,y}$","$a_{dyn,z}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Acceleration [$m/s^2$]",'FontSize',14)
grid minor

figure(4)
plot(wVect_(:,1), wVect_(:,2)*180/pi,'r', ...
     wVect_(:,1), wVect_(:,3)*180/pi,'b', ...
     wVect_(:,1), wVect_(:,4)*180/pi,'g','LineWidth',1)

title("Angular velocity from dynamics",'FontSize',14)
legend("$\omega_{dyn,x}$","$\omega_{dyn,y}$","$\omega_{dyn,z}$",'FontSize',14) 
xlabel("Time [s]")
ylabel("Angular velocity [$^\circ/s$]",'FontSize',14)
grid minor



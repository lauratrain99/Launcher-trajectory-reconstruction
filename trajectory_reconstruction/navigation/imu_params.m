% This script contains the parameters to run the IMU.slx
% Autors: Laura Train, Sergi Escribano, Max Kronwitter

% Common
freq = 2000; %Hz
Ts = 1/freq;
ang0 = [1/1000, 1/1000, 1/1000]; %rad
R_pitch  = [cos(ang0(2))  0  sin(ang0(2)); 0 1 0; -sin(ang0(2))  0   cos(ang0(2))]; 
R_roll   = [1 0 0; 0 cos(ang0(1)) -sin(ang0(1)); 0 sin(ang0(1)) cos(ang0(1))];
R_yaw   = [cos(ang0(3)) -sin(ang0(3)) 0; sin(ang0(3)) cos(ang0(3)) 0; 0 0 1];
misalignment = R_yaw*R_pitch*R_roll;

% misalignment in reference system
mis_ref = [2/180 * pi, 2/180 * pi, 2/180 * pi]; %rad
R_pitch_ref  = [cos(mis_ref(2))  0  sin(mis_ref(2)); 0 1 0; -sin(mis_ref(2))  0   cos(mis_ref(2))]; 
R_roll_ref   = [1 0 0; 0 cos(mis_ref(1)) -sin(mis_ref(1)); 0 sin(mis_ref(1)) cos(mis_ref(1))];
R_yaw_ref   = [cos(mis_ref(3)) -sin(mis_ref(3)) 0; sin(mis_ref(3)) cos(mis_ref(3)) 0; 0 0 1];
misalignment_ref = R_yaw_ref*R_pitch_ref*R_roll_ref;

% misalignment in imu platform and body axes
mis_ax = [2/180 * pi, 2/180 * pi, 2/180 * pi]; %rad
R_pitch_ax  = [cos(mis_ax(2))  0  sin(mis_ax(2)); 0 1 0; -sin(mis_ax(2))  0   cos(mis_ax(2))]; 
R_roll_ax   = [1 0 0; 0 cos(mis_ax(1)) -sin(mis_ax(1)); 0 sin(mis_ax(1)) cos(mis_ax(1))];
R_yaw_ax   = [cos(mis_ax(3)) -sin(mis_ax(3)) 0; sin(mis_ax(3)) cos(mis_ax(3)) 0; 0 0 1];
misalignment_ax = R_yaw_ref*R_pitch_ax*R_roll_ax;

% Gyroscope
gyro.fullrange = 400; %deg/s
gyro.scalefactor = 500/1000000; %ppm
gyro.nonlinearity_200 = 15/1000000; %ppm
gyro.nonlinearity_400 = 20/1000000; %ppm
gyro.xspline = [-gyro.fullrange, -gyro.fullrange/2, 0, gyro.fullrange/2, gyro.fullrange];
gyro.yspline = [-gyro.fullrange - gyro.nonlinearity_400*gyro.fullrange, ...
                -gyro.fullrange/2 - gyro.nonlinearity_200*gyro.fullrange/2, ...
                 0, ...
                 gyro.fullrange/2 - gyro.nonlinearity_200*gyro.fullrange/2, ...
                 gyro.fullrange - gyro.nonlinearity_400*gyro.fullrange];
gyro.biasinst = [0.3, 0.3, 0.3]/3600; %deg/s
gyro.noisepower = ([0.15, 0.15, 0.15]/60).^2; %(deg/s)^2/Hz

gyro.biasrunrun_min = 0; %deg/s
gyro.biasrunrun_max = 4/3600; %deg/s
gyro.biasrunrun = gyro.biasrunrun_min + (gyro.biasrunrun_max-gyro.biasrunrun_min).*rand(3,1);%deg/s

gyro.biastrimoffset_max = 1; %deg/s
gyro.biastrimoffset_min = -1; %deg/s
rng(1)
gyro.biastrimoffset = gyro.biastrimoffset_min + (gyro.biastrimoffset_max-gyro.biastrimoffset_min).*rand(3,1);


gyro.quantization = 2*gyro.fullrange/2^24;

% Accelerometer
acc.scalefactor = 200/1000000; %ppm
acc.nonlinearity = 100/1000000; %ppm
acc.fullrange = 5; %g
acc.xspline = [-acc.fullrange,  0,  acc.fullrange];
acc.yspline = [-acc.fullrange - acc.nonlinearity*acc.fullrange, ...
                0, ...
                acc.fullrange - acc.nonlinearity*acc.fullrange];
acc.biasinst = [0.02, 0.02, 0.02]/1000; %g
acc.noisepower = ([0.03, 0.03, 0.03]/60/9.81).^2; %g^2/Hz


acc.biasswitchonoff_min = -0.38/1000; %g
acc.biasswitchonoff_max = 0.38/1000; %g
acc.biasswitchonoff = acc.biasswitchonoff_min + (acc.biasswitchonoff_max-acc.biasswitchonoff_min).*rand(3,1);%deg/s

acc.biastrimoffset_max = -50/1000; %g
acc.biastrimoffset_min = 50/1000; %g
rng(2)
acc.biastrimoffset = acc.biastrimoffset_min + (acc.biastrimoffset_max-acc.biastrimoffset_min).*rand(3,1);

acc.quantization = 2*acc.fullrange/2^24;


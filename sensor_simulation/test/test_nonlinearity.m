% This script tests the non-linearity implementation


clear;clc;

x = linspace(-4,4,100);

gyro_nonx = [-4,-2,0,2,4];
gyro_nony = [gyro_nonx(1) + gyro_nonx(1)*0.2, gyro_nonx(2) + gyro_nonx(2)*0.1,0,gyro_nonx(4) - gyro_nonx(4)*0.1,gyro_nonx(5) - gyro_nonx(5)*0.2];
yy = spline(gyro_nonx,gyro_nony,x);
figure(1)
plot(x,x,'r')
hold on
plot(x,yy,'b--')
hold on
plot(gyro_nonx,gyro_nony,'*r')
xlim([-4,4])
legend("Ideal","Non-linear","Interpolation points")
xlabel("Input")
ylabel("Output")
hold off
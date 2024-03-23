% This script tests the scale factor accuracy implementation


clear;clc;

x = linspace(-4,4,100);
y1 = x - x*0.1;
y2 = x + x*0.1;

figure(1)
plot(x,x,'r')
hold on
plot(x,y1,'b--')
hold on
plot(x,y2,'b--')
xlim([-4,4])
legend("Ideal","Scale factor")
xlabel("Input")
ylabel("Output")
hold off
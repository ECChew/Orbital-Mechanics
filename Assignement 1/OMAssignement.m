%Orbital Mechanics assignement
%Julien Huynh 10/19/2019 
clear all; clc;
[t,x1] = ode45(@f,[0:10:86400], [7078 0 0 0 0 7.5044]);
[t,x2] = ode45(@f,[0:10:86400], [7078 0 0 0 0 9.5044]);
[t,x3] = ode45(@f,[0:10:86400], [7078 0 0 0 0 11.5044]);
[t,x4] = ode45(@f,[0:10:86400], [7078 0 0 0 0 15.5044]);
%2D Plots
figure(1);
%First orbit
subplot(2,2,1);plot(x1(:,1),x1(:,5));
xlabel('x(km)');ylabel('z(km)');
axis([-1 1 -1 1]*1e4);
title('First orbit');
%Second orbit
subplot(2,2,2);plot(x2(:,1),x2(:,5));
xlabel('x(km)');ylabel('z(km)');
title('Second orbit');
%Third orbit
subplot(2,2,3);plot(x3(:,1),x3(:,5));
xlabel('x(km)');ylabel('z(km)');
title('Third orbit');
axis([-1 1 -0.25 2]*1e4);
%Fourth orbit
subplot(2,2,4);plot(x4(:,1),x4(:,5));
xlabel('x(km)');ylabel('z(km)');
title('Fourth orbit');
axis([-1 1 -0.25 2]*1e4);
suptitle('2-D Plots in the x-z plane');
%3D Plots
figure(2);
%First orbit
subplot(2,2,1);plot3(x1(:,1),x1(:,3),x1(:,5));
xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('First orbit');
%Second orbit
subplot(2,2,2);plot3(x2(:,1),x2(:,3),x2(:,5));
xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('Second orbit');
%Third orbit
subplot(2,2,3);plot3(x3(:,1),x3(:,5),x3(:,3));
xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('Third orbit');
axis([-1 1 -1 1 -0.25 2]*1e4);
%Fourth orbit
subplot(2,2,4);plot3(x4(:,1),x4(:,5),x4(:,3));
xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('Fourth orbit');
axis([-1 1 -1 1 -0.25 2]*1e4);
suptitle('3-D Plots');
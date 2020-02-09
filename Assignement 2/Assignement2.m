%Orbital Mechanics assignement
%Julien Huynh 14/11/2019 
clear all; clc;
global mu;
mu = 398600;
rTol = 1e-9; %Relative tolerance
aTol = 1e-9; %Absolute tolerance
t = [0:10:86400]; %Timespam
maxIt = 500; % Maximum number iteration to avoid running into an inf loop
%First orbit
r0 = [7078; 0; 0]; v0 = [0; 0; 7.5044]; 
[a, nE, i, omega, Omega, theta0] = OrbElements(r0, v0);
M0 = 2 * atan(sqrt((1 - nE) / (1 + nE)) * tan(theta0 / 2)) - nE * (sqrt(1 - nE ^ 2) * sin(theta0)) / (1 + nE * cos(theta0));
nMM = sqrt(mu / a ^ 3); % Mean Motion

%Second orbit
r2 = [7078; 0; 0]; v2 = [0; 0; 9.5044]; 
[a2, nE2, i2, omega2, Omega2, theta02] = OrbElements(r2, v2);
M2 = 2 * atan(sqrt((1 - nE2) / (1 + nE2)) * tan(theta02 / 2)) - nE2 * (sqrt(1 - nE2 ^ 2) * sin(theta02)) / (1 + nE2 * cos(theta02));
nMM2 = sqrt(mu / a2 ^ 3); % Mean Motion

% Initialization of first orbit
theta = theta0;
M = M0;
E = M + nE / 2; % M0 = 0
rc = a * (1 - nE * cos(E(1)));
rx = rc * cos(theta);
ry = rc * sin(theta);
rz = 0;

vx = (sqrt(mu * a) / rc) * (- sin(E));
vy = (sqrt(mu * a) / rc) * (sqrt(1 - nE ^ 2) * cos(E));
vz = 0;

% Initialization of second orbit
theta2 = theta02;
E2 = M2 + nE2 / 2; % M2 = 0
rc2 = a2 * (1 - nE2 * cos(E2));
rx2 = rc2 * cos(theta2);
ry2 = rc2 * sin(theta2);
rz2 = 0;

vx2 = (sqrt(mu * a2) / rc2) * (- sin(E2));
vy2 = (sqrt(mu * a2) / rc2) * (sqrt(1 - nE2 ^ 2) * cos(E2));
vz2 = 0;

%Compute position in inertial frame
[rxi1, ryi1, rzi1] = toInertial(i, omega, Omega, rx, ry); %First orbit
[rxi2, ryi2, rzi2] = toInertial(i2, omega2, Omega2, rx2, ry2); %Second orbit

%Compute velocity in inertial frame
[vxi1, vyi1, vzi1] = toInertial(i, omega, Omega, vx, vy); %First orbit
[vxi2, vyi2, vzi2] = toInertial(i2, omega2, Omega2, vx2, vy2); %Second orbit


%Looping
for j = 2:length(t)
   cpt = 0; % Count how many iterations are needed 
   
   %Compute Mean Anomaly
   M = horzcat(M, nMM * t(j)); %First orbit
   M2 = horzcat(M2, nMM2 * t(j)); %Second orbit
   
   %Compute Eccentric anomalies  
   Et = eccentricAnomaly(M, nE, j);
   Et2 = eccentricAnomaly(M2, nE2, j);
   
   while 1
       %Compute E_{j+1}
       Etp = Et - (Et - nE * sin(Et) - M(j)) / (1 + nE * cos(Et));
       Etp2 = Et2 - (Et2 - nE2 * sin(Et2) - M2(j)) / (1 + nE2 * cos(Et2));
       %Reassign E_{j+1} to E_j
       Et = Etp;
       Et2 = Etp2;
       %Increment counter
       cpt = cpt + 1;
       %Test condition
       if (abs(Et - Etp) / min(abs(Et), abs(Etp)) < rTol) && (abs(Et - Etp) < aTol) || (cpt >= maxIt)
           break;
       end
   end
   %Concantenate eccentric anomalies
   E = horzcat(E, Et);
   E2 = horzcat(E2, Et2);
   
   %Concatenate true anomaly
   theta = horzcat(theta, 2 * atan(sqrt((1 + nE) / (1 - nE)) * tan(Et / 2)));
   theta2 = horzcat(theta2, 2 * atan(sqrt((1 + nE2) / (1 - nE2)) * tan(Et2 / 2)));
   
   %Compute distance
   rc = horzcat(rc, a * (1 - nE * cos(E(j))));
   rc2 = horzcat(rc2, a2 * (1 - nE2 * cos(E2(j))));
   
   %Compute pos non inertial
   rx = horzcat(rx, rc(j) * cos(theta(j)));
   ry = horzcat(ry, rc(j) * sin(theta(j)));
   rz = horzcat(rz, 0);
   
   rx2 = horzcat(rx2, rc2(j) * cos(theta2(j)));
   ry2 = horzcat(ry2, rc2(j) * sin(theta2(j)));
   rz2 = horzcat(rz2, 0);
   
   % Compute velocity
   
   vx = horzcat(vx, (sqrt(mu * a) / rc(j)) * (- sin(E(j))));
   vy = horzcat(vy, (sqrt(mu * a) / rc(j)) * (sqrt(1 - nE ^ 2) * cos(E(j))));
   vz = horzcat(vz, 0);
   
   vx2 = horzcat(vx2,(sqrt(mu * a2) / rc2(j)) * (- sin(E2(j))));
   vy2 = horzcat(vy2, (sqrt(mu * a2) / rc2(j)) * (sqrt(1 - nE2 ^ 2) * cos(E2(j))));
   vz2 = horzcat(vz2, 0);
   
   
   % Inertial frame
   % Position
   [x, y, z] = toInertial(i, omega, Omega, rx(j), ry(j));
   rxi1 = horzcat(rxi1, x);
   ryi1 = horzcat(ryi1, y);
   rzi1 = horzcat(rzi1, z);
   
   [x2, y2, z2] = toInertial(i2, omega2, Omega2, rx2(j), ry2(j));
   rxi2 = horzcat(rxi2, x2);
   ryi2 = horzcat(ryi2, y2);
   rzi2 = horzcat(rzi2, z2);
   
   % Velocity
   [velx, vely, velz] = toInertial(i, omega, Omega, vx(j), vy(j));
   vxi1 = horzcat(vxi1, velx);
   vyi1 = horzcat(vyi1, vely);
   vzi1 = horzcat(vzi1, velz);
   
   [velx2, vely2, velz2] = toInertial(i2, omega2, Omega2, vx2(j), vy2(j));
   vxi2 = horzcat(vxi2, velx2);
   vyi2 = horzcat(vyi2, vely2);
   vzi2 = horzcat(vzi2, velz2);
   
   
end
figure(1)
plot3(rxi1, ryi1, rzi1);xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
axis([-1 1 -1 1 -1 1]*1e4);
title('First orbit');
figure(2)
plot3(rxi2, ryi2, rzi2);xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
axis([-1 1 -1 1 -1 1]*3e4);
title('Second orbit');
figure(3)
plot3(vxi1, vyi1, vzi1);xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('First orbit - Velocity');
figure(4)
plot3(vxi2, vyi2, vzi2);xlabel('x(km)');ylabel('y(km)');zlabel('z(km)');
title('Second orbit - Velocity');
% It seems like there is a problem with both the position and the velocity
% of the second orbit but I could not find why, it might be that one of the
% horzcat is concatenating as one point on the x axis is going to the
% opposite of the value it should be. However, the orbits are polar as
% expected since the inclination is pi.

%Orbital Mechanics assignement
%Julien Huynh 05/12/2019 
clear all; clc;
mu = 398600;
rE = 6378;


%Initial Orbit 
a1 = 29500;
e1 = 0.3;
b1 = a1 * sqrt( 1 - e1^2);
omega1 = 0;
p1 = a1 * (1 - e1^2);
%Final Orbit
a2 = 64000;
e2 = 0.5;
b2 = a2 * sqrt( 1 - e2^2);
omega2 = 0;
p2 = a2 * (1 - e2^2);


theta = 0: pi/180 : 2*pi;

%Velocities in first orbit
vr1 = sqrt(mu / p1) * e1 * sin(theta);
vtheta1 = sqrt(mu / p1) * (1 + e1 * cos(theta));
%[vir1, vit1] = toInertial2(omega1, vr1, vtheta1);


%Velocities in second orbit
vr2 = sqrt(mu / p2) * e2 * sin(theta);
vtheta2 = sqrt(mu / p2) * (1 + e2 * cos(theta));
%[vir2, vit2] = toInertial2(omega2, vr2, vtheta2);


%Transfer orbit
b = [- p1 ; - p2];
for j = 1:length(theta)

cpt = 0;
    %Anomaly in initial orbit
    for k = 1:length(theta)
       %Anomaly in final orbit      
       for l = 1:length(theta)
          %Anomaly in the transfer orbit at its cross with final orbit
          dTheta = omega2 - omega1 + theta(k) - theta(j);
          thetaT2 = dTheta + theta(l) ;
          r1 = p1 / (1+e1*cos(theta(j)));
          r2 = p2 / (1+e2*cos(theta(k)));
          eT = (r1 - r2)/(r2 * cos(thetaT2) - r1 * cos(theta(l)));
          if (eT<1) && (eT>0)
              aT = r1 * (1 + eT * cos(theta(l))) / (1 - eT ^2);
              pT = aT * (1-eT ^2);
              vinitT =  sqrt(mu / pT) * [eT * sin(theta(l)); (1 + eT * cos(theta(l)))];
              vfinT =  sqrt(mu / pT) *[eT * sin(thetaT2); (1 + eT * cos(thetaT2))];
              dV1 = [abs(vr1(j) - vinitT(1)); abs(vtheta1(j) - vinitT(2))];
              dV2 = [abs(vr2(j) - vfinT(1)); abs(vtheta2(j) - vfinT(2))];
              deltaV = dV1 + dV2;
              dV = norm(deltaV);
              if (cpt == 0) %Initialize
                      cpt = cpt+1;
                      t1 = theta(1);
                      tt1 = theta(1);
                      tt2 = thetaT2;
                      t2 = theta(1);
                      mindV = dV;
                      sDV = deltaV;
              end
              if (dV < mindV)
                          t1 = theta(j);
                          t2 = theta(k);
                          tt1 = theta(l);
                          tt2 = thetaT2;
                          mindV = dV; %Minimal delta v
                          sDV = deltaV;
                          seT = eT;
                          saT = aT;
                          spT = pT;
              end
          end
       end
    end
end
     
spT / (1 + seT * cos(tt1)) - p1 / (1+e1*cos(t1)) %Checking radius equality
spT / (1 + seT * cos(tt2)) - p2 / (1+e2*cos(t2)) %Checking radius equality
sbT = saT * sqrt( 1 - seT^2); 

omegaT = omega1 + t1 - tt1;
%Rotation matrices
R = [cos(omegaT) -sin(omegaT); ...
      sin(omegaT)  cos(omegaT)];
R1 = [cos(omega1) -sin(omega1); ...
      sin(omega1)  cos(omega1)];
R2 = [cos(omega2) -sin(omega2); ...
      sin(omega2)  cos(omega2)];
  
%Getting coordinates with respect to apse line

xt = -saT * seT + saT * cos(theta); %Transfer orbit
yt = sbT * sin(theta);%Transfer orbit
x1 = -a1 * e1 +a1 * cos(theta); %Orbit 1
y1 = b1 * sin(theta);%Orbit 1
x2 = -a2 * e2 +a2 * cos(theta);%Orbit 2
y2 = b2 * sin(theta);%Orbit 2
%Burn points
xB1 = -a1 * e1 +a1 * cos(t1); %Orbit 1
yB1 = b1 * sin(t1);%Orbit 1
xB2 = -saT * seT +saT * cos(tt2);%Orbit 2
yB2 = sbT * sin(tt2);%Orbit 2

%Getting them with respect to the default apse line at omega=0
rCT = R*[xt ; yt];   
xrt = rCT(1,:)';      
yrt = rCT(2,:)';
rC1 = R1*[x1 ; y1];   
xr1 = rC1(1,:)';      
yr1 = rC1(2,:)';
rC2 = R2*[x2 ; y2];   
xr2 = rC2(1,:)';      
yr2 = rC2(2,:)';

%Burns
rCB1 = R1*[xB1 ; yB1];   
xCb1 = rCB1(1,:)';      
yCb1 = rCB1(2,:)';

rCB2 = R*[xB2 ; yB2];   
xCb2 = rCB2(1,:)';      
yCb2 = rCB2(2,:)';

figure();
plot(xr1,  yr1);
hold on

scatter(xCb1, yCb1);
scatter(xCb2, yCb2);
plot(xr2,  yr2);
axis([-1 1 -1 1]*1e5);
plot(xrt,  yrt);
dim = [.2 .5 .3 .3];
str = 'Total delta v is ' + string(mindV) + 'km/s';
anomalies = 'Theta1 = ' + string(t1) + ', Theta2 = ' + string(t2)+ ', Theta T1 = '+string(tt1) + ', Theta T2 = ' + string(tt2)
top = 'Transfer orbit parameters : aT = ' + string(saT) + ' eT = ' + string(seT) + ', omegaT = ' + string(omegaT)
annotation('textbox',dim,'String',str,'FitBoxToText','on');
hold off
legend('Transfer orbit','Initial orbit', 'Final orbit', 'First maneuver', 'Second Maneuver')

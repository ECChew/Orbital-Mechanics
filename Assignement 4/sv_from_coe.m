function [r, v] = sv_from_coe(coe)
global muS
h = coe(1);
e = coe(2);
RA = coe(3);
incl = coe(4);
w = coe(5);
TA = coe(6);
%...Equations 4.37 and 4.38 (rp and vp are column vectors):
rp = (h^2/muS) * (1/(1 + e*cos(TA))) * (cos(TA)*[1;0;0] ...
+ sin(TA)*[0;1;0]);
vp = (muS/h) * (-sin(TA)*[1;0;0] + (e + cos(TA))*[0;1;0]);
%...Equation 4.39:
R3_W = [ cos(RA) sin(RA) 0
-sin(RA) cos(RA) 0
0 0 1];
%...Equation 4.40:
R1_i = [1 0 0
0 cos(incl) sin(incl)
0 -sin(incl) cos(incl)];
%...Equation 4.41:
R3_w = [ cos(w) sin(w) 0
-sin(w) cos(w) 0
0 0 1];

Q_pX = R3_W'*R1_i'*R3_w';

r = (Q_pX*rp)';
v = (Q_pX*vp)';


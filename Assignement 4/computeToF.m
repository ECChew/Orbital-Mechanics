function [eE, eM, mE, mM, tofCalc] = computeToF(eT, aT, theta0, alpha)
global muS
%Compute Eccentric anomalies in transfer orbit
%At departure
eE = atan2(sqrt(1 - eT ^2) * sin(theta0)/(1 + eT * cos(theta0)), (eT + cos(theta0))/(1 + eT * cos(theta0)));
%At arrival
eM = atan2(sqrt(1 - eT ^2) * sin(theta0 + alpha)/(1 + eT * cos(theta0 + alpha)), (eT + cos(theta0 + alpha))/(1 + eT * cos(theta0 + alpha)));
%Compute corresponding mean anomalies in transfer orbit
%At departure
mE = eE - eT * sin(eE);
%At arrival
mM = eM - eT * sin(eM);
%Compute corresponding time of flight
if (mM - mE) < 0
    tofM = sqrt(aT ^ 3 / muS) * mM;
    tofE = sqrt(aT ^ 3 / muS) * mE;
    T = sqrt(aT ^ 3 / muS) * 2 * pi;
    tofCalc = tofM + (T - tofE); 
else
    tofCalc = sqrt(aT ^ 3 / muS) *(mM-mE); 
end
end


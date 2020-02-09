%Orbital Mechanics assignement
%Julien Huynh & Clément Chandon 29/01/2020
clearvars; clc;
global muS;
muS = 1.327124e11;
muE = 398600;

%Creating date times
t1 = datetime(2005,6,20,12,0,0);
t2 = datetime(2005,11,10,12,0,0);
dd = t1:t2; % Departure date
t1 = datetime(2006,1,1,12,0,0);
t2 = datetime(2007,3,1,12,0,0);
da = t1:t2; % Arrival date




 
%Initialize distance and anomaly table
[rE, LE] = ephem(2005, 6, 20, 12, 0, 0, 3);
[rM, LM] = ephem(2006, 1, 1, 12, 0, 0, 4);
[iE, wE, thetaE] = marsAngles(2006, 1, 1);
[iM, wM, thetaM] = marsAngles(2006, 1, 1);

for j = 2:length(dd)
    [rEtemp, LEtemp] = ephem(dd(j).Year, dd(j).Month, dd(j).Day, dd(j).Hour, dd(j).Minute, dd(j).Second, 3);
    rE = horzcat(rE, rEtemp);
    LE = horzcat(LE, LEtemp);
    [iEtemp, wEtemp, thetaEtemp] = marsAngles(da(j).Year, da(j).Month, da(j).Day);
    iE = horzcat(iE, iEtemp);
    wE = horzcat(wE, wEtemp);
    thetaE = horzcat(thetaE, thetaEtemp);
end
%Getting distances for each arrival date
for j = 2:length(da)
    [rMtemp, LMtemp] = ephem(da(j).Year, da(j).Month, da(j).Day, da(j).Hour, da(j).Minute, da(j).Second, 4);
    rM = horzcat(rM, rMtemp);
    LM = horzcat(LM, LMtemp);
    [iMtemp, wMtemp, thetaMtemp] = marsAngles(da(j).Year, da(j).Month, da(j).Day);
    iM = horzcat(iM, iMtemp);
    wM = horzcat(wM, wMtemp);
    thetaM = horzcat(thetaM, thetaMtemp);
end

LE = LE .* pi / 180;
LM = LM .* pi / 180;
iM = iM .* pi / 180;
wM = wM .* pi / 180;
thetaM = thetaM .* pi / 180;

counter = 0;
for j = 1:length(dd)
    for k = 1:length(da)
        tof = seconds(da(k) - dd(j));
        alpha = mod(LM(k) - LE(j), 2*pi);
        theta0 = 0;
        iT = asin(sin(wM(k) + thetaM(k)) * sin(iM(k))/sin(alpha));
        cpt = 0;
        while(theta0<2*pi)
            eT = (rM(k) - rE(j)) / (rE(j) * cos(theta0) - rM(k) * cos(theta0 + alpha)); %Eccentricity
            if (eT<0.5)&&(eT>0)
                pT = rE(j) * (1 + eT * cos(theta0)); %Parameter
                aT = pT / (1 - eT ^ 2); %Semi major axis
                [eE, eM, mE, mM, tofCalc] = computeToF(eT, aT, theta0, alpha);
                if cpt == 0
                   s = tofCalc/86400;
                   t = theta0;
                   diff = abs(tof/86400 - tofCalc/86400);
                   swT = thetaE(j) + wE(j) - theta0;
                   siT = iT;
                   seT = eT;
                   spT = pT; 
                   cpt = 1;
                else
                    t = horzcat(t, theta0);
                    s = horzcat(s, tofCalc/86400);
                    swT = horzcat(swT, thetaE(j) + wE(j) - theta0);
                    siT = horzcat(siT, iT);
                    seT = horzcat(seT, eT);
                    spT = horzcat(spT, pT);
                    diff = horzcat(diff, abs(tof/86400-tofCalc/86400));
                end
            else
                t = horzcat(t, NaN);
                diff = horzcat(diff, NaN);
                swT = horzcat(swT, NaN);
                siT = horzcat(siT, NaN);
                seT = horzcat(seT, NaN);
                spT = horzcat(spT, NaN);
            end
         
                theta0 = theta0+pi/180;
        end
        [m, i] = min(diff);
        if m <= 1
            theta = t(i);
            eT2 = seT(i);
            pT2 = spT(i);
            iT2 = siT(i);
            iT2 = siT(i);
            wT2 = swT(i);
            vTE = sqrt(muS / pT2) * [eT2 * sin(theta); (1 + eT2 * cos(theta)); 0];
            Rw = [cos(wT2) sin(wT2) 0;
                -sin(wT2) cos(wT2) 0;
                0 0 1];
            Ri = [1 0 0;
                0 cos(iT2) sin(iT2);
                0 -sin(iT2) cos(iT2)];
            Q = Rw * Ri;
            vTE = Q * vTE;
            gammaTE = atan2(eT2 * sin(theta), 1 + eT2 * cos(theta));
            [vE, gammaE] = EarthVandGamma3D(dd(j).Year,dd(j).Month , dd(j).Day);
            vinf = sqrt(norm(vE)^2 + norm(vTE)^2 - 2 * norm(vE) * norm(vTE) *cos(gammaTE - gammaE) * cos(iT2));
            if vinf < 7
                if counter == 0
                    date = [dd(j), da(k)];
                    nvinf = vinf;
                    counter = 1;
                else
                    nvinf = horzcat(nvinf, vinf);
                    date = vertcat(date,[dd(j), da(k)]);
                end
            end
        end
    end
end
scatter3(date(:, 1), date(:, 2), nvinf);
dlmwrite('dataincscatter.txt',nvinf);
xlabel('Departure date');ylabel('Arrival date');zlabel('vinf');

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
ddformat = datenum(dd);
t1 = datetime(2006,1,1,12,0,0);
t2 = datetime(2007,3,1,12,0,0);
da = t1:t2; % Arrival date
daformat = datenum(da);

ddmat = repmat(ddformat,size(daformat,2),1);
damat = repmat(daformat',1,size(ddformat,2));

 
%Initialize distance and anomaly table
[rE, LE] = ephem(2005, 6, 20, 12, 0, 0, 3);
[rM, LM] = ephem(2006, 1, 1, 12, 0, 0, 4);

for j = 2:length(dd)
    [rEtemp, LEtemp] = ephem(dd(j).Year, dd(j).Month, dd(j).Day, dd(j).Hour, dd(j).Minute, dd(j).Second, 3);
    rE = horzcat(rE, rEtemp);
    LE = horzcat(LE, LEtemp);
end
%Getting distances for each arrival date
for j = 2:length(da)
    [rMtemp, LMtemp] = ephem(da(j).Year, da(j).Month, da(j).Day, da(j).Hour, da(j).Minute, da(j).Second, 4);
    rM = horzcat(rM, rMtemp);
    LM = horzcat(LM, LMtemp);
end
LE = LE .* pi / 180;
LM = LM .* pi / 180;

counter = 0;
for j = 1:length(dd)
    for k = 1:length(da)
        tof = seconds(da(k) - dd(j));
        alpha = mod(LM(k) - LE(j), 2*pi);
        theta0 = 0;
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
                   seT = eT;
                   spT = pT;
                   cpt = 1;
                else
                    t = horzcat(t, theta0);
                    s = horzcat(s, tofCalc/86400);
                    seT = horzcat(seT, eT);
                    spT = horzcat(spT, pT);
                    diff = horzcat(diff, abs(tof/86400-tofCalc/86400));
                end
            else
                t = horzcat(t, NaN);
                diff = horzcat(diff, NaN);
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
            vTE = sqrt(muS / pT2) * [eT2 * sin(theta); (1 + eT2 * cos(theta))];
            gammaTE = atan2(eT2 * sin(theta), 1 + eT2 * cos(theta));
            [vE, gammaE] = EarthVandGamma(dd(j).Year,dd(j).Month , dd(j).Day);
            vinf = vTE - vE;
            if norm(vinf) < 7
                if counter == 0
                    %date = [dd(j), da(k)];
                    nvinf = norm(vinf);
                    counter = 1;
                else
                    %date = vertcat(date, [dd(j), da(k)]);
                    nvinf = horzcat(nvinf, norm(vinf));
                end
            else
                if counter == 0
                    nvinf = NaN;
                    counter = 1;
                else
                    nvinf = horzcat(nvinf, NaN);
                end
            end
        else
            if counter == 0
                    nvinf = NaN;
                    counter = 1;
            else
                    nvinf = horzcat(nvinf, NaN);
            end
        end
        if j == 1 && k == length(da)
            v = nvinf;
            counter = 0;
        elseif k == length(da)
            v = vertcat(v, nvinf);
            counter = 0; %Reset so that we can reinitialize the nvinf1 after concatenating
        end
    end
end
contourf(ddmat,damat,v');
dlmwrite('data.txt',v);
xlabel('Departure date');ylabel('Arrival date');zlabel('vinf');

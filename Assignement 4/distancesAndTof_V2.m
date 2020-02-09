function [A] = distancesAndTof(dd, md, yd, dlen, da, ma, ya, alen)

m31 = [1, 3, 5, 7, 8, 10];
DD = zeros(dlen, 3);
hour = 0;
minute = 0;
second = 0;

for i = 1:dlen
    DD(i, 1) = yd;
    DD(i, 2) = md;
    DD(i, 3) = dd;
    if (~ismember(md, m31)) && dd >= 30
        dd = 1;
        md = md + 1;
    elseif (ismember(md, m31)) && dd >= 31
        dd = 1;
        md = md + 1;
    else
        dd = dd + 1;
    end
end

DA = zeros(alen, 3);

for i = 1:alen
    DA(i, 1) = ya;
    DA(i, 2) = ma;
    DA(i, 3) = da;
    if (~ismember(ma, m31)) && da >= 30 && ma~=2 && ma~= 12
        da = 1;
        ma = ma + 1;
    elseif (ismember(ma, m31)) && da >= 31
        da = 1;
        ma = ma + 1;
    elseif (ma == 2 && da >= 28)
        da = 1;
        ma = ma + 1;
    elseif (ma == 12 && da >= 31)
        ya = ya + 1;
        ma = 1;
        da = 1;
    else
        da = da + 1; 
    end
end

A = zeros(dlen * alen, 15);
cnt = 1;
for i = 1:dlen
    year_d = DD(i, 1);
    month_d = DD(i, 2);
    day_d = DD(i, 3);
    
    [coe_Ed, r_Ed, vE, ~] = planet_elements_and_sv ...
                    (3, year_d, month_d, day_d, hour, minute, second);
    [~, r_Md, ~, ~] = planet_elements_and_sv ...
                    (4, year_d, month_d, day_d, hour, minute, second);
%     alpha_d = coe_Ed(6)-coe_Md(6);

    for j = 1:alen
        year_a = DA(j, 1);
        month_a = DA(j, 2);
        day_a = DA(j, 3);
        
        [~, r_Ea, ~, ~] = planet_elements_and_sv ...
                    (3, year_a, month_a, day_a, hour, minute, second);
        [coe_Ma, r_Ma, ~, ~] = planet_elements_and_sv ...
                        (4, year_a, month_a, day_a, hour, minute, second);
        alpha = coe_Ed(6)-coe_Ma(6);
        
        A(cnt, 1) = datenum(year_d, month_d, day_d);
        A(cnt, 2) = datenum(year_a, month_a, day_a);
        A(cnt, 3) = norm(r_Ed);
        A(cnt, 4) = alpha*180/pi;
        A(cnt, 5) = norm(r_Ma);
%         A(cnt, 8) = alpha;
        A(cnt, 6) = norm(vE);
        A(cnt, 7) = coe_Ed(2); % eccentricity of Earth at departure
        A(cnt, 8) = coe_Ed(6)*180/pi; % true anomaly of Earth at departure
        A(cnt, 9) = coe_Ma(6)*180/pi; % true anomaly of Mars at arrival
        A(cnt, 10) = coe_Ma(5); % argument of perihelion of Mars at arrival
        A(cnt, 11) = coe_Ma(4); % inclination of Mars at arrival
        
        cnt = cnt + 1;
    end
end
end

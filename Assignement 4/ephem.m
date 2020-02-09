function [nr, L] = ephem(y, m, d, h, min, s, planet_id)
[coe, r, v, jd] = planet_elements_and_sv(planet_id, y, m, d, h, min, s);
nr = norm(r); %Distance
L = mod(coe(5) + coe(6) + coe(3), 360); %Absolute anomaly
end


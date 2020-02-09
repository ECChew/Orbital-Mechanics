function [inc, w, theta] = EarthAnomalies(y, m, d)
[coe, r, v, jd] = planet_elements_and_sv(3, y, m, d, 12, 0, 0);
inc = coe(4); %Inclination in degrees
w = coe(5); %Argument of perihelion
theta = coe(6); %True anomaly in degrees
end


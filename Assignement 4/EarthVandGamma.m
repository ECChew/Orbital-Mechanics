function [vE, gammaE] = EarthVandGamma(y, m, d)
[coe, r, v, jd] = planet_elements_and_sv(3, y, m, d,12, 0, 0);
vE = [v(1); v(2)];
e = coe(2); 
theta = coe(6);
gammaE = atan2(e * sin(theta), 1 + e * cos(theta));
end
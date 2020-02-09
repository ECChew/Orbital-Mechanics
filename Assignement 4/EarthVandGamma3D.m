function [vE, gammaE] = EarthVandGamma3D(y, m, d)
[coe, r, vE, jd] = planet_elements_and_sv(3, y, m, d,12, 0, 0);
e = coe(2); 
theta = coe(6);
gammaE = atan2(e * sin(theta), 1 + e * cos(theta));
end
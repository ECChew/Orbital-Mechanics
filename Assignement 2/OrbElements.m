function [a, nE, i, omega, Omega, theta0] = OrbElements(r, v)
global mu;
h = cross(r, v);
nR = norm(r); nV = norm(v) ;% Amplitudes of r0 and v0
a = (2 /nR - nV ^ 2 / mu) ^(-1);
e = 1 / mu * (cross(v, h) - mu * r ./ nR);
nE = norm(e);
i = acos(h(3) / norm(h));
k = [0; 0; 1];
n = cross(k,h);
if (n(2) < 0)
    Omega = 2 * pi - arccos(n(1)/norm(n));
else
    Omega = acos(n(1)/norm(n));
end
if (e(3) < 0)
    omega = 2 * pi - acos(dot(n, e)/(norm(n)* norm(e)));
else
    omega = acos(dot(n, e)/(norm(n)* norm(e)));
end
if (dot(r, v) < 0)
    theta0 = 2 * pi - acos(dot(e, r)/(norm(e)* norm(r0)));
else
    theta0 = acos(dot(e, r)/(norm(e)* norm(r)));
end
end


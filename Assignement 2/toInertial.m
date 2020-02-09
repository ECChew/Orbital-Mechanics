function [rxi, ryi, rzi] = toInertial(i, omega, Omega, rx, ry)
rxi = rx * (cos(omega) * cos(Omega) - sin(omega) * cos(i) * sin(Omega)) - ry * (sin(omega) * cos(Omega) + cos(omega) * cos(i) * sin(Omega));
ryi = rx * (cos(omega) * sin(Omega) + sin(omega) * cos(i) * cos(Omega)) + ry * (cos(omega) * cos(i) * cos(Omega) - sin(omega) * sin(Omega));
rzi = rx * (sin(omega) * sin(i)) + ry * (cos(omega) * sin(i));
end


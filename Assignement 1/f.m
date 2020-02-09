function [Xdot] = f(t, X)
    mu = 398600; x = X(1); xd = X(2); y = X(3); yd = X(4); z = X(5); zd = X(6);
    r = (x^2+y^2+z^2)^0.5; xdd = -1*mu*x/(r^3); ydd = -1*mu*y/(r^3); zdd = -1*mu*z/(r^3);
    Xdot = [xd;xdd;yd;ydd;zd;zdd];
    end

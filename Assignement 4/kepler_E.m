function E = kepler_E(e, M)
tol = 1e-6;
%Initialize
if M < pi
    E = M + e/2;
else
    E = M - e/2;
end
%Solve with Newton iterative method
while 1
    if abs(E - e*sin(E) - M)/(1 - e*cos(E)) < tol
        break
    else
        E = E - (E - e*sin(E) - M)/(1 - e*cos(E));
    end
end
end

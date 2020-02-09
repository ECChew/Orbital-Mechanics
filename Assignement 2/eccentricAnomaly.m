function Et = eccentricAnomaly(M, nE, j)
if (M(j - 1) < pi) && (M(j - 1) > 0)
        Et = M(j) + nE / 2;
   else
        Et = M(j) - nE / 2;
   end
end


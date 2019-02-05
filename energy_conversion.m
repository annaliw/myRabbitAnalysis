function energy = energy_conversion(xvar, t0, param)
    energy = param(2)*1./(xvar-t0).^2 + param(1); 
end
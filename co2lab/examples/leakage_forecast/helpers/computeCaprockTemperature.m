function t = computeCaprockTemperature( Gt, seafloor_temp, seafloor_depth, temp_gradient )
% NB: returns value in Celsius, not Kelvin

    t = seafloor_temp + (Gt.cells.z - seafloor_depth) .* temp_gradient ./ 1000;
    
end


function [rho, mu, sr, sw] = getValuesSPE134891()

    % Following values from p. 13-14 of SPE 134891 ("Reservoir modeling of CO2
    % plume behavior calibrated against monitoring data from Sleipner, Norway").

    rho = [1020 760] .* kilogram/meter^3;  % [brine, CO2]
    mu  = [8e-4 6e-2*milli] * Pascal*second; % [brine, CO2]
    sr  = 0.21; % Residual gas saturation 
    sw  = 0.11; % Residual brine (oil) saturation

end

function h = convertMassesToHeight(masses, CO2obj, p, tinfo, rhoW, slope, ...
                                   z, areas, poro, approximate, fully_compressible)
    N = numel(masses);
    h = zeros(N,1);
    
    % compute caprock temperature
    T = tinfo{1} + (z - tinfo{2}) .* (tinfo{3}/1000);

    for i = 1:N
        if masses(i)~= 0
            h(i) = columnHeight(p(i), T(i), tinfo{3}/1000, CO2obj, rhoW, slope, masses(i), ...
                                areas(i), poro(i), approximate, fully_compressible);
        else
            h(i) = 0;
        end
    end
end

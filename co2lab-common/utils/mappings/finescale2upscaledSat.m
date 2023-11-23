function [S, Smax] = finescale2upscaledSat(sg, res_wat, res_gas, Gt, poro)

    cix = Gt.columns.cells;
    pvol = poro(cix) .* Gt.parent.cells.volumes(cix);
    gvol = pvol .* sg(cix);
    
    cmap = rldecode((1:Gt.cells.num)', diff(Gt.cells.columnPos));
    col_pvol = accumarray(cmap, pvol);
    col_gvol = accumarray(cmap, gvol);
    
    S = col_gvol ./ col_pvol; % volume of gas in a column as a fraction of
                              % the total porevolume of the column
    
    % to determine Smax, we fill all cells with gas saturation >= res_gas 
    % completely with CO2, and compute S again
    threshold = 0.9; % allow residual saturation to be slightly lower than res_gas
    sgmax = sg;
    % full saturation in all cells with saturation about residual threshold
    sgmax(sgmax > threshold * res_gas) = 1 - res_wat; 
    
    gmaxvol = pvol .* sgmax(cix);
    col_gmaxvol = accumarray(cmap, gmaxvol);
    Smax = col_gmaxvol ./ col_pvol;
    
end


sGmax or smax
s
pressure
h

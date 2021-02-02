function sG = getGasSaturationFromHeight(T, t, B, b, h)
    % Internal function: Reconstruct gas saturation from VE assumption
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.
    d_local = t + h;
    sG = (d_local - T)./(B-T);
    sG = max(sG, 0); % Gas < 0 means that the gas-liquid interface is actually above the cell top
    sG = min(sG, 1); % Gas > 1 means that the gas-liquid interface is below the cell bottom
end

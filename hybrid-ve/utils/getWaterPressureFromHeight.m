function pW_local = getWaterPressureFromHeight(pt, t, B, b, h, pW, g, rhow, rhog)
    % Internal function: Reconstruct water pressure from VE assumption
    % T is top of column for each cell, t is top of cells to be
    % reconstructed, with b/B having the same interpretation for the
    % respective cell bottoms. h is the height of the gas plume in each
    % cell.
    % pW is the water pressure at the bottom of the cell, g the gravity
    % magnitude and rhow/rhog the water and gas densities.

    dz_below = (b - pt);

    h_below = max(t + h - pt, 0);
    pW_local = pW - ((dz_below - h_below).*rhow + h_below.*rhog).*g;
end
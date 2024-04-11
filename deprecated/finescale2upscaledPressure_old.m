function P = finescale2upscaledPressure(p, Gt, fluid)
 
    % index of uppermost cell in each column
    cells_upper = Gt.columns.cells(Gt.cells.columnPos(1:end-1));
    
    p_upper = p(cells_upper);
    
    dz_upper = Gt.columns.dz(cells_upper);
    
    rho = fluid.rhoWS .* fluid.bW(p_upper);
    
    P = p_upper - rho .* dz_upper/2 * norm(gravity());
    
end

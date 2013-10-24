function [alpha, intAlpha, intAlpha2] = setupCorrectionPolys(fluid, p, t, theta)

    % if strcmpi(fluid.compressible, 'FULL')
    %     % We use a polynomial expansion to correct for density changes in
    %     % vertical direction
    %     VEP = VEpolys;
    %     rho_top = fluid.rho  (p, t);
    %     beta_top = fluid.beta(p, t);
    %     if isa(p, 'ADI')
    %         % avoid making bder_top an ADI variable, as the 'CO2props' object
    %         % cannot supply a matching derivative function
    %         bder_top = fluid.bder(p.val, t);
    %     else
    %         bder_top = fluid.bder(p, t);
    %     end

    %     alpha     = @(h) VEP.alpha     (h, rho_top, beta_top, bder_top, theta);
    %     intAlpha  = @(h) VEP.intAlpha  (h, rho_top, beta_top, bder_top, theta);
    %     intAlpha2 = @(h) VEP.intSqAlpha(h, rho_top, beta_top, bder_top, theta);
    % else
        % No correction for density changes should be done in the vertical
        % direction.  
        alpha     = @(h) 1;
        intAlpha  = @(h) h;
        intAlpha2 = @(h) h; 
    % end

end

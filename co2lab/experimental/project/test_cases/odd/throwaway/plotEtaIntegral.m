function h = plotEtaIntegral(p_range, t_range, CO2, theta, G, height, choice)
%
% Choice: 1 - Ieta; 2 - INupEta; 3 - INugEta; 4 - Ieta2

    res = 300;
    
    [h, pvals, tvals, pgrid, tgrid] = rhoColorPlot(CO2, p_range, t_range, res, ...
                                                   res, 20);
    
    gct = norm(gravity) * cos(theta);
    Gct = G/1000 * cos(theta);

    EOS.rho    = @CO2.rho;
    EOS.beta   = @CO2.beta;     
    EOS.gamma  = @CO2.gamma;        
    EOS.beta2  = @(p,t) CO2.rhoDPP(p,t) ./ CO2.rho(p,t);
    EOS.gamma2 = @(p,t) CO2.rhoDTT(p,t) ./ CO2.rho(p,t);        
    EOS.chi    = @(p,t) CO2.rhoDPT(p,t) ./ CO2.rho(p,t);
    EOS.compressible = 'full';
    
    [Ieta, INupEta, INugEta, Ieta2, Eta] = etaIntegrals(EOS, pgrid(:), tgrid(:), Gct, gct);

    switch choice
      case 1
        fun = Ieta;
      case 2
        fun = INupEta;
      case 3
        fun = INugEta;
      case 4
        fun = Ieta2;
      case 5
        fun = Eta;
      otherwise
        fun = Ieta;
    end
    
    %vals = reshape(fun(height*ones(res*res, 1)), res, res) - 1;
    vals = reshape(fun(height*ones(res*res, 1)), res, res);
    
    surf(vals, 'edgecolor','none');
    %cvals = [-0.01, 0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.11, 0.2, 0.5, 1.0];
    %cvals = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02];
    %cvals = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0];
    % cvals = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, ...
    %          0.04, 0.05, 0.1, 0.2, 0.5, 1.0];
    cvals = [-0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0, 0.01, 0.02, 0.03, ...
             0.04, 0.05, 0.1, 0.2, 0.3, 0.4];
    cvals = cvals+1;
    
    [c, h] = contour(tvals, pvals, vals, cvals, 'k');
    clabel(c, h, 'FontSize', 15);
    hold off;
    
end

% hold on;
    % [c, h] = contour(gvals, dvals, m, cvals, 'k');
    % clabel(c, h);

function compareEtas

    gravity on;
    g = norm(gravity);
    h = 100;
    G = 45; % temperature gradient
    RES = 300;
            
% Loading CO2
    CO2 = CO2props('rho_big_trunc', '');
    
% define test range
    pspan = [5 9] * 1e6;
    tspan = [300 315];
    
% define grid
   pvals = linspace(pspan(1), pspan(2), RES);
   tvals = linspace(tspan(1), tspan(2), RES);
   
   [tgrid, pgrid] = meshgrid(tvals, pvals);
   
   eta1 = eta(g, CO2, pgrid(:), tgrid(:), G/1000, h, 2);
   fp_eta1 = fp_eta(g, CO2, pgrid(:), tgrid(:), G/1000, h, 2);

% define alternative eta function

   EOS.rho          = @CO2.rho;
   EOS.beta         = @CO2.beta;
   EOS.gamma        = @CO2.gamma; 
   EOS.beta2        = @(p,t) CO2.rhoDPP(p,t)./CO2.rho(p,t);
   EOS.gamma2       = @(p,t) CO2.rhoDTT(p,t)./CO2.rho(p,t);
   EOS.chi          = @(p,t) CO2.rhoDPT(p,t)./CO2.rho(p,t);
   EOS.compressible = 'full';
  
   [~, ~, ~, ~, etafun, fpfun, ~, fpetafun] = etaIntegrals(EOS, pgrid(:), tgrid(:), G, g);

   eta2 = etafun(h);
   fp2  = fpfun(h);
   fp_eta2 = fpetafun(h);
   
   eta1 = reshape(eta1, RES, RES);
   eta2 = reshape(eta2, RES, RES);
   
   surf(eta1-eta2, 'edgecolor','none');
   keyboard;
    
end


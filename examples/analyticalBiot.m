function [u, f] = analyticalBiot(d, params)
    
% d = spatial dimension

    mu = params.mu;
    lamb = params.lamb;
    alpha = params.alpha;
    K = params.K;
    tau = params.tau;
    rho = params.rho; 

    % displacement
    u = cell(d, 1);
    % force
    f = cell(d, 1);

    switch d 
      case 2
        u{1} = @(x, y) (y.*(1 - x).*sin(2.*pi.*x.*y));
        u{2} = @(x, y) (y.^2.*cos(2.*pi.*x));
        u{3} = @(x, y) (u{1}(x, y));
        f{1} = @(x, y) (-alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) - 4.*pi.*(lamb.*y.*(-pi.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x)) - 2.*mu.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) + mu.*(-pi.*x.^2.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*(x - 1).*cos(2.*pi.*x.*y) + y.*sin(2.*pi.*x))));
        f{2} = @(x, y) (-alpha.*(x - 1).*(2.*pi.*x.*y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) - 2.*pi.*lamb.*y.*(-2.*pi.*x.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*cos(2.*pi.*x.*y) + (x - 1).*cos(2.*pi.*x.*y)) - lamb.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y) - 2.*cos(2.*pi.*x)) - mu.*(-4.*pi.^2.*x.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + 2.*pi.*x.*y.*cos(2.*pi.*x.*y) + 4.*pi.^2.*y.^2.*cos(2.*pi.*x) + 4.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 4.*mu.*cos(2.*pi.*x));
        f{3} = @(x, y) (4.*pi.*K.*tau.*x.*(x - 1).*(pi.*x.*y.*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) + 4.*pi*K*tau*y.^2*(pi*y*(x - 1)*sin(2*pi*x*y) - cos(2*pi*x*y)) - alpha*y*(2*pi*y*(x - 1)*cos(2*pi*x*y) + sin(2*pi*x*y) - 2*cos(2*pi*x)) - rho*y*(x - 1)*sin(2*pi*x*y));
      case 3
        
      otherwise
        error('d not recognized')

    end
    
    
    
end


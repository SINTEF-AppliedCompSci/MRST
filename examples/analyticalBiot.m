function [u, p, force, src] = analyticalBiot(d, params)
    
    % d = spatial dimension
    mu = params.mu;
    lambda = params.lambda;
    alpha = params.alpha;
    K = params.K;
    tau = params.tau;
    rho = params.rho; 

    % Displacement
    u = cell(d, 1);
    % Volumetric force
    force = cell(d, 1);

    switch d 
      case 2
        u{1} = @(x, y) (y.*(1 - x).*sin(2.*pi.*x.*y));
        u{2} = @(x, y) (y.^2.*cos(2.*pi.*x));
        p = @(x, y) (u{1}(x, y));
        force{1} = @(x, y) (-alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 4.*pi.*(lambda.*y.*(-pi.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x)) - 2.*mu.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) + mu.*(-pi.*x.^2.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*(x - 1).*cos(2.*pi.*x.*y) + y.*sin(2.*pi.*x))));
        force{2} = @(x, y) (-alpha.*(x - 1).*(2.*pi.*x.*y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 2.*pi.*lambda.*y.*(-2.*pi.*x.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*cos(2.*pi.*x.*y) + (x - 1).*cos(2.*pi.*x.*y)) + lambda.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y) - 2.*cos(2.*pi.*x)) + mu.*(-4.*pi.^2.*x.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + 2.*pi.*x.*y.*cos(2.*pi.*x.*y) + 4.*pi.^2.*y.^2.*cos(2.*pi.*x) + 4.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) - 4.*mu.*cos(2.*pi.*x));
        src = @(x, y) (-4.*pi.*K.*tau.*x.*(x - 1).*(pi.*x.*y.*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - 4.*pi.*K.*tau.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y) - 2.*cos(2.*pi.*x)) - rho.*y.*(x - 1).*sin(2.*pi.*x.*y));
      case 3
        error('not yet implemented');
      otherwise
        error('d not recognized');

    end
    
    
    
end


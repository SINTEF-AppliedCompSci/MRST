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
    % Stress
    stress = cell(d*(d + 1)/2, 1);    
    switch d 
      case 2
        % Displacement
        u{1} = @(x, y) (y.*(1 - x).*sin(2.*pi.*x.*y));
        u{2} = @(x, y) (y.^2.*cos(2.*pi.*x));
        % Pressure
        p = @(x, y) (u{1}(x, y));
        % Volumetric mechanical force
        force{1} = @(x, y) (-alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 4.*pi.*(lambda.*y.*(-pi.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x)) - 2.*mu.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) + mu.*(-pi.*x.^2.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*(x - 1).*cos(2.*pi.*x.*y) + y.*sin(2.*pi.*x))));
        force{2} = @(x, y) (-alpha.*(x - 1).*(2.*pi.*x.*y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 2.*pi.*lambda.*y.*(-2.*pi.*x.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*cos(2.*pi.*x.*y) + (x - 1).*cos(2.*pi.*x.*y)) + lambda.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y) - 2.*cos(2.*pi.*x)) + mu.*(-4.*pi.^2.*x.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + 2.*pi.*x.*y.*cos(2.*pi.*x.*y) + 4.*pi.^2.*y.^2.*cos(2.*pi.*x) + 4.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) - 4.*mu.*cos(2.*pi.*x));
        % Fluid source term
        src = @(x, y) (-4.*pi.*K.*tau.*x.*(x - 1).*(pi.*x.*y.*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - 4.*pi.*K.*tau.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y) - 2.*cos(2.*pi.*x)) - rho.*y.*(x - 1).*sin(2.*pi.*x.*y));
        % Stress (we do not follow here the convention that the third should be multiplied by two)
        %
        % stress tensor =  | stress{1}, stress{3} |  
        %                  | stress{3}, stress{2} |
        %
        % first Voigt component
        stress{1} = lambda.*y.*(-2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) - sin(2.*pi.*x.*y) + 2.*cos(2.*pi.*x)) + 1.0.*mu.*(4.*pi.*y.^2.*(1 - x).*cos(2.*pi.*x.*y) - 2.*y.*sin(2.*pi.*x.*y));
        % second Voigt component
        stress{2} = lambda.*y.*(-2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) - sin(2.*pi.*x.*y) + 2.*cos(2.*pi.*x)) + 4.0.*mu.*y.*cos(2.*pi.*x);
        % third  Voigt component
        stress{3} = 1.0.*mu.*(2.*pi.*x.*y.*(1 - x).*cos(2.*pi.*x.*y) - 2.*pi.*y.^2.*sin(2.*pi.*x) + (1 - x).*sin(2.*pi.*x.*y));
      case 3
        error('not yet implemented');
      otherwise
        error('d not recognized');

    end
    
    
    
end


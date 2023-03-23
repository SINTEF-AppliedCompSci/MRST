function output = analyticalBiot(d, params)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


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
        % Voigt components
        stress{1} = @(x, y) (lambda.*y.*(-2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) - sin(2.*pi.*x.*y) + 2.*cos(2.*pi.*x)) + 1.0.*mu.*(4.*pi.*y.^2.*(1 - x).*cos(2.*pi.*x.*y) - 2.*y.*sin(2.*pi.*x.*y)));
        stress{2} = @(x, y) (lambda.*y.*(-2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) - sin(2.*pi.*x.*y) + 2.*cos(2.*pi.*x)) + 4.0.*mu.*y.*cos(2.*pi.*x));
        stress{3} = @(x, y) (1.0.*mu.*(2.*pi.*x.*y.*(1 - x).*cos(2.*pi.*x.*y) - 2.*pi.*y.^2.*sin(2.*pi.*x) + (1 - x).*sin(2.*pi.*x.*y)));

      case 3
        % Displacement
        u{1} = @(x, y, z) (y.*(1 - x).*sin(2.*pi.*x.*y));
        u{2} = @(x, y, z) (z.*y.^2.*cos(2.*pi.*x));
        u{3} = @(x, y, z) (x.*y.*z);

        % Pressure
        p = @(x, y, z) (u{1}(x, y, z));

        % Volumetric mechanical force
        force{1} = @(x, y, z) (-alpha.*y.*(2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) - lambda.*y.*(4.*pi.^2.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) - 4.*pi.*y.*cos(2.*pi.*x.*y) - 4.*pi.*z.*sin(2.*pi.*x) + 1) - 8.*pi.*mu.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - mu.*y + 4.*pi.*mu.*(-pi.*x.^2.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*(x - 1).*cos(2.*pi.*x.*y) + y.*z.*sin(2.*pi.*x)));
        force{2} = @(x, y, z) (-alpha.*(x - 1).*(2.*pi.*x.*y.*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)) + 2.*pi.*lambda.*y.*(-2.*pi.*x.*y.*(x - 1).*sin(2.*pi.*x.*y) + x.*cos(2.*pi.*x.*y) + (x - 1).*cos(2.*pi.*x.*y)) - lambda.*(x - 2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + 2.*z.*cos(2.*pi.*x) - sin(2.*pi.*x.*y)) - mu.*x - 4.*mu.*z.*cos(2.*pi.*x) + mu.*(-4.*pi.^2.*x.*y.^2.*(x - 1).*sin(2.*pi.*x.*y) + 2.*pi.*x.*y.*cos(2.*pi.*x.*y) + 4.*pi.^2.*y.^2.*z.*cos(2.*pi.*x) + 4.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + sin(2.*pi.*x.*y)));
        force{3} = @(x, y, z) (-2.*y.*(lambda + mu).*cos(2.*pi.*x));

        % Fluid source term
        src = @(x, y, z) (-4.*pi.*K.*tau.*x.*(x - 1).*(pi.*x.*y.*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) - 4.*pi.*K.*tau.*y.^2.*(pi.*y.*(x - 1).*sin(2.*pi.*x.*y) - cos(2.*pi.*x.*y)) + alpha.*y.*(x - 2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + 2.*z.*cos(2.*pi.*x) - sin(2.*pi.*x.*y)) - rho.*y.*(x - 1).*sin(2.*pi.*x.*y));

        % Stress (we do not follow here the convention that the third should be multiplied by two)
        %
        %            | stress{1}, stress{6}, stress{5} |  
        %  stress =  | stress{6}, stress{2}, stress{4} |
        %            | stress{5}, stress{4}, stress{3} |

        % Voigt components
        stress{1} = @(x, y, z) (lambda.*y.*(x - 2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + 2.*z.*cos(2.*pi.*x) - sin(2.*pi.*x.*y)) + 1.0.*mu.*(4.*pi.*y.^2.*(1 - x).*cos(2.*pi.*x.*y) - 2.*y.*sin(2.*pi.*x.*y)));
        stress{2} = @(x, y, z) (lambda.*y.*(x - 2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + 2.*z.*cos(2.*pi.*x) - sin(2.*pi.*x.*y)) + 4.0.*mu.*y.*z.*cos(2.*pi.*x));
        stress{3} = @(x, y, z) (lambda.*y.*(x - 2.*pi.*y.*(x - 1).*cos(2.*pi.*x.*y) + 2.*z.*cos(2.*pi.*x) - sin(2.*pi.*x.*y)) + 2.0.*mu.*x.*y);
        stress{4} = @(x, y, z) (1.0.*mu.*(x.*z + y.^2.*cos(2.*pi.*x)));
        stress{5} = @(x, y, z) (1.0.*mu.*y.*z);
        stress{6} = @(x, y, z) (1.0.*mu.*(2.*pi.*x.*y.*(1 - x).*cos(2.*pi.*x.*y) - 2.*pi.*y.^2.*z.*sin(2.*pi.*x) + (1 - x).*sin(2.*pi.*x.*y)));

      otherwise
        error('d not recognized');
    end

    output.u_fun      = u;
    output.p_fun      = p;
    output.force_fun  = force;
    output.src_fun    = src;
    output.stress_fun = stress;
end

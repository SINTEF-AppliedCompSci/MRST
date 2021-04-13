function [u, f, mu] = analyticalReferencePaper(d, kappa, mu0, lamb0)
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
%
% Reference solution from paper
%  title={Finite volume methods for elasticity with weak symmetry},
%  author={Keilegavlen, Eirik and Nordbotten, Jan Martin},
%  journal={International Journal for Numerical Methods in Engineering},
%  volume={112},
%  number={8},
%  pages={939--962},
%  year={2017},
%  publisher={Wiley Online Library}
%
%  The displacement u which is returned by this function has not been
%  multiplied by the indicator function, as it appears in the paper. 


    % displacement
    u = cell(d, 1);
    % force
    f = cell(d, 1);

    switch d 
      case 2
        mu = @(x, y) (1 - Xindfunction2D(x, y)) + kappa*Xindfunction2D(x, y);
        u{1} = @(x, y) 1./mu(x, y).*((x - 0.5).^2.*(y - 0.5).^2);
        u{2} = @(x, y) 1./mu(x, y).*(-2/3*(x - 0.5).*(y - 0.5).^3);

        f{1} = @(x, y) mu0*((2*x - 1).^2 + (2*y - 1).^2)/2;
        f{2} = @(x, y) -mu0*(2*x - 1).*(2*y - 1);

      case 3
        mu = @(x, y, z) (1 - Xindfunction3D(x, y, z)) + kappa*Xindfunction3D(x, ...
                                                          y, z);

        u{1} = @(x,y,z) 1./mu(x, y, z).*(((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2);
        u{2} = @(x,y,z) 1./mu(x, y, z).*(((x - 0.5).^2).*(y - 0.5).^2.*(z - 0.5).^2);
        u{3} = @(x,y,z) 1./mu(x, y, z).*(-2/3*((x - 0.5).*(y - 0.5).^2 + (x - 0.5).^2.*(y - 0.5)).*(z - 0.5).^3);

        f{1} = @(x, y, z) (mu0*((2*x - 1).*(2*z - 1).^2.*(2*x + 4*y - 3) + 2*(2*y - 1).^2.*(2*z - 1).^2 + (2*y - 1).*((2*x - 1).^2.*(2*y - 1) + (2*z - 1).^2.*(-4*x - 2*y + 3)))/8);
        f{2} = @(x, y, z) (mu0*(2*(2*x - 1).^2.*(2*z - 1).^2 + (2*x - 1).*((2*x - 1).*(2*y - 1).^2 + (2*z - 1).^2.*(-2*x - 4*y + 3)) + (2*y - 1).*(2*z - 1).^2.*(4*x + 2*y - 3))/8);
        f{3} = @(x, y, z) (mu0*(2*z - 1).*((1 - 2*x).*(2*z - 1).^2 + (1 - 2*y).*(2*z - 1).^2 + (2*x - 1).^2.*(6*y - 3) - 12*(2*x - 1).*(2*y - 1).*(x + y - 1) + (6*x - 3).*(2*y - 1).^2)/12);

      otherwise
        error('d not recognized')
    end
end

function u = Xindfunction2D(x, y)
    ix = (x > 0.5) & (y > 0.5);
    u = 0*x;
    u(ix) = 1;
end

function u = Xindfunction3D(x, y, z)
    ix = (x > 0.5) & (y > 0.5) & (z > 0.5);
    u = 0*x;
    u(ix) = 1;
end

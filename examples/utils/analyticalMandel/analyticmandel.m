function p = analyticmandel(t, x, params, varargin)
% Compute analytical solution for pressure in normalized coordinates

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

    
    opt = struct('num_modes', 100);
    opt = merge_options(opt, varargin{:});
    
    N = opt.num_modes;
    
    nu = params.nu;
    cW = params.cW;
    
    % only incompressible case for the moment
    assert(cW == 0, 'only consider incompressible case for the moment');
    eta = (1 - nu)/(1 - 2*nu);
    
    % compute values of alphas
    alphas = computealphas(eta, N);

    docheck = false;
    if docheck
        % check alphas
        figure
        hold on
        a = (0 : 1e-2: 10)';
        plot(a/pi, tan(a)/pi);
        plot(a/pi, 2*eta*a/pi);
        axis([0, 10/pi, 0, 2*eta*10/pi]);
    end
    
    p = 0*x;
    casename = 'coussy';
    switch casename
      case 'coussy'
        for i = 1 : N
            % formula from Coussy
            an = alphas(i);
            cn = cos(an*x) - cos(an);
            cn = (sin(an)/(an - sin(an)*cos(an)))*cn;
            cn = exp(-an^2*t)*cn;
            p = p + 2*cn;
        end
      case 'verruijt'
        for i = 1 : N
            % formula from Verruijt (We have checked that they are equivalent to Coussy - up a constant coefficient)
            an = alphas(i);
            cn = cos(an*x) - cos(an);
            cn = cos(an)/(1 - 2*eta*(cos(an).^2))*cn;
            cn = cn*exp(-an^2*t);
            p = p + 2*eta*cn;
        end
    end

end

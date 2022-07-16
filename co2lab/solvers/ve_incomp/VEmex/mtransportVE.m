function varargout = mtransportVE(sol, Gt, dT, rock, fluid, varargin)  
% Wrapper for the mex function for explicit transport for VE
%
% SYNOPSIS:
% mode 1: delete VE transport solver
%  mtransportVE();  
% mode 2: run VE transport solver
% [state.h, state.h_max] = mtransportVE(state, G_top, tf, rock, ...
%                                          fluid, 'pn', pv1);  
%
% DESCRIPTION:
%   Function mtransportVE solves the Buckley-Leverett transport
%   equation
%
%        h_t + f(h)_x = q
%
%   using a first-order upwind discretisation in space and a forward Euler
%   discretisation in time.  The transport equation is solved on the time
%   interval [0,tf].
%
%   The upwind forward Euler discretisation of the Buckley-Leverett model
%   for the Vertical Equilibrium model can be written as:
%
%     h^(n+1) = h^n - (dt./pv)*((H(h^n) - max(q,0) - min(q,0)*f(h^n))
%
%   where
%        H(h) = f_up(h)(flux + grav*lam_nw_up*(z_diff+rho_diff*h_diff(h)))
%   
%   z_diff, h_diff are two point approximations to grad_x z, and grad_x h,
%   f_up and lam_nw_up are the Buckely-Leverett fractional flow function
%   and the mobility for the non-wetting phase, respectively, evaluated for
%   upstream mobility:
%
%        f_up = *A_w*lam_w(h)./(A_w*lam_w(h)+A_nw*lam_nw(h))
%        lam_nw_up = diag(A_nw*lam_nw(h)
%
%   pv is the porevolume, lam_x is the mobility for phase x, while A_nw
%   and A_w are index matrices that determine the upstream mobility.
%
%
% ASSUMPTIONS
%  * injection is CO2 only
%  * cell relperm: - computed using vertical integration
%                  - does not account for residual CO2 in water phase 
%  * poro: computed as average of column in z-direction
%  * dt is estimated using a method proposed by Coats
%
% PARAMETERS:
%   state   - Reservoir solution structure containing valid water
%             saturation state.h(:,1) with one value for each cell
%             in the grid.
%
%   G_top   - Grid data structure discretising the top surface of the 
%             reservoir model, as defined by function 'topSurfaceGrid'.
%
%   tf      - End point of time integration interval (i.e., final time),
%             measured in units of seconds.
%
%   fluid   - Data structure as defined by function 'initVEFluid'.
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%
%   wells, src, bc - Source terms
%
%   verbose       - Whether or not time integration progress should be
%                   reported to the screen.
%                   Default value: verbose = false.
%
%   gravity       - Gravity acceleration strength. Default value = 0.0
%
%
% RETURNS:
%   h     - Thickness of CO2.
%   
%   h_max - Maximal thickness of CO2.   

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

if nargin == 0
   VETransportCPU();
   varargout = {};
else
   [h, max_h] = VETransportCPU(sol, Gt, dT, rock, fluid, varargin{:});  
   varargout{1} = h;
   varargout{2} = max_h;   
end

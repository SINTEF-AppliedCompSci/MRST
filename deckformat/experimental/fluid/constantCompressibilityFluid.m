function fluid = initCompressibleFluid(rho0, p0, c0, mu0, varargin)
%Construct PVT and relperm functions, constant compr., single phase.
%
% SYNOPSIS:
%   fluid = initCompressibleFluid(rho0, p0, c0, mu0)
%
% PARAMETERS:
%   rho0 - Fluid reference density, in units of kg/m^3.
%
%   p0   - Pressure, in units of Pascal, at which fluid reference
%          parameters (density, compressibility, and viscosity) are
%          defined.
%
%   c0   - Fluid reference compressibility at pressure 'p0'.
%
%   mu0  - Fluid reference viscosity, in units of Pa s, at pressure 'p0'.
%
% RETURNS:
%   fluid - A 1-by-1 structure array of PVT and relative permeability
%           evaluation functions.  The individual evaluation functions
%           (structure fields) are,
%             - pvt -- Evaluates pvt data for the single phase compressible
%                      fluid.  Specifically, given an n-by-1 vector of
%                      pressures p and an n-by-1 vector of surface volumes
%                      z, the call
%
%                         [c, rho, mu, u, B, R] = fluid.pvt(p, z)
%
%                      computes, respectively, n-by-1 phase
%                      compressibilities (c), densities (rho), viscosities
%                      (mu), and reservoir volumes (u).
%
%              - relperm --
%                      Evaluates one-phase relative permeability curves.
%                      Specifically, given an n-by-1 vector of fluid
%                      saturations, s (implicitly assumed s = ones([n,1])),
%                      the call
%
%                         kr = fluid.relperm(s)
%
%                      computes an n-by-1 array of relative permeability
%                      values.  As 'fluid' is a single phase fluid model,
%                      all relative permeability values are one, i.e.,
%
%                         kr = ones(size(s))
%
%
%              - surfaceDensity--
%
%              - info --
%                      Human-readable summary for fluid.
%
% EXAMPLE:
%
% SEE ALSO:
%   `initBlackoilFluid`.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


   function [c, rho, mu, u, B, R] = pvtfun(varargin)
      % Poor man's class replacement
      assert(nargin==2, 'fluid.pvt takes p and z as arguments.')
      p = varargin{1};
      z = varargin{2};

      assert(numel(p) == size(z, 1), ...
         'In fluid.pvt, p and z must have same number of rows');
      assert(size(z, 2) ==1, ...
         'In fluid.pvt,  z must have one column.' )

      c    = repmat(c0,  [numel(p), 1]);
      rho  = rho0 .* exp(c .* (p - p0));
      mu   = repmat(mu0, [numel(p), 1]);
      u    = rho0 .* z ./ rho;
      B    = u./z;
      R    = zeros([numel(p), 1]);
   end

   fluid.pvt            = @pvtfun;
   fluid.relperm        = @(s, varargin) ones([size(s,1), 1]);
   fluid.info           = 'single compressible fluid';
   fluid.surfaceDensity = rho0;
end

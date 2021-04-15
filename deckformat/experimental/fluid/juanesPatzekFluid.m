function fluid = initSimpleThreephaseCompressibleFluid(rho0, p0, c0, mu0, varargin)
%Construct PVT and relperm functions, constant compressibility three-phase
%mode a la Juanes&Patzek
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
%                         [c, rho, mu, u] = fluid.pvt(p, z)
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


   function [c, rho, mu, u, R] = pvtfun(p, z)
      if numel(p) ~= size(z, 1),
         error('There must be one pressure for each row in mass table');
      end

      np   = numel(rho0);
      c    = repmat(c0,  [numel(p), 1]);
      rho  = bsxfun(@times, rho0, exp(bsxfun(@times, c0, (p - p0))));
      mu   = repmat(mu0, [numel(p), 1]);
      u    = bsxfun(@times, rho0, z ./ rho);
      R    = zeros([numel(p), np]);
   end

   function [kr, dkr] = relpermfun(s)
      % Juanes&Patzek "Three-Phase Displacement Theory: An Improved
      % Description of Relative Permeabilities", eqs. 32--34.
      krw  = @(sw)     sw.*sw;
      dkrw = @(sw)    [2.*sw, -2.*sw, 0.*sw];
      krg  = @(sg)     sg.*sg;
      dkrg = @(sg)    [0.*sg, -2.*sg, 2.*sg];
      kro  = @(sw,sg)  (1-sw-sg).*(1-sw).*(1-sg);
      dkro = @(sw,sg) [-(2-2*sw-sg).*(1-sg), ...
                       -4+5.*sw+5.*sg-4.*sw.*sg-sw.^2-sg.^2, ...
                       -(2-2.*sg-sw).*(1-sw)];
     %{
      krw  = @(sw)     sw;
      dkrw = @(sw)    repmat([1, -1, 0], [numel(sw), 1]);
      krg  = @(sg)     sg;
      dkrg = @(sg)    repmat([0, -1, 1], [numel(sg), 1]);
      kro  = @(sw,sg)  (1-sw-sg);
      dkro = @(sw,sg) repmat([-1, 2, -1], [numel(sw), 1]);
      %}

      kr  = [ krw(s(:,1)),  kro(s(:,1), s(:,3)),  krg(s(:,3))];
      dkr = [dkrw(s(:,1)), dkro(s(:,1), s(:,3)), dkrg(s(:,3))];
      dkr = spblockdiag(dkr', [3,3]);
   end



   fluid.pvt     = @pvtfun;
   fluid.relperm = @relpermfun;
   fluid.surfaceDensity = reshape(rho0, 1, []);
   fluid.viscosity = ones([1, numel(rho0)]) * centi * poise();
   fluid. miscible = false([1, numel(rho0)]);
   fluid.info = ['(Almost) Juanes&Patzek "Three-Phase Displacement Theory: An Improved ', ...
                 'Description of Relative Permeabilities", eqs. 32--34.'];
end



% -------------------------------------------------------------------------

function S = spblockdiag(v, blocksize)
   % Construct sparse block diagonal matrix  with (uniform) block size
   % from values v.
   assert (numel(blocksize)==2);

   [i,j] = ndgrid(1:blocksize(1), 1:blocksize(2));

   n     = prod(blocksize);
   N     = numel(v)/n;

   i = repmat(i, [1,N]) + kron(blocksize(1)*(0:N-1), ones(blocksize));
   j = repmat(j, [1,N]) + kron(blocksize(2)*(0:N-1), ones(blocksize));
   S = sparse(i, j, v);
end


function fluid = initEclipseFluid(deck, varargin)
%Construct an MRST fluid object from an ECLIPSE input deck
%
% SYNOPSIS:
%   fluid = initEclipseFluid(deck)
%
% PARAMETERS:
%   deck - An ECLIPSE input deck as defined by function 'readEclipseDeck'.
%
% RETURNS:
%   fluid - A 1-by-1 structure array of PVT and relative permeability
%           evaluation functions.  The individual evaluation functions
%           (structure fields) are,
%             - pvt -- Evaluates pvt data, i.e.,
%
%                         [c, rho, mu, u, R, B] = fluid.pvt(p, z)
%
%                      computes, respectively, n-by-np phase
%                      compressibilities (c), densities (rho), viscosities
%                      (mu), and reservoir volumes (u), dissolution and
%                      vaporization ratios (R) and formation volume
%                      factors (B).
%
%             - relperm --
%                      Relative permeability function,
%
%                          kr       = fluid.relperm(s, state)
%                         [kr, dkr] = fluid.relperm(s, state)
%
%             - surfaceDensity --
%                      Surface densities of primary components (standard
%                      units, kilogram/meter^3).
%
%             - info --
%                      Human-readable summary for fluid.
%
%             - names --
%                      Textual representation of phases.
%
% EXAMPLE:
%   deck  = readEclipseDeck(fn);
%   deck  = convertDeckUnits(deck);
%   fluid = initEclipseFluid(deck);
%
% SEE ALSO:
%   `initCompresibleFluid`.

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


   opt = struct('verbose', mrstVerbose);
   opt = merge_options(opt, varargin{:});

   [pvtfun, rho_s, n, info1, incomp] = eclipsePhaseProperties(deck);
   [relperm, pcap, info2]            = eclipseRelperm        (deck);

   fluid.names = n;
   fluid.info  = [info1, info2];
   if isfield(deck.RUNSPEC, 'TITLE'),
      fluid.info = [sprintf('Title of parameter file: ''%s''\n\n',...
                            deck.RUNSPEC.TITLE), ...
                    fluid.info];
   end

   if size(rho_s, 1) > 1,
      % Multiple PVT regions not supported.  Pick region one.
      rho_s     = rho_s(1, :);
      func      = cellfun(@(f) f{1}, pvtfun, 'UniformOutput', false);
      func      = clear_absent(func, n);
      fluid.pvt = @(p, z) blackOilPVT(func, rho_s, p, z);
      clear func

      dispif(mrstVerbose, ...
            ['Multiple PVT regions not supported.\n', ...
             '\t-> Using region one (1) only.\n\n']);
   else
      pvtfun = clear_absent(pvtfun, n);
      fluid.pvt = @(p, z) blackOilPVT(pvtfun, rho_s, p, z);
   end

   fluid.relperm        = relperm;
   fluid.pc             = pcap;

   %fluid.pcinv          = pcinv;

   fluid.surfaceDensity = rho_s;

   fluid.properties = @(state) properties(state, incomp, fluid.pvt);
   fluid.saturation = @(state) saturation(state, fluid.pvt);

   if opt.verbose,
      fprintf('----------------- Fluid summary ---------------------\n\n');
      if isfield(deck.RUNSPEC, 'TITLE'),
         fprintf('Title of parameter file: ''%s''\n\n', deck.RUNSPEC.TITLE);
      end
      fprintf('PVT properties:\n%s\n', info1);
      fprintf('Relative permeability:\n%s', info2);
      fprintf('-----------------------------------------------------\n\n');
   end
end

%--------------------------------------------------------------------------

function pvtfun = clear_absent(pvtfun, present)
   phases = { 'WATER', 'OIL', 'GAS' };

   [i, j] = blockDiagIndex(numel(phases), numel(present));

   m = strcmpi(reshape(phases (i), [], 1), ...
               reshape(present(j), [], 1));

   pvtfun(~ any(reshape(m, numel(phases), []), 2)) = { [] };
end

%--------------------------------------------------------------------------

function [mu, rho] = properties(x, incomp, pvt)
   if ~all(incomp),
      error('This fluid does not support incompressible ''properties''.');
   end

   [rho, rho, mu] = pvt(1, zeros([1, size(x.s, 2)]));                  %#ok
end

%--------------------------------------------------------------------------

function s = saturation(state, pvt)
   if isfield(state, 'z'),
      [u, u, u, u] = pvt(state.pressure, state.z);                     %#ok
      s            = bsxfun(@rdivide, u, sum(u,2));
   else
      s = state.s;
   end
end

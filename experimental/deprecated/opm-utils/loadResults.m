function res = loadResults(odir, np, varargin)
%Load time-stepping results from simulation
%
% SYNOPSIS:
%   res = loadResults(output_dir, np)
%   res = loadResults(output_dir, np, 'pn1', pv1, ...)
%
% PARAMETERS:
%   output_dir - Master directory containing all simulation results.
%                Individual result vectors are assumed to be contained in
%                named directories directly below 'output_dir'.
%
%                The following result vectors are currently supported
%
%                    concentration, flux, pressure, reorder_it,
%                    residualcounts, saturation, surfvolume
%
%   np         - Number of fluid phases.  Needed to import saturations and
%                surface volumes
%
%   'pn'/pv    - List of 'key'/value pairs defining optional parameters.
%                The supported options are:
%                  - 'steps' -- List of time steps for which to load
%                               simulation results.  Must be numeric and
%                               valid indices between one (1) and the total
%                               number of simulation steps (inclusive).
%                               Default value: steps = [] (Load *all*
%                               simulation results).
%
% RETURNS:
%   res - Simulation results.  Structure containing the result vectors
%         detected during directory scanning of 'output_dir'.  The elements
%         of this structure are cell arrays, the elements of which contain
%         the actual simulation results--possibly restricted to those
%         results pertaining to a specified time step sequence.
%
% EXAMPLE:
%   % Load all simulation results in 'test-output' from a three-phase
%   % simulation of a 60-by-220-by-5 Cartesian grid model, and plot the
%   % fifth pressure field as well as the saturation of the second
%   % component of the twentieth time step:
%
%   G   = cartGrid([60, 220, 5]);
%   res = loadResults('test-output', 3);
%
%   figure
%   plotCellData(G, res.pressure{5}, ...
%                'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.85)
%   view(3)
%
%   figure
%   plotCellData(G, res.saturation{20}(:,2), ...
%                'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceAlpha', 0.85)
%   view(3)

%{
Copyright 2009-2015 SINTEF ICT, Applied Mathematics.

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

% $Date: 2012-09-14 20:40:02 +0200 (Fri, 14 Sep 2012) $
% $Revision: 9694 $

   opt = struct('steps', []);
   opt = merge_options(opt, varargin{:});

   steps = opt.steps;
   if ~isempty(steps),
      assert (isnumeric(steps), '''STEPS'' must be numeric step numbers.');

      steps = reshape(steps, 1, []);
   end

   id     = @(x) x;
   pvec   = @(x) reshape(x(:), np, []) .';
   rcount = @(x) accumarray([x(:,2) + 1, 2 - x(:,1)], 1);

   res = struct();
   elm = dir(odir);
   nm  = { elm.name };

   for e = nm([ elm.isdir ] & ~ strncmp(nm, '.', 1)),
      assert (~isfield(res, e{1}), ...
              'Vector ''%s'' encountered multiple times?', e{1});

      switch e{1},
         case { 'concentration', 'flux', 'pressure', 'reorder_it' },
            t = id;

         case { 'saturation', 'surfvolume' },
            t = pvec;

         case 'residualcounts',
            t = rcount;

         otherwise
            dispif(mrstVerbose, ...
                   'Unknown result vector ''%s''. Ignored.\n', e{1});

            continue
      end

      res.(e{1}) = load_vector(fullfile(odir, e{1}), t, steps);
   end
end

%--------------------------------------------------------------------------

function v = load_vector(d, transform, steps)
   f = dir(d);
   i = sortrows([vertcat(f.datenum), (1 : numel(f)) .']);
   p = ~ vertcat(f(i(:,2)).isdir);  % Non-directories

   fn = { f(i(p,2)).name };
   
   if ~isempty(steps),
      assert (all(1 <= steps) && all(steps <= numel(fn)), ...
              '''STEPS'' must be valid step index numbers.');

      fn = fn(steps);
   end

   v = cell(size(fn));
   k = 1;

   for s = fn,
      v{k} = transform(load(fullfile(d, s{1})));
      k    = k + 1;
   end
end

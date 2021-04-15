function ExTable = expandMiscibleTable(T, varargin)
%Expand Miscible PVT Table Definition to Format Suitable for Visualization
%
% SYNOPSIS
%   ExTable = expandMiscibleTable(T)
%   ExTable = expandMiscibleTable(T, 'pn1', pv1, ...)
%
% PARAMETERS:
%   T      - Miscible table, or cell-array of miscible tables, as defined
%            by function readEclipseDeck and, possibly, function
%            convertDeckUnits.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are
%
%               - inner -- Name (string) of inner table range.
%                          Default: 'inner'.
%
%               - outer -- Name (string) of outer table range.
%                          Default: 'outer'.
%
% RETURNS:
%   ExTable - Structure of expanded miscible property table definitions.
%             One element for each input table (typically: PVT region)
%             identified by the argument 'T'.
%
% EXAMPLE:
%   % Input simulation deck, convert to MRST's strict SI unit conventions.
%   deck = convertDeckUnits(readEclipsDeck(file))
%
%   % Expand PVTO and PVTG Tabulated (Miscible) Property Functions.
%   ExPVTO = expandMiscibleTable(deck.PROPS.PVTO, ...
%                                'inner', 'Po', 'outer', 'Rs')
%   ExPVTG = expandMiscibleTable(deck.PROPS.PVTG, ...
%                                'inner', 'Rv', 'outer', 'Pg')
%
%   %----------------------------------------------------------------
%
%   % Visualize First Region's 'Bo' Factor (i.e., the Oil FVF).
%   [Po, Rs] = deal(convertTo(ExPVTO(1).Po, barsa), ExPVTO(1).Rs);
%
%   figure
%   surf(Po, Rs, ExPVTO(1).B, ...
%        'EdgeColor', 'k', 'FaceColor', 'interp', ...
%        'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
%
%   % Add Contour Plot of 'Bo' for Additional Visual Aid
%   hold on, contour(Po, Rs, ExPVTO(1).B, 101)
%
%   title('Oil Formation Volume Factor (B_o) [rm^3/Sm^3]')
%   xlabel('Oil Pressure (P_o) [Bar]')
%   ylabel('Dissolved Gas-Oil Ratio (R_s) [Sm^3/Sm^3]')
%
%   %----------------------------------------------------------------
%
%   % Visualize Last Region's Tabulated Gas Viscosity Function
%   [Rv, Pg, mu_g] = deal(ExPVTG(end).Rv, ...
%                         convertTo(ExPVTG(end).Pg, barsa), ...
%                         convertTo(ExPVTG(end).mu, centi*poise));
%
%   figure
%   surf(Pg, Rv, mu_g, ...
%        'EdgeColor', 'k', 'FaceColor', 'interp', ...
%        'Marker', '*', 'MarkerEdgeColor', 'r', 'MarkerSize', 5)
%
%   % Add Contour Plot of 'mu_g' for Additional Visual Aid
%   hold on, contour(Pg, Rv, mu_g, 101)
%
%   title('Gas Viscosity (\mu_g) [cP]')
%   xlabel('Gas Pressure (P_g) [Bar]')
%   ylabel('Vaporized Oil-Gas Ratio (R_v) [Sm^3/Sm^3]')
%
% SEE ALSO:
%   readEclipseDeck, convertDeckUnits.

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

   if isstruct(T)
      T = { T };
   end

   if ~iscell(T)
      error('Argument:Invalid', ...
            'Input Must Be Structure or Cell-Array of Structures');
   end

   check_input(T);

   [outer, inner, B, mu] = cellfun(@expand, T, 'UniformOutput', false);

   [oname, iname] = get_axis_names(varargin{:});

   ExTable = reshape(struct(oname, outer, ...
                            iname, inner, ...
                            'B', B, 'mu', mu), size(T));
end

%--------------------------------------------------------------------------

function check_input(T)
   is_valid = ...
      @(t) isstruct(t) && (numel(t) == 1) && ...
           all(isfield(t, { 'key', 'pos', 'data' }));

   assert (all(cellfun(is_valid, T)), ...
           'Inputs Must Be Scalar, Miscible Table Structures');
end

%--------------------------------------------------------------------------

function [outer, inner, B, mu] = expand(T)
   n = reshape(diff(T.pos), [], 1);
   m = max(n);

   % NaN => Undefined values don't show up in plots.
   [outer, inner, B, mu] = deal(NaN([m, numel(n)]));

   start = reshape((0 : (numel(n) - 1)) .* m, [], 1);
   ix    = mcolon(start + 1, start + reshape(n, size(start)));

   outer(ix) = rldecode(reshape(T.key, [], 1), n);
   inner(ix) = T.data(:,1);
   B    (ix) = T.data(:,2);
   mu   (ix) = T.data(:,3);
end

%--------------------------------------------------------------------------

function [oname, iname] = get_axis_names(varargin)
   opt = struct('inner', 'inner', 'outer', 'outer');
   opt = merge_options(opt, varargin{:});

   oname = opt.outer;
   iname = opt.inner;
end

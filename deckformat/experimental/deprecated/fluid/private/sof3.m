function [krow, krog, Soco, Socr, Sorw, Sorg, Somax, krsmax] = ...
      sof3(T, varargin)
%Construct oil-water and oil-gas relperm eval. functions from SOF3 table.
%
% SYNOPSIS:
%   [krow, krog, Soco, Socr, Sorw, Sorg, Somax, krSomax] = sof3(table)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SOF3'.  May be a numeric array or a a cell array of
%           such arrays.
%
% RETURNS:
%   In the following, 'ntab' refers to the number of input tables
%   represented by the input parameter 'table'.  If 'table' is a single
%   numeric array, then ntab = 1, otherwise (i.e., when 'table' is a cell
%   array of numeric arrays) ntab = numel(table).
%
%   In an ECLIPSE input deck, 'ntab' typically corresponds to the 'NTSFUN'
%   parameter specified in item one of the 'TABDIMS' keyword.
%
%   krow    - A 1-by-ntab cell array of function handles such that the call
%             'krow(so)' computes the relative permeability of oil at the
%             *oil* saturation 'so'.  When called with two output
%             parameters, i.e., [ko{1:2}] = krow(so), then the second
%             output parameter is the derivative of the oil relative
%             permeability with respect to the *oil* saturation at 'so'.
%
%   krog    - A 1-by-ntab cell array of function handles such that the call
%             'krog(so)' computes the relative permeability of oil at the
%             *oil* saturation 'so'.  When called with two output
%             parameters, i.e., [ko{1:2}] = krog(so), then the second
%             output parameter is the derivative of the oil relative
%             permeability with respect to the *oil* saturation at 'so'.
%
%   Soco    - A 1-by-ntab numeric array of connate oil saturations.
%
%   Socr    - A 1-by-ntab numeric array of critical oil saturations.  This
%             is the largest oil saturation for which *both* krow(So)=0 and
%             krog(So)=0.
%
%   Sorw    - A 1-by-ntab numeric array of Residual oil saturation in the
%             two-phase oil-water system.  This is the saturation at which
%             the oil relative permeability krow(so) becomes zero.
%
%   Sorg    - A 1-by-ntab numeric array of Residual oil saturation in the
%             two-phase oil-gas system.  This is the saturation at which
%             the oil relative permeability krog(so) becomes zero.
%
%   Somax   - A 1-by-ntab numeric array of maximum water saturations.
%
%   krSomax - A 1-by-ntab numeric array of oil relative permeabilities at
%             corresponding maximum oil saturations, Somax.
%
% SEE ALSO:
%   `readRelPermTable`, `swfn`, `sgfn`, `sof2`, `swof`, `sgof`.

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


   assert (isnumeric(T) || (iscell(T) && all(cellfun(@isnumeric, T))),  ...
          ['Parameter ''table'' must be a numeric array or cell array', ...
           ' of numeric arrays.']);

   if isnumeric(T), T = { T }; end

   ntab = numel(T);

   krow  = cell ([1, ntab]);
   krog  = cell ([1, ntab]);
   Soco  = zeros([1, ntab]);
   Socr  = zeros([1, ntab]);
   Somax = zeros([1, ntab]);
   Sorw  = zeros([1, ntab]);
   Sorg  = zeros([1, ntab]);
   krsmax = zeros([1, ntab]);

   for k = 1 : ntab,
      [Soco(k), Socr(k), Sorw(k), Sorg(k), Somax(k), krsmax(k)] = ...
         table_points(T{k});

      if Soco(k) > 0, T{k} = [0, T{k}(1, 2:end); T{k}]; end

      krow{k} = @(so) interp_water(so, T{k}, Soco(k));
      krog{k} = @(so) interp_water(so, T{k}, Soco(k));
   end
end

%--------------------------------------------------------------------------

function [Soco, Sorw, Sorg, Socr, Somax, krsmax] = table_points(T)
   % Basic validation of input data.
   assert (~ any(any(T < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~ (T(1, 2) > 0), 'krow(Soco) must be zero (SOF3).');
   assert (~ (T(1, 3) > 0), 'krog(Soco) must be zero (SOF3).');

   assert (all (diff(T(:,1)) > 0), ...
          ['Oil saturation must be strictly increasing down ', ...
           'column (SOF3).']);
   assert (~any(diff(T(:,2)) < 0), ...
          ['Oil relative permeability in oil/water system must ', ...
           'level or increase down column (SOF3).']);
   assert (~any(diff(T(:,3)) < 0), ...
          ['Oil relative permeability in oil/gas system must ', ...
           'level or increase down column (SOF3).']);

   assert (T(end,2) == T(end,3), ...
          ['Oil relative permeability in oil/water and oil/gas ', ...
           'systems must be equal at maxmimum oil (connate water) ', ...
           'saturation (SOF3).']);

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate oil saturation (>= 0) is first So encountered in table.
   %
   Soco = T(1,1);  assert (~ (Soco < 0));

   % 2) Residual oil saturation in oil-water system.
   i    = find(T(:,2) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Sorw = T(i - 1, 1);                   assert (~ (Sorw < Soco));

   % 3) Residual oil saturation in oil-gas system.
   i    = find(T(:,3) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Sorg = T(i - 1, 1);                   assert (~ (Sorg < Soco));

   % 4) Critical water saturation (>= 0) is largest So for which both
   %    krow(So)=0 and krog(So)=0.  This is the minimum of 'Sorw' and
   %    'Sorg'.
   Socr = min(Sorw, Sorg);

   % 5) Maximum oil saturation in table.
   Somax = T(end,1);

   % 6) Oil rel-perm at Somax.
   krsmax = T(end,2);
end

%--------------------------------------------------------------------------

% Water relative permeability as a function of water saturation, sw.
% Note: Water saturation is always >= Swco.
function varargout = interp_water(so, T, Soco)
   so = max(so, Soco);

   varargout{1}    = interpTable(T(:,1), T(:,2), so);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), so);
   end
end

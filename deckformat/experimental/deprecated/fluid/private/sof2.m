function [kro, Soco, Socr, Somax, krsmax] = sof2(T, varargin)
% Construct oil-water or oil-gas relperm eval. functions from SOF2 table.
%
% SYNOPSIS:
%   [kro, Soco, Socr, Somax, krSomax] = sof2(table)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SOF2'.  May be a numeric array or a a cell array of
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
%   kro     - A 1-by-ntab cell array of function handles such that the call
%             'kro(so)' computes the relative permeability of oil at the
%             *oil* saturation 'so'.  When called with two output
%             parameters, i.e., [ko{1:2}] = kro(so), then the second output
%             parameter is the derivative of the oil relative permeability
%             with respect to the *oil* saturation at 'so'.
%
%   Soco    - A 1-by-ntab numeric array of connate oil saturations.
%
%   Socr    - A 1-by-ntab numeric array of critical oil saturations.  This
%             is the largest oil saturation for which kro(So)=0.
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

   kro   = cell ([1, ntab]);
   Soco  = zeros([1, ntab]);
   Socr  = zeros([1, ntab]);
   Somax = zeros([1, ntab]);
   krsmax = zeros([1, ntab]);

   for k = 1 : ntab,
      [Soco(k), Socr(k), Somax(k), krsmax(k)] = table_points(T{k});

      if Soco(k) > 0, T{k} = [0, T{k}(1, 2:end); T{k}]; end

      kro{k} = @(so) interp_oil(so, T{k}, Soco(k));
   end
end

%--------------------------------------------------------------------------

function [Soco, Socr, Somax, krsmax] = table_points(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:2) < 0)), ...
      'Saturation and rel-perm values must be non-negative.');

   assert (~ (T(1, 2) > 0), 'kro(Soco) must be zero (SOF2).');

   assert (all (diff(T(:,1)) > 0), ...
          ['Oil saturation must be strictly increasing down ', ...
           'column (SOF2).']);
   assert (~any(diff(T(:,2)) < 0), ...
          ['Oil relative permeability in oil/water system must ', ...
           'level or increase down column (SOF2).']);

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate oil saturation (>= 0) is first So encountered in table.
   %
   Soco = T(1,1);  assert (~ (Soco < 0));

   % 2) Residual oil saturation
   i    = find(T(:,2) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Socr = T(i - 1, 1);                   assert (~ (Socr < Soco));

   % 3) Maximum oil saturation in table.
   Somax = T(end,1);

   % 4) Oil rel-perm at Somax
   krsmax = T(end,2);
end

%--------------------------------------------------------------------------

% Oil relative permeability as a function of oil saturation, so.
% Note: Oil saturation is always >= Soco.
function varargout = interp_oil(so, T, Soco)
   so = max(so, Soco);

   varargout{1}    = interpTable(T(:,1), T(:,2), so);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), so);
   end
end

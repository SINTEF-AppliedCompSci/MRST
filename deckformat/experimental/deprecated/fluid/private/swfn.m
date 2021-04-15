function [krw, Swco, pc, Swcr, Swmax, pcinv] = swfn(T, varargin)
%Construct water relperm evaluation functions from SWFN table.
%
% SYNOPSIS:
%   [krw, Swco, pc, Swcr, Swmax, pcinv] = swfn(table)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SWFN'.  May be a numeric array or a a cell array of
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
%   krw   - A 1-by-ntab cell array of function handles such that the call
%           'krw(sw)' computes the relative permeability of water at the
%           *water* saturation 'sw'.  When called with two output
%           parameters, i.e., [kw{1:2}] = krw(sw), then the second output
%           parameter is the derivative of the water relative permeability
%           with respect to the *water* saturation at 'sw'.
%
%   Swco  - A 1-by-ntab numeric array of connate water saturations.  Oil
%           saturation in run must not exceed 1-Swco.
%
%   pc    - A 1-by-ntab cell array of function handles such that the call
%           'pc(sw)' computes the capillary pressure. When called with two
%           output parameters, i.e., [pc{1:2}] = pc(sw), then the second
%           output parameter is the derivative of capillary pressure.
%
%   Swcr  - A 1-by-ntab numeric array of critical water saturations.  This
%           is the largest water saturation for which krw(Sw)=0.
%
%   Swmax - A 1-by-ntab numeric array of maximum water saturations.
%
%   pcinv - A 1-by-ntab cell array of function handles such that the call
%           'pcinv(dp)' computes the inverse of the capillary pressure.
%
% SEE ALSO:
%   `readRelPermTable`, `sgfn`, `swof`, `sgof`.

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

   krw   = cell ([1, ntab]);
   pc    = cell ([1, ntab]);
   pcinv = cell ([1, ntab]);
   Swco  = zeros([1, ntab]);
   Swcr  = zeros([1, ntab]);
   Swmax = zeros([1, ntab]);

   for k = 1 : ntab,
      [Swco(k), Swcr(k), Swmax(k)] = table_points(T{k});

      if Swco(k) > 0, T{k} = [0, T{k}(1, 2:end); T{k}]; end

      krw  {k} = @(sw) interp_water(sw, T{k}, Swco(k));
      pc   {k} = @(sw) interp_pc   (sw, T{k}, Swco(k));
      pcinv{k} = @(dp) interp_pcinv(dp, T{k});
   end
end

%--------------------------------------------------------------------------

function [Swco, Swcr, Swmax] = table_points(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:2) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~ (T(1, 2) > 0), 'krw(Swco) must be zero (SWFN).');
   assert (all (diff(T(:,1)) > 0), ...
          ['Water saturation must be strictly increasing down ', ...
           'column (SWFN).']);
   assert (~any(diff(T(:,2)) < 0), ...
          ['Water relative permeability must level or increase ', ...
           'down column (SWFN).']);

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate water saturation (>= 0) is first Sw encountered in table.
   %
   Swco = T(1,1);  assert (~ (Swco < 0));

   % 2) Critical water saturation (>= 0) is highest Sw for which krw(Sw)=0.
   %
   i    = find(T(:,2) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Swcr = T(i - 1, 1);                   assert (~ (Swcr < Swco));

   % 3) Maximal water saturation in table
   Swmax = max(T(:,1));
end

%--------------------------------------------------------------------------

% Water relative permeability as a function of water saturation, sw.
% Note: Water saturation is always >= Swco.
function varargout = interp_water(sw, T, Swco)
   sw = max(sw, Swco);

   varargout{1}    = interpTable(T(:,1), T(:,2), sw);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), sw);
   end
end

%--------------------------------------------------------------------------

% Capillary pressure and derivative as a function of water saturation
function varargout = interp_pc(sw, T, Swco)
   sw = max(sw, Swco);

   varargout{1}    = interpTable(T(:,1), T(:,3), sw);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,3), sw);
   end
end

%--------------------------------------------------------------------------

% Inverse of capillary pressure function
function varargout = interp_pcinv(dp, T)
   [i, i] = unique(T(:,3));                                            %#ok
   varargout{1} = interpTable(T(i,3), T(i,1), dp);
end

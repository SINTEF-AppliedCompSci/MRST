function [krw, kro, krocw, Swco, Swcr, Sowcr, Swmax, pc, pcinv] = swof(T)
%Construct water/oil rel-perm evaluation functions from SWOF table.
%
% SYNOPSIS:
%   [krw, kro, krocw, Swco, Swcr, Sowcr, Swmax, pc, pcinv] = swof(table)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SWOF'.  May be a numeric array or a cell array of such
%           arrays.
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
%   kro   - A 1-by-ntab cell array of function handles such that the call
%           'kro(so)' computes the relative permeability of oil at the
%           *oil* saturation 'so'.  When called with two output parameters,
%           i.e., [ko{1:2}] = kro(so), then the second output parameter is
%           the derivative of the oil relative permeability with respect to
%           the *oil* saturation at 'so'.
%
%   krocw - A 1-by-ntab numeric array of oil relative permeability at
%           So = 1 - Swco (i.e., at connate water).
%
%   Swco  - A 1-by-ntab numeric array of connate water saturations.  Oil
%           saturation in run must not exceed 1-Swco.
%
%   Swcr  - A 1-by-ntab numeric array of critical water saturations.  This
%           is the largest water saturation for which krw(Sw)=0.
%
%   Sowcr - A 1-by-ntab numeric array of critical oil saturations.  This is
%           the largest oil saturation for which kro(So)=0.
%
%   Swmax - A 1-by-ntab numeric array of maximum water saturations.
%
%   pc    - A 1-by-ntab cell array of function handles such that the call
%           'pc(sw)' computes the capillary pressure.  When called with two
%           output parameters, i.e., [pcow{1:2}] = pc(sw), then the second
%           output parameter is the derivative of capillary pressure with
%           respect to the *water* saturation at 'sw'.
%
%   pcinv - A 1-by-ntab cell array of function handles such that the call
%           'pcinv(sw)' computes the inverse of the capillary pressure.
%
% SEE ALSO:
%   `readRelPermTable`, `sgof`.

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
   kro   = cell ([1, ntab]);
   pc    = cell ([1, ntab]);
   pcinv = cell ([1, ntab]);
   Swco  = zeros([1, ntab]);
   Swcr  = zeros([1, ntab]);
   Sowcr = zeros([1, ntab]);
   Swmax = zeros([1, ntab]);
   krocw = zeros([1, ntab]);

   for k = 1 : ntab,
      [Swco(k), Swcr(k), Sowcr(k), Swmax(k), krocw(k)] = table_points(T{k});

      if Swco(k) > 0,
         T{k} = [0, T{k}(1, 2:end); ...
                 T{k}            ];
      end

      if Swmax(k) < 1,
         T{k} = [T{k}                      ; ...
                 1, 1, 0, T{k}(end, 4:end)]; %may not the eclipse truth
      end

      krw  {k} = @(sw) interp_water(sw, T{k}, Swco(k));
      kro  {k} = @(so) interp_oil  (so, T{k}, Swco(k));
      pc   {k} = @(sw) interp_pc   (sw, T{k}, Swco(k));
      pcinv{k} = @(dp) interp_pcinv(dp, T{k}         );
   end
end

%--------------------------------------------------------------------------

function [Swco, Swcr, Sowcr, Swmax, krocw] = table_points(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:3) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~ (T( 1 ,2) > 0), ...
           'Water must be immobile at connate water saturation.');
   assert (~ (T(end,3) > 0), ...
           'Oil must be immobile at maximum water saturation.');

   assert (all (diff(T(:,1)) > 0), ...
           'Water saturation must be monotonically increasing.');
   assert (~any(diff(T(:,2)) < 0), ...
           'Water rel-perm must be level or increasing down column.');
   assert (~any(diff(T(:,3)) > 0), ...
           'Oil rel-perm must be level or decreasing down column.');

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate water saturation (>= 0) is first Sw encountered in table.
   %
   Swco = T(1,1);

   % 2) Critical water saturation (>= 0) is highest Sw for which krw(Sw)=0.
   %
   i = find(T(:,2) > 0, 1, 'first');
   assert (~isempty(i), ...
           'Table must define a critical water saturation, Swcr.');

   Swcr = T(i - 1, 1);

   % 3) Critical oil saturation (>= 0) is highest So for which kro(So)=0.
   %
   %      Sowcr = max { S_o | k_{row}(S_o) = 0 },  S_o \in [0,1]
   %
   %    This is 1-lowest(Sw) at which the same criterion is satisfied.
   %
   i     = find(~(T(:,3) > 0), 1, 'first');
   Sowcr = 1 - T(i, 1);

   % 4) Maximal water saturation in table.
   Swmax = T(end, 1);

   % 5) Oil relperm at Sw = Swco (=> So = 1 - Swco)
   %
   krocw = T(1, 3);
end

%--------------------------------------------------------------------------

% Water relative permeability as a function of water saturation, sw.
% Note: Water saturation is always >= Swco.
function varargout = interp_water(sw, T, Swco)
   sw = max(sw, Swco);

   varargout{1}    =  interpTable(T(:,1), T(:,2), sw);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), sw);
   end
end

%--------------------------------------------------------------------------

% Oil relative permeability as a function of oil saturation, so.
% Note: Oil saturation is always <= 1 - Swco.
function varargout = interp_oil(so, T, Swco)
   sw = max(1 - so, Swco);

   varargout{1}    =   interpTable(T(:,1), T(:,3), sw);

   if nargout > 1,
      varargout{2} = -dinterpTable(T(:,1), T(:,3), sw);
   end
end

%--------------------------------------------------------------------------

% Capillary pressure and derivative as a function of water saturation
function varargout = interp_pc(sw, T, Swco)
   sw = max(sw, Swco);

   varargout{1}    =  interpTable(T(:,1), T(:,4), sw);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,4), sw);
   end
end

%--------------------------------------------------------------------------

% Inverse of capillary pressure function
function varargout = interp_pcinv(dp, T)
   [i, i, j]    = unique(T(:,4));                                      %#ok
   varargout{1} = interpTable(T(i,4), T(i,1), dp);
end

function [krg, kro, krocg, Sgco, Sgcr, Sogcr, Sgmax, pc, pcinv] = ...
      sgof(T, varargin)
%Construct gas/oil rel-perm evaluation functions from SGOF table.
%
% SYNOPSIS:
%   [krg, kro, krocg, ...
%    Sgco, Sgcr, Sogcr, Sgmax, pc, pcinv] = sgof(table)
%   [krg, kro, krocg, ...
%    Sgco, Sgcr, Sogcr, Sgmax, pc, pcinv] = sgof(table, Swco)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SGOF'.  May be a numeric array or a cell array of such
%           arrays.
%
%   Swco  - Connate water saturation.  One scalar, numeric value for each
%           input table represented by 'table'.
%           OPTIONAL.  Zero connate water saturation will be used if this
%           parameter is unspecified.
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
%   krg   - A 1-by-ntab cell array of function handles such that the call
%           'krg(sg)' computes the relative permeability of gas at the
%           *gas* saturation 'sg'.  When called with two output
%           parameters, i.e., [kg{1:2}] = krg(sg), then the second output
%           parameter is the derivative of the gas relative permeability
%           with respect to the *gas* saturation at 'sg'.
%
%   kro   - A 1-by-ntab cell array of function handles such that the call
%           'kro(so)' computes the relative permeability of oil at the
%           *oil* saturation 'so'.  When called with two output
%           parameters, i.e., [ko{1:2}] = kro(so), then the second output
%           parameter is the derivative of the oil relative permeability
%           with respect to the *oil* saturation at 'so'.
%
%   krocg - A 1-by-ntab numeric array of oil relative permeability at
%           So = 1 - Sgco (i.e., at connate gas).
%
%   Sgco  - A 1-by-ntab numeric array of connate gas saturations.  Oil
%           saturation in run must not exceed 1-Sgco.
%
%   Sgcr  - A 1-by-ntab numeric array of critical gas saturations.  This
%           is the largest gas saturation for which krg(sg)=0.
%
%   Sogcr - A 1-by-ntab numeric array of critical oil saturations.  This is
%           the largest oil saturation for which kro(So)=0.
%
%   Sgmax - A 1-by-ntab numeric array of maximum gas saturations.
%
%   pc    - A 1-by-ntab cell array of function handles such that the call
%           'pc(sg)' computes the capillary pressure. When called with two
%           output parameters, i.e., [pcog{1:2}] = pc(sg), then the second
%           output parameter is the derivative of capillary pressure with
%           respect to the *gas* saturation at 'sg'.
%
%   pcinv - A 1-by-ntab cell array of function handles such that the call
%           'pcinv(sg)' computes the inverse of the capillary pressure.
%
% SEE ALSO:
%  `readRelPermTable`, `swof`.

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

   % Set connate water saturation (if water is present).
   Swco = zeros([1, ntab]);
   if nargin > 1 && isnumeric(varargin{1}),
      if numel(varargin{1}) == 1,
         Swco(:) = varargin{1};
      else
         assert (numel(varargin{1}) == ntab, ...
                 'Connate water must be specified in each region');
         Swco = varargin{1};
      end

      assert (~(any(Swco < 0) || any(Swco > 1)), ...
              'Connate water outside permissible region [0,1].');
   end

   krg   = cell ([1, ntab]);
   kro   = cell ([1, ntab]);
   pc    = cell ([1, ntab]);
   pcinv = cell ([1, ntab]);
   krocg = zeros([1, ntab]);
   Sgco  = zeros([1, ntab]);
   Sgcr  = zeros([1, ntab]);
   Sogcr = zeros([1, ntab]);
   Sgmax = zeros([1, ntab]);

   for k = 1 : ntab,
      [krocg(k), Sgco(k), Sgcr(k), Sogcr(k), Sgmax(k)] = ...
         table_points(T{k}, Swco(k));

      krg  {k} = @(sg) interp_gas  (sg, T{k}, Swco(k));
      kro  {k} = @(so) interp_oil  (so, T{k}, Swco(k));
      pc   {k} = @(sg) interp_pc   (sg, T{k}, Swco(k));
      pcinv{k} = @(dp) interp_pcinv(dp, T{k}         );
   end
end

%--------------------------------------------------------------------------

function [krocg, Sgco, Sgcr, Sogcr, Sgmax] = table_points(T, Swco)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:3) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~(T(1,1) > 0), ...
           'Minimum gas saturation in table (Sgco) must be zero.');
   assert (~(T(1,2) > 0), ...
           'Gas must be immobile at minimum gas saturation.');
   assert (~(T(end,3) > 0), ...
           'Oil must be immobile at maximum gas saturation.');

   assert (all(diff(T(:,1)) > 0), ...
           'Gas saturation must be monotonically increasing.');
   assert (~any(diff(T(:,2)) < 0), ...
           'Gas rel-perm must be level or increasing down column.');
   assert (~any(diff(T(:,3)) > 0), ...
           'Oil rel-perm must be level or decreasing down column.');

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate gas saturation
   %
   Sgco = T(1,1);

   % 2) Maximal gas saturation
   %
   Sgmax = T(end,1);

   if     Sgmax > 1 - Swco,

      warning (msgid('Sgmax:NoSwco'), ...
              ['Maximum gas saturation in ''SGOF'' table (=%.2f) ', ...
               'does not account for presence of connate water.'], Sgmax);

   elseif Sgmax < 1 - Swco,

      warning (msgid('Sgmax:LessSwco'), ...
              ['Maximum gas saturation in ''SGOF'' table (=%.2f) is ' , ...
               'less than maximum value derived from connate water.\n', ...
               'The maximum declared gas saturation should normally ' , ...
               'be 1-(connate water saturation) (=%.2f).']            , ...
               Sgmax, 1 - Swco);

   end

   % 3) Critical gas saturation (>= 0) is highest Sg for which krg(Sg) = 0.
   %
   i = find(T(:,2) > 0, 1, 'first');
   assert (~isempty(i), ...
           'Table must define a critical gas saturation, Sgcr.');

   Sgcr = T(i - 1, 1);

   % 4) Oil relative permeability at connate gas.
   krocg = T(1, 3);

   % 5) Critical oil saturation (>= 0) is highest So for which kro(So)=0.
   %
   %      Sogcr = max { S_o | k_{rog}(S_o) = 0 },  S_o \in [0,1]
   %
   %    This is 1-lowest(Sg) at which the same criterion is satisfied.
   %
   i     = find(~(T(:,3) > 0), 1, 'first');
   Sogcr = 1 - T(i, 1);
end

%--------------------------------------------------------------------------

% Gas relative permeability as a function of gas saturation, sg.
% Note: Gas saturation is always <= 1-Swco.
function varargout = interp_gas(sg, T, Swco)
   sg = min(sg, 1 - Swco);

   varargout{1}    = interpTable(T(:,1), T(:,2), sg);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,2), sg);
   end
end

%--------------------------------------------------------------------------

% Oil relative permeability as a function of oil saturation, so.
% Note: Oil saturation is always <= 1-Swco.
function varargout = interp_oil(so, T, Swco)
   sg = max(1 - Swco - so, 0);

   varargout{1}    =   interpTable(T(:,1), T(:,3), sg);

   if nargout > 1,
      varargout{2} = -dinterpTable(T(:,1), T(:,3), sg);
   end
end

%--------------------------------------------------------------------------

% Oil-gas capillary pressure and derivative as a function of gas saturation
function varargout = interp_pc(sg, T, Swco)
   sg = min(sg, 1 - Swco);

   varargout{1}    =  interpTable(T(:,1), T(:,4), sg);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,4), sg);
   end
end

%--------------------------------------------------------------------------

% Inverse of capillary pressure function
function varargout = interp_pcinv(dp, T)
   [i, i, j]    = unique(T(:,4));                                      %#ok
   varargout{1} = interpTable(T(i,4), T(i,1), dp);
end

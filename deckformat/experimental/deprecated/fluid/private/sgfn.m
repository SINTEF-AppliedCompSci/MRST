function [krg, Sgco, Sgcr, Sgmax, pc, pcinv] = sgfn(T, varargin)
%Construct gas relperm evaluation functions from SGFN table.
%
% SYNOPSIS:
%   [krg, Sgco, Sgcr, Sgmax, pc, pcinv] = sgfn(table)
%   [krg, Sgco, Sgcr, Sgmax, pc, pcinv] = sgfn(table, Swco)
%
% PARAMETERS:
%   table - Saturation function table as defined through the input deck
%           keyword 'SGFN'.  May be a numeric array or a a cell array of
%           such arrays.
%
%   Swco  - Connate water saturation.  One scalar value (in [0,1]) for each
%           saturation region.  Applicable to water-oil-gas systems.
%           Default value: Swco = 0 (no connate water in system).
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
%           parameters, i.e., [krg{1:2}] = krg(sg), then the second output
%           parameter is the derivative of the gas relative permeability
%           with respect to the *gas* saturation at 'sg'.
%
%   Sgco  - A 1-by-ntab numeric array of connate gas saturations.
%
%   Sgcr  - A 1-by-ntab numeric array of critical gas saturations.  This
%           is the largest gas saturation for which krg(sg)=0.
%
%   Sgmax - A 1-by-ntab numeric array of maximum gas saturations.
%
%   pc    - A 1-by-ntab cell array of function handles such that the call
%           'pc(sw)' computes the capillary pressure. When called with two
%           output parameters, i.e., [pc{1:2}] = pc(sw), then the second
%           output parameter is the derivative of capillary pressure.
%
%   pcinv - A 1-by-ntab cell array of function handles such that the call
%           'pcinv(dp)' computes the inverse of the capillary pressure.
%
% SEE ALSO:
%  `readRelPermTable`, `swfn`, `sgofn`, `swof`.

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

   krg   = cell ([1, ntab]);
   pc    = cell ([1, ntab]);
   pcinv = cell ([1, ntab]);
   Sgcr  = zeros([1, ntab]);
   Sgmax = zeros([1, ntab]);
   Sgco  = zeros([1, ntab]);

   % Set connate water saturation (if water is present).
   Swco = zeros([1, ntab]);
   if nargin > 1 && isnumeric(varargin{1}) && ...
         any(numel(varargin{1}) == [1, ntab]),

      Swco = varargin{1};

      assert (all(~ ((Swco < 0) | (Swco > 1)) ));
   end

   if numel(Swco) == 1, Swco = repmat(Swco, [1, ntab]); end

   for k = 1 : ntab,
      [Sgmax(k), Sgcr(k), Sgco(k)] = table_points(T{k}, Swco(k));

      krg  {k} = @(sg) interp_gas  (sg, T{k}, Swco(k));
      pc   {k} = @(sg) interp_pc   (sg, T{k}, Swco(k));
      pcinv{k} = @(dp) interp_pcinv(dp, T{k});
   end
end

%--------------------------------------------------------------------------

function [Sgmax, Sgcr, Sgco] = table_points(T, Swco)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:2) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');

   assert (~ (T(1, 2) > 0), 'krg(Sgco) must be zero (SGFN).');
   assert (all (diff(T(:,1)) > 0), ...
          ['Gas saturation must be strictly increasing down ', ...
           'column (SGFN).']);
   assert (~any(diff(T(:,2)) < 0), ...
          ['Gas relative permeability must level or increase ', ...
           'down column (SGFN).']);

   %-----------------------------------------------------------------------
   % Data OK.  Extract critical values from table.
   %

   % 1) Connate gas saturation
   %
   Sgco = min(T(:,1));  assert (~ (Sgco < 0));

   % 2) Maximal gas saturation
   %
   Sgmax = max(T(:,1));

   if     Sgmax > 1 - Swco,

      warning(msgid('Sgmax:NoSwco'), ...
             ['Maximum gas saturation in ''SGFN'' table (=%f) does ', ...
              'not account for presence of connate water.'], Sgmax);

   elseif Sgmax < 1 - Swco,

      warning(msgid('Sgmaxt:LessSwco'), ...
             ['Maximum gas saturation in ''SGFN'' table (=%f) is ',  ...
              'less than connate water saturation.\n',                ...
              'The maximum gas saturation should normally be 1-Swco', ...
              ' (=%f).'], Sgmax, 1 - Swco);

   end

   % 3) Critical gas saturation (>= 0) is highest Sg for which krg(Sg) = 0.
   %
   i    = find(T(:,2) > 0, 1, 'first');  assert (~isempty(i) && (i > 1));
   Sgcr = T(i - 1, 1);                   assert (~ (Sgcr < 0));
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

% Oil-gas capillary pressure and derivative as a function of gas saturation
function varargout = interp_pc(sg, T, Swco)
   sg = min(sg, 1 - Swco);

   varargout{1} = interpTable(T(:,1), T(:,3), sg);

   if nargout > 1,
      varargout{2} = dinterpTable(T(:,1), T(:,3), sg);
   end
end

%--------------------------------------------------------------------------

% Inverse of capillary pressure function
function varargout = interp_pcinv(dp, T)
   [i, i] = unique(T(:,3));                                            %#ok
   varargout{1} = interpTable(T(i,3), T(i,1), dp);
end

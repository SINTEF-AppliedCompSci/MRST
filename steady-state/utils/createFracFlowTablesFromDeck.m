function T = createFracFlowTablesFromDeck(deck)
%Undocumented Utility Function

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

assert(isfield(deck, 'PROPS'), 'Invalid deck. Missing ''PROPS'' field.');
assert(isfield(deck.PROPS, 'SWOF'), 'Missing property ''SWOF''.');
assert(isfield(deck.PROPS, 'PVTW'), 'Missing property ''PVTW''.');

% PVTW
% Columns: [REF.PRES., REF.FVF, COMPRESSIBILITY, REF.VISC., VISCOSIBILITY]
muW   = deck.PROPS.PVTW(4); % water viscosity
if deck.PROPS.PVTW(5) > 0
   warning(['There is nonzero water viscosibility. Only the reference '...
      'viscosity is used.']);
end

% PVCDO
if isfield(deck.PROPS, 'PVCDO')
   % Columns: [REF.PRES., REF.FVF, COMPRESSIBILITY, REF.VISC., VISCOSIBILITY]
   muO   = deck.PROPS.PVCDO(4); % oil viscosity
   if deck.PROPS.PVCDO(5) > 0
      warning(['There is nonzero water viscosibility. '...
          'Only the reference viscosity is used.']);
   end
elseif isfield(deck.PROPS, 'PVDO')
    muO = deck.PROPS.PVDO{1}(1,3); % oil viscosity
elseif isfield(fluid, 'muO')
   %muO = fluid.muO(100*barsa
   error('Have to be implemented')
end

% SWOF
% Columns: [Sw, Krw, Krow, Pcow]
% Fractional flow = (krW/muW) / ( (krW/muW) + (krO/muO) )
T  = cellfun(@(x) [( x(:,2)/muW )./( (x(:,2)/muW) + (x(:,3)/muO) ), ...
   x(:,1)], deck.PROPS.SWOF, 'UniformOutput', false);

end

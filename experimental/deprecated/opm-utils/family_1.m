function f1 = family_1(family_2)
%Convert ECLIPSE sat-func family II (SWFN &c) to family I (SWOF &c)
%
% SYNOPSIS:
%   f1 = familiy_1(f2)
%
% PARAMETERS:
%   f2 - Structure describing relative permeability (and capillary
%        pressure) functions using "family II" keywords (i.e., SWFN, SGFN,
%        SOF2 and/or SOF3).  Indidvidual fields corresponding to the actual
%        keyword data.  Fields may be cell arrays of tables to support
%        multiple regions (e.g., SATNUM and/or IMBNUM).
%
%        Typically corresponds to the 'deck.PROPS' sub-structure input
%        using the 'readEclipseDeck' function from the 'deckformat' module.
%
% RETURNS:
%   f1 - Structure describing the same saturation functions using
%        "family I" keywords (i.e., SWOF and/or SGOF).
%
% SEE ALSO:
%   readEclipseDeck, writeDeck.

%{
Copyright 2009-2016 SINTEF ICT, Applied Mathematics.

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

% $Date: 2012-10-03 15:11:46 +0200 (Wed, 03 Oct 2012) $
% $Revision: 9993 $

   wat = isfield(family_2, 'SWFN');
   gas = isfield(family_2, 'SGFN');

   assert (wat || gas, ...
           'Must specify at least WATER/OIL or GAS/OIL system.');

   if xor(wat, gas),
      % Two-phase system
      assert (isfield(family_2, 'SOF2'), ...
              'Must specify ''SOF2'' in two-phase systems.');
   else
      % Three-phase system
      assert (isfield(family_2, 'SOF3'), ...
              'Must specify ''SOF3'' in three-phase systems.');
   end

   if wat, f1.SWOF = wat_oil(family_2, gas); end
   if gas, f1.SGOF = gas_oil(family_2, wat); end
end

%--------------------------------------------------------------------------

function swof = wat_oil(f2, gas)
   if ~ iscell(f2.SWFN), f2.SWFN = { f2.SWFN }; end

   ntab = numel(f2.SWFN);
   swof = cell([1, ntab]);

   swfn = f2.SWFN;
   if gas,
      if ~ iscell(f2.SOF3), f2.SOF3 = { f2.SOF3 }; end

      sof2 = cellfun(@(t) t(:,1:2), f2.SOF3, 'UniformOutput', false);
   else
      if ~ iscell(f2.SOF2), f2.SOF2 = { f2.SOF2 }; end

      sof2 = f2.SOF2;
   end

   assert (numel(swfn) == numel(sof2));

   for k = 1 : ntab,
      swof{k} = wat_oil_kernel(swfn{k}, sof2{k});
   end

   if ntab == 1, swof = swof{1}; end
end

%--------------------------------------------------------------------------

function sgof = gas_oil(f2, wat)
   if ~ iscell(f2.SGFN), f2.SGFN = { f2.SGFN }; end

   ntab = numel(f2.SGFN);
   sgof = cell([1, ntab]);

   sgfn = f2.SGFN;
   if wat,
      if ~ iscell(f2.SOF3), f2.SOF3 = { f2.SOF3 }; end
      if ~ iscell(f2.SWFN), f2.SWFN = { f2.SWFN }; end

      assert (numel(f2.SWFN) == numel(f2.SOF3));

      sof2 = cellfun(@(t) t(:,[1,3]), f2.SOF3, 'UniformOutput', false);
      swco = cellfun(@(t) t(1,1)    , f2.SWFN);
   else
      if ~ iscell(f2.SOF3), f2.SOF3 = { f2.SOF3 }; end

      sof2 = f2.SOF2;
      swco = zeros(size(f2.SOF2));
   end

   assert (numel(sgfn) == numel(sof2));
   assert (numel(sgfn) == numel(swco));

   for k = 1 : ntab,
      sgof{k} = gas_oil_kernel(sgfn{k}, sof2{k}, swco(k));
   end

   if ntab == 1, sgof = sgof{1}; end
end

%--------------------------------------------------------------------------

function swof = wat_oil_kernel(swfn, sof2)
   gpp = @(x,y) interp1(x, y, 'linear', 'pp');

   ppw = gpp(swfn(:,1), swfn(:, 2:end));
   ppo = gpp(flipud(1 - sof2(:,1)), flipud(sof2(:,2)));

   sw   = usat(swfn(:,1), sof2(:,1), 1e9);
   swof = ppval(ppw, sw);
   krow = ppval(ppo, sw);

   swof = [sw, swof(:,1), krow, swof(:,2)];
end

%--------------------------------------------------------------------------

function sgof = gas_oil_kernel(sgfn, sof2, swco)
   gpp = @(x,y) interp1(x, y, 'linear', 'pp');

   ppg = gpp(sgfn(:,1), sgfn(:,2:end));
   ppo = gpp(flipud(1 - (sof2(:,1) + swco)), flipud(sof2(:,2)));

   sg   = usat(sgfn(:,1), sof2(:,1) + swco, 1e9);
   sgof = ppval(ppg, sg);
   krog = ppval(ppo, sg);

   sgof = [sg, sgof(:,1), krog, sgof(:,2)];
end

%--------------------------------------------------------------------------

function s = usat(s1, s2, adj)
   s = unique(fix(adj * [s1 ; 1 - s2]) / adj);
end

function [pvtfuns, surface_density, names, info, incomp] = ...
      eclipsePhaseProperties(deck, varargin)
%Construct PVT evaluators for each active phase of an ECLIPSE/FrontSim deck
%
% SYNOPSIS:
%   [pvtfuns, rho_s, names, info, incomp] = eclipsePhaseProperties(deck)
%
% PARAMETERS:
%   deck    - An ECLIPSE/FrontSim input deck as defined by function
%             readEclipseDeck.
%
% RETURNS:
%   pvtfuns - A three-element cell array of (cell arrays of) function
%             handles.  The elements of the return value correspond to the
%             PVT/property functions of
%
%                 Aqua (Water) -- pvtfuns{1}
%                 Liquid (Oil) -- pvtfuns{2}
%                 Vapour (Gas) -- pvtfuns{3}
%
%             As a special case, the i'th element of 'pvtfuns' is empty
%             (i.e., ISEMPTY(pvtfuns{i})==TRUE) if the corresponding phase
%             is not declared active in the input 'deck'.
%
%   rho_s   - An m-by-np matrix of surface densities as specified by
%             DENSITY keyword.  Here, 'm' is the number of PVT regions
%             (often one) and 'np' is the number of declared phases (often
%             two or three).
%
%   names   - Cell array of names of fluid components (WATER, OIL, GAS)
%             declared in deck.
%
%   info    - Human readable summary of fluid properties.
%
%   incomp  - Boolean vector specifying incomressible components.  Water is
%             considered incompressible if the compressibility is zero.
%             Similarly, oil is considered incompressible if all rows of
%             columns 2:3 in PVDO are identical or if zero compressibility
%             is specified in PVCDO.  Gas is considered incompressible if
%             all rows of columns 2:3 in PVDG are identical.
%
% SEE ALSO:
%   `readEclipseDeck`.

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

   rspec = deck.RUNSPEC;
   props = deck.PROPS;

   verify_pvt_specification(rspec, props);

   ntpvt = 1;
   if isfield(rspec, 'TABDIMS')
      ntpvt = rspec.TABDIMS(2);
   end

   % Hard-wired phase sequence:
   %     Aqua (Water) <-> 1
   %     Liquid (Oil) <-> 2
   %     Vapour (Gas) <-> 3
   %
   [A, L, V] = deal(1, 2, 3);

   % Determine active phases in run/input deck.
   names = { 'water', 'oil', 'gas' };
   phase = isfield(rspec, upper(names));
   names = names(phase);

   pos   = cumsum(double(phase));

   pvtfuns         = repmat({ cell([1, ntpvt]) }, [1, 3]);
   miscible        = false([1, 3]);
   incomp          = false([1, 3]);
   surface_density = zeros([ntpvt, pos(end)]);

   rhoS = extract_surface_density(props);

   if phase(A)
      [pvtfuns{A}, miscible(A), incomp(A)] = ...
         wat_functions(ntpvt, props);

      surface_density(:, pos(A)) = rhoS.water;
   end

   if phase(L)
      [pvtfuns{L}, miscible(L), incomp(L)] = ...
         oil_functions(ntpvt, props, pos, L, V);

      surface_density(:, pos(L)) = rhoS.oil;
   end

   if phase(V)
      [pvtfuns{V}, miscible(V), incomp(V)] = ...
         gas_functions(ntpvt, props, pos, L, V);

      surface_density(:, pos(V)) = rhoS.gas;
   end

   % Miscible data requires presence of both solvent and solute.
   if ~(miscible(L) == phase(V) || phase(V))
      error('Cannot specify miscible liquid phase in absence of vapour.');
   end
   if ~(miscible(V) == phase(L) || phase(L))
      error('Cannot specify miscible vapour phase in absence of liquid.');
   end
   if ~(miscible(A) == phase(V) || phase(V))
      error('Cannot specify miscible aquaic phase in absence of vapour.');
   end

   info = summarise(pos, surface_density, phase, ...
                    miscible, incomp, A, L, V);

   incomp = incomp(phase);
end

%--------------------------------------------------------------------------

function verify_pvt_specification(rspec, props)
% Ensure minimally consistent PVT specifcation in run
%
% Rules:
%
%   Run declares        Then the following conditions must be met
%   ----------------    ---------------------------------------------
%   Dissolved gas    => Both OIL and GAS are active phases
%                       OIL PVT properties entered using PVTO or PVCO
%                       Rs cannot be constant (RSCONST, RSCONSTT)
%
%   Vapourised oil   => Both OIL and GAS are active phases
%                       GAS PVT properties entered using PVTG
%                       Rv cannot be constant (RVCONST, RVCONSTT)
%
%   PVTO, PVCO       => Both OIL and GAS are active phases
%                       Rs cannot be constant (RSCONST, RSCONSTT)
%
%   PVTG             => Both OIL and GAS are active phases
%                       Rv cannot be constant (RVCONST, RVCONSTT)
%
%   PVDO, PVDG       => OIL, GAS active phases, respectively
%
%   PVDO + const. Rs => GAS is not active
%
%   PVDG + const. Rv => OIL is not active

   if isfield(rspec, 'DISGAS'),
      assert (all(isfield(rspec, {'OIL', 'GAS'})),             ...
             ['Dissolved gas (keyword ''DISGAS'') cannot be ', ...
              'declared unless run declares active OIL and GAS phases.']);

      assert (sum(isfield(props, {'PVTO', 'PVCO'})) == 1,            ...
             ['Dissolved gas declaration must be accompanied by a ', ...
              'PVTO or PVCO table.']);

      assert (~any(isfield(props, {'RSCONST', 'RSCONSTT'}))      , ...
             ['Dissolved gas implies an active gas phase.  The ' , ...
              'dissolved gas/oil ratio (Rs) therefore cannot be ', ...
              'constant as implied by ''RSCONST'' or ''RSCONSTT''.']);
   end

   if isfield(rspec, 'VAPOIL'),
      assert (all(isfield(rspec, {'OIL', 'GAS'})),              ...
             ['Vapourised oil (keyword ''VAPOIL'') cannot be ', ...
              'declared unless run declares active OIL and GAS phases.']);

      assert (isfield(props, 'PVTG'), ...
             ['Vapourised oil declaration must be accompanied by a ', ...
              'PVTG table.']);

      assert (~any(isfield(props, {'RVCONST', 'RVCONSTT'}))       , ...
             ['Vapourised oil implies an active oil phase.  The ' , ...
              'Vapourised oil/gas ratio (Rv) therefore cannot be ', ...
              'constant as implied by ''RVCONST'' or ''RVCONSTT''.']);
   end

   if any(isfield(props, {'PVTO', 'PVCO'})),
      assert (all(isfield(rspec, {'OIL', 'GAS'})), ...
             ['Live oil properties cannot be defined unless ', ...
              'both OIL and GAS are declared active phases.']);

      assert (~any(isfield(props, {'RSCONST', 'RSCONSTT'})), ...
             ['live oil models preclude constant dissolved ', ...
              'gas/oil ratio (Rs).']);
   end

   if isfield(props, 'PVTG'),
      assert (all(isfield(rspec, {'OIL', 'GAS'})), ...
             ['Wet gas properties cannot be defined unless ', ...
              'both OIL and GAS are declared active phases.']);

      assert (~any(isfield(props, {'RVCONST', 'RVCONSTT'})), ...
             ['Wet gas models preclude constant vapourised ', ...
              'oil/gas ratio (Rv).']);
   end

   if isfield(props, 'PVDO'),
      assert (isfield(rspec, 'OIL'), ...
              'Dead oil runs must define an active OIL phase.');

      if any(isfield(props, {'RSCONST', 'RSCONSTT'})),
         assert (~isfield(rspec, 'GAS'), ...
                ['Dead oil runs with constant dissolved gas/oil ', ...
                 'ratio cannot declare an active GAS phase.']);
      end
   end

   if isfield(props, 'PVDG'),
      assert (isfield(rspec, 'GAS'), ...
              'Dry gas runs must declare an active GAS phase.');

      if any(isfield(props, {'RVCONST', 'RVCONSTT'})),
         assert (~isfield(rspec, 'OIL'), ...
                ['Dry gas runs with constant vapourised oil/gas ', ...
                 'ratio cannot declare an active GAS phase.']);
      end
   end
end

%--------------------------------------------------------------------------

function [f, miscible, incomp] = oil_functions(ntab, props, pos, L, V)
   kws = { 'PVCDO', 'PVDO', 'PVTO' };
   fn  = fieldnames(props);

   [i, j]     = blockDiagIndex(numel(fn), numel(kws));
   tmp        = [reshape(fn(i), [], 1), reshape(kws(j), [], 1)];
   kw_present = accumarray(j, strcmp(tmp(:,1), tmp(:,2))) > 0;

   nkw = sum(kw_present);

   if nkw == 0,
      error('Liquid properties specified through unsupported keyword.');
   elseif nkw > 1,
      error('Liquid properties specified through more than one keyword.');
   end

   f = cell([1, ntab]);
   if isfield(props, 'PVCDO'),

      assert (size(props.PVCDO, 1) == ntab);

      for t = 1 : ntab,
         f{t} = @(p, z) pvcdo(props.PVCDO(t,:), p);
      end

      miscible = false;
      incomp   = ~any(abs(props.PVCDO(:,3)) > 0);

   elseif isfield(props, 'PVDO'),

      assert (numel(props.PVDO) == ntab);

      incomp = true;
      for t = 1 : ntab,
         f{t} = @(p, z) pvdx(props.PVDO{t}, p);

         % Constant 'B' or single line table
         incomp = incomp && (numel(rlencode(props.PVDO{t}(:,2:3))) == 2);
      end

      miscible = false;

   else

      % PVTO as verified by 'nkw' tests above.
      assert (numel(props.PVTO) == ntab);

      for t = 1 : ntab,
         tab = [nan([size(props.PVTO{t}.data, 1), 1]), props.PVTO{t}.data];
         tab(props.PVTO{t}.pos(1 : end - 1), 1) = props.PVTO{t}.key;

         f{t} = @(p, z) pvto(tab, p, z(:, pos(V)) ./ (z(:, pos(L))+eps));
      end

      miscible = true;
      incomp   = false;

   end

   if ntab == 1,
      % Optimize for common case (single PVT region).
      f = f{1};
   end
end

%--------------------------------------------------------------------------

function [f, miscible, incomp] = gas_functions(ntab, props, pos, L, V)
   kws = { 'PVDG', 'PVTG' };
   fn  = fieldnames(props);

   [i, j]     = blockDiagIndex(numel(fn), numel(kws));
   tmp        = [reshape(fn(i), [], 1), reshape(kws(j), [], 1)];
   kw_present = accumarray(j, strcmp(tmp(:,1), tmp(:,2))) > 0;

   nkw = sum(kw_present);

   if nkw == 0,
      error('Vapour properties specified through unsupported keyword.');
   elseif nkw > 1,
      error('Vapour properties specified through more than one keyword.');
   end

   f = cell([1, ntab]);
   if isfield(props, 'PVDG'),

      assert (numel(props.PVDG) == ntab);

      incomp = true;
      for t = 1 : ntab,
         f{t} = @(p, z) pvdx(props.PVDG{t}, p);

         % Constant 'B' or single line table
         incomp = incomp && (numel(rlencode(props.PVDG{t}(:,2:3))) == 2);
      end

      miscible = false;

   else

      % PVTG as verified by 'nkw' tests above.
      assert (numel(props.PVTG) == ntab);

      for t = 1 : ntab,
         tab = [nan([size(props.PVTG{t}.data, 1), 1]), props.PVTG{t}.data];
         tab(props.PVTG{t}.pos(1 : end - 1), 1) = props.PVTG{t}.key;

         f{t} = @(p, z) pvtg(tab, p, z(:, pos(L)) ./ (z(:, pos(V))+eps));
      end

      miscible = true;
      incomp   = false;

   end

   if ntab == 1,
      % Optimize for common case (single PVT region).
      f = f{1};
   end
end

%--------------------------------------------------------------------------

function [f, miscible, incomp] = wat_functions(ntab, props)
   if ~isfield(props, 'PVTW'),
      error('Water keywords other than ''PVTW'' currently not supported.');
   end

   assert (size(props.PVTW, 1) == ntab);

   f = cell([1, ntab]);
   for t = 1 : ntab,
      f{t} = @(p, z) pvtw(props.PVTW(t,:), p);
   end
   miscible = false;
   incomp   = ~any(abs(props.PVTW(:,3)) > 0);

   if ntab == 1,
      % Optimize for common case (single PVT region).
      f = f{1};
   end
end

%--------------------------------------------------------------------------

function rhoS = extract_surface_density(props)

   if isfield(props, 'DENSITY')
      rhoS = extract_DENSITY(props);
   elseif isfield(props, 'GRAVITY')
      rhoS = extract_GRAVITY(props);
   else
      error('Density:Missing', 'Fluid Density Data Misssing');
   end
end

%--------------------------------------------------------------------------

function rhoS = extract_DENSITY(props)
   rhoS = struct('water', props.DENSITY(:,2), ...
                 'oil',   props.DENSITY(:,1), ...
                 'gas',   props.DENSITY(:,3));
end

%--------------------------------------------------------------------------

function rhoS = extract_GRAVITY(props)
   rho0 = struct('wat', 1000*kilogram/meter^3, ...
                 'air', 1.22*kilogram/meter^3);  % E100 reference values

   rho_oil = rho0.wat .* (141.5 ./ (props.GRAVITY(:,1) + 131.5));

   rhoS = struct('water', rho0.wat .* props.GRAVITY(:, 2), ...
                 'oil',   rho_oil, ...
                 'gas',   rho0.gas .* props.GRAVITY(:, 3));
end

%--------------------------------------------------------------------------

function info = summarise(pos, rhoS, phase, miscible, incomp, A, L, V)
   info = '';

   if phase(A),
      info = [ info, densstring('water', '', rhoS(:, pos(A))) ];
   end

   if phase(L),
      info = [ info, densstring('oil', '  ', rhoS(:, pos(L))) ];
   end

   if phase(V),
      info = [ info, densstring('gas', '  ', rhoS(:, pos(V))) ];
   end

   info = [ info, newline ];

   if phase(A),
      info = [ info, mcstring('Aquaic', miscible(A), incomp(A)) ];
   end

   if phase(L),
      info = [ info, mcstring('Liquid', miscible(L), incomp(L)) ];
   end

   if phase(V),
      info = [ info, mcstring('Vapour', miscible(V), incomp(V)) ];
   end
end

%--------------------------------------------------------------------------

function sr = densstring(ph, blank, rhoS)
   sr = sprintf(' %8.2f', reshape(rhoS, 1, []));
   sr = [ 'Surface density for ', ph, ':', blank, sr, newline ];
end

%--------------------------------------------------------------------------

function mc = mcstring(ph, misc, incomp)
   mc = [ ph, ' phase:  ' ];

   if misc,   mc = [ mc, 'miscible'  ];
   else       mc = [ mc, 'immiscible']; end

   mc = [ mc, ' and ' ];

   if incomp, mc = [ mc, 'incompressible' ];
   else       mc = [ mc, 'compressible'   ]; end

   mc = [ mc, newline ];
end

%--------------------------------------------------------------------------

function s = newline
   s = sprintf('\n');
end

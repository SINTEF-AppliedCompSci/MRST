function [pvtfun, surface_density, info, names] = processBlackOilDec(deck, varargin)
%Construct Black Oil PVT model evaluation function.
%

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



   error('You should consider joining the new decade...');

   % Hard-wired sequence, aqua (water), liquid (oil), vapor (gas)
   A = 1;   L = 2;    V = 3;


   check = @(deck, f) isfield(deck,f);


   % Declaration of which phases are present.
   phase    = [isfield(deck.RUNSPEC, 'WATER'), ...
               isfield(deck.RUNSPEC, 'OIL'), ...
               isfield(deck.RUNSPEC, 'GAS')];

   % User readable phase names
   names = {'water','oil','gas'};
   names = names(phase);


   incomp   = false(3, 1);
   miscible = false(3, 1);
   pvtfun   = cell(3,1);


   fprintf('Using keywords\n');


   %% Water phase pvt function.
   if check(deck.PROPS, 'PVTW'),
      if ~phase(A),
         error(id('WaterPhaseData:PhaseNotPresent'), ...
            'Data for water phase specified even though water is not present.');
      elseif isempty(pvtfun{A}),
         fprintf(' PVTW\n');
         incomp(A) = size(deck.PROPS.PVTW(1,:), 1) == 1;
         pvtfun{A} = @(varargin)pvtw(deck.PROPS.PVTW(1,:), varargin{:});
      else
         error(id('WaterPhaseData:AlreadySet'), ...
            'Data for water phase specified more than once.\n');
      end
   end

   %% Oil phase pvt function
   if check(deck.PROPS,'PVTO'),
      if ~phase(L),
         error(id('WaterPhaseData:PhaseNotPresent'), ...
            'Data for oil phase specified even though oil is not present.');

      elseif isempty(pvtfun{L}),
         fprintf(' PVTO\n');
         t=1;%use first table
         props=deck.PROPS;
         tab = [nan([size(props.PVTO{t}.data, 1), 1]), props.PVTO{t}.data];
         tab(props.PVTO{t}.pos(1 : end - 1), 1) = props.PVTO{t}.key;

         pvtfun{L}=@(varargin) pvto(tab,varargin {:});
         miscible(L) = true;
      else
         error('Data for liquid phase specified more than once.');
      end

   elseif check(deck.PROPS,'PVDO'),
      if ~phase(L),
         error(id('LiquidPhaseData:PhaseNotPresent'), ...
            'Data for oil phase specified even though oil is not in present.');
      elseif isempty(pvtfun{L}),
         fprintf(' PVDO\n');
         incomp(L) = size(deck.PROPS.PVDO{1}, 1) == 1;
         pvtfun{L}=@(varargin) pvdx(deck.PROPS.PVDO{1}, varargin{:});
      else
         error('Data for liquid phase specified more than once.');
      end
   elseif check(deck.PROPS,'PVCDO'),
      if ~phase(L),
         error(id('LiquidPhaseData:PhaseNotPresent'), ...
            'Data for oil phase specified even though oil is not in present.');
      elseif isempty(pvtfun{L}),
         fprintf(' PVCDO\n');
         pvtfun{L}=@(varargin) pvcdo(deck.PROPS.PVCDO(1,:), varargin{1});
      else
         error('Data for liquid phase specified more than once.');
      end
   end

   %% Gas phase pvt function
   if check(deck.PROPS,'PVTG'),
      if ~phase(V),
         error(id('VaporPhaseData:PhaseNotPresent'), ...
            'Data for vapor phase specified even though vapor is not present.');
      elseif isempty(pvtfun{V}),
         fprintf(' PVTG\n');
         t=1;%use first table
         props=deck.PROPS;
         tab = [nan([size(props.PVTG{t}.data, 1), 1]), props.PVTG{t}.data];
         tab(props.PVTG{t}.pos(1 : end - 1), 1) = props.PVTG{t}.key;

         pvtfun{V}=@(varargin)pvtg(tab, varargin{:});
         miscible(V) = true;
      else
         error('Data for vapor phase specified more than once.');
      end

   elseif check(deck.PROPS,'PVDG'),
      if ~phase(V),
         error(id('VaporPhaseData:PhaseNotPresent'), ...
            'Data for vapor phase specified even though vapor is not in present.');
      elseif isempty(pvtfun{V}),
         fprintf(' PVDG\n');
         incomp(V) = size(deck.PROPS.PVDG{1}, 1) == 1;
         pvtfun{V}=@(varargin) pvdx(deck.PROPS.PVDG{1}, varargin{:});
      else
         error('Data for vapor phase specified more than once.');
      end
   end


   if sum(phase) < 1,
      error('Huh!?, No phases?');
   end


   %% Micible data requires two phases to be present...
   if ~(miscible(L) == phase(V) || phase(V)),
      error('Cannot specify miscible liquid phase without vapor phase present');
   end
   if ~(miscible(V) == phase(L) || phase(L)),
      error('Cannot specify miscible vapor phase without liquid phase present');
   end

   if ~(miscible(A) == phase(V) || phase(V)),
      error('Cannot specify miscible aquaic phase without vapor phase present');
   end

   %%

   if ~isempty(pvtfun{A}), pvtfun{A} = @(p,z) pvtfun{A}(p); end
   if ~isempty(pvtfun{L}), pvtfun{L} = @(p,z) pvtfun{L}(p, z(:,V)./z(:,L)); end
   if ~isempty(pvtfun{V}), pvtfun{V} = @(p,z) pvtfun{V}(p, z(:,L)./z(:,V)); end

   pos = cumsum(phase);
   for i=1:numel(phase),
      if phase(i),
         assert(~isempty(pvtfun{i}), ...
         'Invalid deck specified: PVT data not specified for the %s phase.',...
          names{pos(i)});
      end
   end
   clear i pos
   %% Surface densities

   pos = cumsum(phase);
   if check(deck.PROPS, 'DENSITY'), d = deck.PROPS.DENSITY(1,:);
      assert(numel(d)==3);
      if phase(A), surface_density(pos(A)) = d(2);end
      if phase(L), surface_density(pos(L)) = d(1);end
      if phase(V), surface_density(pos(V)) = d(3);end
   else
      error (['PVT table must contain keyword "density" with one table ',...
         'per pvt region']);
   end

   %%
   info = '';
   if phase(A), info = [info,sprintf('Surface density for water: %8.2f\n',...
         surface_density(pos(A)))];
   end
   if phase(L), info = [info,sprintf('Surface density for oil:   %8.2f\n',...
         surface_density(pos(L)))];
   end
   if phase(V), info = [info,sprintf('Surface density for gas:   %8.2f\n\n',...
         surface_density(pos(V)))];
   end

   if phase(A),
      if miscible(A),info = [info,sprintf('Aquaic phase:  miscible ')];
      else           info = [info,sprintf('Aquaic phase:  immiscible ')]; end
      if incomp(A),  info = [info,sprintf('and incompressible\n')];
      else           info = [info,sprintf('and compressible\n')]; end
   end


   if phase(L),
      if miscible(L),info = [info,sprintf('Liquid phase:  miscible ')];
      else           info = [info,sprintf('Liquid phase:  immiscible ')]; end
      if incomp(L),  info = [info,sprintf('and incompressible\n')];
      else           info = [info,sprintf('and compressible\n')]; end
   end

   if phase(V),
      if miscible(V),info = [info,sprintf('Vapour phase:  miscible ')];
      else           info = [info,sprintf('Vapour phase:  immiscible ')]; end
      if incomp(V),  info = [info,sprintf('and incompressible\n')];
      else           info = [info,sprintf('and compressible\n')]; end
   end
   disp(' ');
   disp(info);
end

function s = id(s)
   s = ['processBlackOilDeck:', s];
end

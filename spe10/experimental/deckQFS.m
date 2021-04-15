function  deck=deckQFS(cartdims,physdims,varargin)
%Undocumented Utility Function

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

opt = struct('perm', 100,...
             'poro', 0.1,...
             'permcase','uniform',...
             'use_pressure',true);
opt = merge_options(opt, varargin{:});
if(nargin==1)
   myperm='spe10_layer';
   layers=cartdims;
   cartdims = [  60,  220,   numel(layers)];
   %physdims = [1200, 2200, 2*cartdims(end)] .* ft();
   physdims = cartdims.*[20,10,2]*ft;
else
   layers=1:cartdims(3);
   if(strcmp(opt.permcase,'uniform'))
      myperm='uniform';
      perm = opt.perm;
      poro = opt.poro;
   else
      disp('Try to use spe10 perm and poro');
      error('Case still to be defined');
   end
end
nc=prod(cartdims);
% define runspec
deck.RUNSPEC.cartDims=cartdims;
deck.RUNSPEC.DIMENS=cartdims;
deck.RUNSPEC.OIL=1;
deck.RUNSPEC.WATER=1;
%deck.RUNSPEC.FIELD=1;
deck.RUNSPEC.METRIC=1;
deck.RUNSPEC.TABDIMS=[1     1    20    50    20    50     1    20    20     1    10     1    -1     0     1];
deck.RUNSPEC.WELLDIMS=[5 10 2 1 5 10 5 4 3 0 1 1];
deck.RUNSPEC.AQUDIMS=[0 0 0 0 10 10 0 0];
deck.RUNSPEC.START=datenum('01 Jan 2010');
deck.RUNSPEC.NOGRAV=1;
%define solution solution
pres=200;
sat=1;
deck.SOLUTION.PRESSURE=repmat(pres,[nc, 1]);
deck.SOLUTION.SWAT=repmat(sat,[nc, 1]);
deck.SOLUTION.SOIL=repmat(1-sat,[nc, 1]);
% define REGIONS
deck.REGIONS.SATNUM=ones(nc,1);
%define props
s=linspace(0,1,2)';alpha=2;
deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
deck.PROPS.PVTW=[0 1 0 0.5 0];
deck.PROPS.PVCDO=[0 1 0 5 0];

%deck.PROPS.PVDO{1}=[0 1 2;8000 1 2;16000 1 2];
deck.PROPS.ROCK=[4000 3.000000000000000e-06 NaN NaN NaN NaN];
deck.PROPS.DENSITY=[52 64 0.044000000000000];
deck.PROPS.PVTW=[0 1 0 0.500000000000000 0];
% define summary
deck.SUMMARY=[];
% difine SC
%%
%refdepth = 4000*ft;
refdepth = 0 + physdims(3)/(2.0*cartdims(3));
deck.SCHEDULE.control.WELSPECS = ...
   { ...
   'I01', 'W', 1          , 1          , refdepth, 'OIL', 0, 'STD', 'SHUT', 'NO', 0, 'SEG', 0; ...
   'P01', 'W', cartdims(1), cartdims(2), refdepth, 'OIL', 0, 'STD', 'SHUT', 'NO', 0, 'SEG', 0; ...
   };

deck.SCHEDULE.control.COMPDAT = ...
   { ...
   'I01', 1          , 1          , 1, cartdims(3), 'OPEN', 0, 0, 0.125, -1, 0, 'Default', 'Z', -1; ...
   'P01', cartdims(1), cartdims(2), 1, cartdims(3), 'OPEN', 0, 0, 0.125, -1, 0, 'Default', 'Z', -1; ...
   };

if(opt.use_pressure)
   deck.SCHEDULE.control.WCONINJE = ...
      { ...
      'I01', 'WAT', 'OPEN', 'BHP', Inf, Inf, 500, Inf, 0, 0, ...
      };

   deck.SCHEDULE.control.WCONPROD = ...
      { ...
      'P01', 'OPEN', 'BHP', Inf, Inf, Inf, Inf, Inf, 200, 0, 0, 0, ...
      };
else
   % use scaled spe10 rates
   resv = 5000*stb * (numel(layers) / 85);
   deck.SCHEDULE.control.WCONINJE = ...
      { ...
      'I01', 'WAT', 'OPEN', 'RESV', NaN, resv, Inf, Inf, 0, 0, ...
      };

   p_press = convertTo(4000*psia, barsa);
   deck.SCHEDULE.control.WCONPROD = ...
      { ...
      'P01', 'OPEN', 'BHP', Inf, Inf, Inf, Inf, Inf, p_press, 0, 0, 0, ...
      };
end

% define grid
if(false)
   grdecl=grdeclBox(cartdims,physdims,0);

else
   grdecl.cartDims = deck.RUNSPEC.DIMENS;
   dvec=deal(physdims./cartdims);
   myfields={'DXV','DYV','DZV'};
   for i=1:3
      grdecl.(myfields{i})=repmat(dvec(i),[grdecl.cartDims(i),1]);
   end
end

deck.GRID=grdecl;
deck.GRID.ACTNUM=int32(ones(nc,1));

switch myperm
   case 'uniform'
      deck.GRID.PORO=repmat(poro,[nc,1]);
      deck.GRID.PERMX=repmat(perm,[nc,1]);
      deck.GRID.PERMY=repmat(perm,[nc,1]);
      deck.GRID.PERMZ=repmat(perm,[nc,1]);
   case 'spe10_layer'
      assert(cartdims(3)==numel(layers));
      rock= SPE10_rock(layers);
      deck.GRID.PORO=rock.poro;
      deck.GRID.PERMX=rock.perm(:,1);
      deck.GRID.PERMY=rock.perm(:,2);
      deck.GRID.PERMZ=rock.perm(:,3);
   otherwise
      error('No such case')
end

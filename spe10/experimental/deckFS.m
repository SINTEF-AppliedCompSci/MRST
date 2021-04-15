function  deck=deckFS(cartdims,physdims,varargin)
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

opt = struct('Verbose', mrstVerbose,...
             'perm', 100*milli*darcy,...
             'poro', 0.1,...
             'permcase','uniform',...
             'radius',0.125 );
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
      error('Case still to be defined');
      disp('Try to use spe10 perm and poro');
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
deck.RUNSPEC.START=734139;
%define solution solution
%pres=200*barsa;
% define REGIONS
deck.REGIONS.SATNUM=ones(nc,1);
%define props
%s=linspace(0,1,2)';alpha=2;
%deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
chop = @(x) min(max(0,x),1);
s_wc = 0.0; %0.2;
s_or = 0.0; %0.2;
s = linspace(s_wc, 1 - s_or, 41)'; alpha = 2;
s_star = (s - s_wc)/(1 - s_wc - s_or);
swof = [s, chop(s_star.^alpha), chop((1-s_star).^alpha), s*0.0];
swof = [swof; [1.0 1.0 0.0 0.0]];
deck.PROPS.SWOF{1} = swof;

pres = convertTo(convertFrom(6000, psia), barsa);
deck.SOLUTION.PRESSURE=ones(nc,1)*pres;
deck.SOLUTION.SWAT=ones(nc,1)*s_wc;
deck.SOLUTION.SOIL=ones(nc,1)*(1-s_wc);
%deck.PROPS.PVTW=[0 1 0 0.500000000000000 0];
%deck.PROPS.PVDO{1}=[0 1 2;8000 1 2;16000 1 2];
%deck.PROPS.ROCK=[4000 3.000000000000000e-06 NaN NaN NaN NaN];
%deck.PROPS.DENSITY=[900 1000 0.044000000000000];
deck.PROPS.DENSITY = [53*pound/(ft^3) 64*pound/(ft^3) 1];
deck.PROPS.ROCK=[pres convertTo(convertFrom(1e-6, 1/psia), 1/barsa) NaN NaN NaN NaN];
deck.PROPS.PVTW=[pres 1.01 convertTo(convertFrom(3e-6, 1/psia), 1/barsa) 0.300000000000000 0];
deck.PROPS.PVDO{1} = [ [300, 800, 8000]'*psia/barsa, [1.05, 1.02, 1.01]', [2.85, 2.99, 3]' ];
% define summary
deck.SUMMARY=[];
% difine SC
%%
deck.SCHEDULE.control.WELSPECS=...
{...
'I01'    'W'    [  ceil(cartdims(1)/2)]    [ ceil(cartdims(2)/2) ]    [4000*ft()]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
'P01'    'W'    [  1]    [1]                                          [4000*ft()]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
'P02'    'W'    [  cartdims(1)]    [1]                                [4000*ft()]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
'P03'    'W'    [  cartdims(1)]    [cartdims(2)]                      [4000*ft()]    'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
'P04'    'W'    [  1]    [cartdims(2)]                                [4000*ft()]    'OIL'    [0]    'STD'    'SHUT'    'NO'   [0]    'SEG'    [0]...
};
radius=opt.radius;
deck.SCHEDULE.control.COMPDAT=...
{...
 'I01'     [  ceil(cartdims(1)/2)]    [ ceil(cartdims(2)/2) ]   [1]    [cartdims(3)]    'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
 'P01'    [  1]    [1]     [1]    [cartdims(3)]                                         'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
 'P02'    [  cartdims(1)]    [1]   [1]    [cartdims(3)]                                 'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
 'P03'    [  cartdims(1)]    [cartdims(2)]   [1]    [cartdims(3)]                       'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
 'P04'    [  1]    [cartdims(2)]    [1]    [cartdims(3)]                                'OPEN'    [0]    [0]    [radius]    [-1]    [0]    'Default'    'Z'    [-1];...
};
use_pressure=false;
if(use_pressure)
   deck.SCHEDULE.control.WCONINJE=...
      {...
      'I01'  'WAT'  'OPEN'  'BHP'  [Inf]  [Inf]  [500]  [Inf]  [0]  [0]...
      };
   deck.SCHEDULE.control.WCONPROD=...
      {...
      'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [200]  [0]  [0]  [0];...
      'P02'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [200]  [0]  [0]  [0];...
      'P03'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [200]  [0]  [0]  [0];...
      'P04'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [200]  [0]  [0]  [0]
      };
else
   % use scaled spe10 rates

   deck.SCHEDULE.control.WCONINJE=...
      {...
      'I01'  'WAT'  'OPEN'  'RESV'  [NaN]  [5000*stb()*numel(layers)/85]  [Inf]  [Inf]  [0]  [0]...
      };
   p_press=4000*psia()/barsa();
   deck.SCHEDULE.control.WCONPROD=...
      {...
      'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press]  [0]  [0]  [0];...
      'P02'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press]  [0]  [0]  [0];...
      'P03'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press]  [0]  [0]  [0];...
      'P04'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press]  [0]  [0]  [0]
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
      grdecl.(myfields{i})=ones(grdecl.cartDims(i),1)*dvec(i);
   end
end
deck.GRID=grdecl;
deck.GRID.ACTNUM=int32(ones(nc,1));
switch myperm
   case 'uniform'
      deck.GRID.PORO=ones(nc,1)*poro;
      deck.GRID.PERMX=ones(nc,1)*perm;
      deck.GRID.PERMY=ones(nc,1)*perm;
      deck.GRID.PERMZ=ones(nc,1)*perm;
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

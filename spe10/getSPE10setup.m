function [G, W, rock] = getSPE10setup(varargin)
%Initialise properties for Model 2 of tenth SPE Comparative Solution Project
%
% SYNOPSIS:
%   [G, W, rock] = getSPE10setup
%   [G, W, rock] = getSPE10setup(layers)
%   [G, W, rock] = getSPE10setup(layers, wloc)
%
% PARAMETERS:
%   layers - Which of the 85 model layers to include in a specific test.
%            OPTIONAL.  If unspecified or empty, all 85 layers (a total of
%            60-by-220-by-85 == 1,122,000 cells) are included.
%
%            Some possible values are
%               layers = ( 1 : 35).';  %  Tarbert formation
%               layers = (36 : 85).';  %  Upper Ness formation
%
%   wloc   - Location of the five wells. OPTIONAL. If unspecified, we use
%            the default location
%               wloc     = [  1,   60,     1,   60,   30;
%                             1,    1,   220,  220,  110];
%
% RETURNS:
%   G    - MRST grid structure as described in grid_structure.
%
%   W    - Well structure.  Injector at 500 Bar, producers at 200 Bar.
%          Inner product 'ip_tpf'.
%
%   rock - Rock structure having fields 'perm' and 'poro' pertaining to
%          the specified layers.  Permeability values are in strict SI
%          (metres squared).
%
% EXAMPLE:
%   [G, W, rock] = getSPE10setup(85)
%
% SEE ALSO:
%   `getSPE10rock`, `makeSPE10DataAvailable`.

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

   % Define layers
   [cartDims, physDims, layers] = define_layers(varargin{:});

   % Define grid
   try
       mrstModule add libgeometry
       G = mcomputeGeometry(cartGrid(cartDims, physDims));
   catch
       G = computeGeometry(cartGrid(cartDims, physDims));
   end

   % Get reservoir rock properties in given layers
   rock = getSPE10rock(layers);

   % Define wells
   % All wells completed in all layers of (reduced?) model
   W = define_wells(G, rock, varargin{:});
end

%--------------------------------------------------------------------------

function [cartDims, physDims, layers] = define_layers(varargin)
   layers = (1 : 85) .';

   if nargin > 0 && isnumeric(varargin{1}) && ~isempty(varargin{1})
      layers = varargin{1};
   end

   % Define grid dimensions and well locations
   cartDims = [ 60, 220, numel(layers)];
   physDims = cartDims .* [ 20, 10, 2 ]*ft;
end

%--------------------------------------------------------------------------

function W = define_wells(G, rock, varargin)
   [wtype, wtarget, wrad, wloc, wname, sgn] = well_setup(varargin{:});

   W = [];
   for w = 1 : numel(wtype)
      compi = [ 1, 0, 0 ];
      if sgn(w) < 0, compi(:) = 1; end

      W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                       'Type', wtype{w}, 'Val', wtarget(w), ...
                       'Radius', wrad(w), 'Name', wname{w}, ...
                       'Sign', sgn(w), 'InnerProduct', 'ip_tpf', ...
                       'Comp_i', compi);
   end
end

%--------------------------------------------------------------------------

function [wtype, wtarget, wrad, wloc, wname, sgn] = well_setup(varargin)
   wtype   = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
   wtarget = [200,   200,   200,   200,   500] .* barsa();
   wrad    = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;

   wloc    = [  1,   60,     1,   60,   30;
                1,    1,   220,  220,  110];

   wname   = {'P1', 'P2', 'P3', 'P4', 'I1'};
   sgn     = [ -1 ,  -1 ,  -1 ,  -1 ,   1 ];

   if nargin > 1 && isnumeric(varargin{2})
      wloc = varargin{2};
   end
end

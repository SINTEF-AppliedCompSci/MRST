function [G, W, rock] = SPE10_setup(varargin)
%Initialise properties for Model 2 of tenth SPE Comparative Solution Project
%
% SYNOPSIS:
%   [G, W, rock] = SPE10_setup
%   [G, W, rock] = SPE10_setup(layers)
%   [G, W, rock] = SPE10_setup(layers, wloc)
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
%   G    - SAMSIM grid structure as described in grid_structure.
%   W    - Well structure.  Injector at 500 Bar, producers at 200 Bar.
%   rock - Rock structure having fields 'perm' and 'poros' pertaining to
%          the specified layers.
%
% EXAMPLE:
%   [G, W, rock] = SPE10_setup(85)
%
% SEE ALSO:
%   SPE10_rock, computeGeometry.


%% Define layers
layers = (1 : 85) .';
if nargin > 0 && isnumeric(varargin{1}) && ~isempty(varargin{1}),
   layers = varargin{1};
end

%% Define grid dimensions and well locations
cartDims = [  60,  220,   numel(layers)];
physDims = [1200, 2200, 2*cartDims(end)] .* ft();   % ft -> m

wtype    = {'bhp', 'bhp', 'bhp', 'bhp', 'bhp'};
wtarget  = [200,   200,   200,   200,   500] .* barsa();
wrad     = [0.125, 0.125, 0.125, 0.125, 0.125] .* meter;
wloc     = [  1,   60,     1,   60,   30;
              1,    1,   220,  220,  110];
wname    = {'P1', 'P2', 'P3', 'P4', 'I1'};
sgn      = [ -1 ,  -1 ,  -1 ,  -1 ,   1 ];
if nargin > 1 && isnumeric(varargin{2}),
   wloc = varargin{2};
end

%% Get reservoir rock properties in given layers
rock      = SPE10_rock(layers);
rock.perm = convertFrom(rock.perm, milli*darcy);

%% Define grid
if ~readCache({cartDims}, 'verbose', false),
   G = cartGrid(cartDims, physDims);
   G = computeGeometry(G);

   writeCache({cartDims}, {'G'});
end

%% Define wells
% All wells completed in all layers of (reduced?) model

W = [];
for w = 1 : numel(wtype),
   W = verticalWell(W, G, rock, wloc(1,w), wloc(2,w), [], ...
                    'Type', wtype{w}, 'Val', wtarget(w), ...
                    'Radius', wrad(w), 'Name', wname{w}, ...
                    'Sign', sgn(w));
end

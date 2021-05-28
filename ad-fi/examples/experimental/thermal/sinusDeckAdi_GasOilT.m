function [deck,G] = sinusDeckAdi_GasOilT(cartDims, physDims, nsteps, dt, theta,...
                                 depth, phi, perm, rate, p_press)
%Make a GRDECL structure for simple sloping sinus-formed reservoir.
%
% SYNOPSIS:
%   [deck,G] = sinusDeckAdi(n, L, nt, dt, theta, depth, phi, K, r, p)
%
% PARAMETERS:
%   n     - Three-element vector, [nx, ny, nz], specifying the number of
%           cells in the 'x', 'y', and 'z' coordinate directions
%           respectively.
%   L     - Three element vector [Lx, Ly, H], specifying the physical
%           dimensions of the reservoir
%   nt    - Number of pressure steps (scalar)
%   dt    - The length of each time step (scalar)
%   theta - Angle of inclination of the reservoir (scalar)
%   depth - Reference depth for the reservoir (scalar)
%   phi   - Constant porosity (scalar)
%   K     - Constant isotropic permeability (scalar value)
%   r     - Injection rate for CO2
%   p     - Reference pressure (scalar)
%
% RETURNS:
%   deck  - A deck structure suitable for further processing by deck
%           functions like 'initEclipseGrid'
%   G     - Grid structure. Since this is constructed as part of making
%           the deck, we add the optional possibility of reusing it.

% SEE ALSO:
%   compareMethods

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

nc=prod(cartDims);

% define runspec
deck.RUNSPEC.cartDims = cartDims;
deck.RUNSPEC.DIMENS   = cartDims;
deck.RUNSPEC.OIL      = 1;
deck.RUNSPEC.GAS      = 1;
deck.RUNSPEC.METRIC   = 1;
deck.RUNSPEC.TABDIMS  = [1 1 20 50 20 50 1 20 20 1 10  1 -1  0 1];
deck.RUNSPEC.WELLDIMS = [5 10 2 1 5 10 5 4 3 0 1 1];
deck.RUNSPEC.AQUDIMS  = [0 0 0 0 10 10 0 0];
deck.RUNSPEC.START    = 734139;

% one sat num region
deck.REGIONS.SATNUM=ones(nc,1);

%define props
deck.PROPS.DENSITY = [1000 1000 600];
deck.PROPS.ROCK=[100 0 NaN NaN NaN NaN];
deck.PROPS.PVCDO=[100 1.0 0.0 0.40 0];
deck.PROPS.PVDG{1}=[100 1 0.1;1000 1 0.1];
s                  = linspace(0,1,2)';
drho               = (deck.PROPS.DENSITY(1)-deck.PROPS.DENSITY(3));
deck.PROPS.SGOF{1} = [s, s, 1-s, ((1-s)*drho*norm(gravity)*physDims(3))/barsa];

% define summary
deck.SUMMARY=[];

% define grid
deck.GRID = grdeclSloping(cartDims, physDims, ...
   'theta',theta,'amp',physDims(3)/5,'lambda',physDims(1)/2);
deck.GRID.ZCORN  = deck.GRID.ZCORN + min(deck.GRID.ZCORN) + depth;
%deck.GRID.ACTNUM = int32(ones(nc,1));
deck.GRID.PORO   = ones(nc,1)*phi;
deck.GRID.PERMX  = ones(nc,1)*perm;
deck.GRID.PERMY  = ones(nc,1)*perm;
deck.GRID.PERMZ  = ones(nc,1)*perm;

% define solution
deck.SOLUTION.SGAS = zeros(nc,1);
deck.SOLUTION.SOIL = ones(nc,1);

% compute hydrostatic pressure distribution
G = initEclipseGrid(convertDeckUnits(deck));
G = computeGeometry(G);
deck.SOLUTION.PRESSURE = convertTo(...
   p_press*barsa + norm(gravity)*deck.PROPS.DENSITY(2)*G.cells.centroids(:,3),...
   barsa);

% define well schedule
deck.SCHEDULE.control.WELSPECS=...
    {...
    'I01'  'W' [ceil(cartDims(1)/2)]    [ceil(cartDims(2)/2) ]   [1000]        'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    'P01'  'W' [cartDims(1)]            [cartDims(2)]            [1004.03647]  'OIL'    [0]    'STD'    'SHUT'    'NO'    [0]    'SEG'    [0];...
    }; %#ok
radius=0.01;
deck.SCHEDULE.control.COMPDAT=...
    {...
    'I01'  [ceil(cartDims(1)/4)]    [ ceil(cartDims(2)/2) ]   [1]    [cartDims(3)]    'OPEN'  [0]  [0] [radius] [-1] [0]   'Default'  'Z'    [-1];...
    'P01'  [cartDims(1)]            [cartDims(2)]             [1]    [cartDims(3)]    'OPEN'  [0]  [0] [radius] [-1] [0]   'Default'  'Z'    [-1];...
    };%#ok

% use scaled spe10 rates
deck.SCHEDULE.control.WCONINJE=...
    {...
    'I01'  'OIL'  'OPEN'  'BHP'  [rate]  [rate]  [300-89.4018]  [Inf]  [0]  [0] [0] [0] [0] [0] ...
    };%#ok
deck.SCHEDULE.control.WCONPROD=...
    {...
    'P01'  'OPEN'  'BHP'  [Inf]  [Inf]  [Inf]  [Inf]  [Inf]  [p_press-89.4018]  [0]  [0]  [0];...
    };%#ok
deck.SCHEDULE.control = [deck.SCHEDULE.control;deck.SCHEDULE.control];
%%deck.SCHEDULE.control(2).WCONINJE{3} = 'SHUT';
deck.SCHEDULE.step.val     = ones(nsteps,1)*dt/day;
%deck.SCHEDULE.step.val     = [deck.SCHEDULE.step.val,deck.SCHEDULE.step.val*40];
deck.SCHEDULE.step.control = [ones(nsteps,1)];%,2*ones(nsteps,1)];

% write deck to file
%writeDeck(deck, fullfile(datadir, mfilename));

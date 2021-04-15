function deck = initializeDeck()
%Initialise simple structure of ECLIPSE format type
%
% SYNOPSIS:
%   deck = initializeDeck()
%
% RETURNS
%   deck - Output structure containing the known, which is empty,
%          This structure contains the following fields, corresponding to
%          separate deck sections documented in the ECLIPSE or FrontSim
%          manuals:
%
%              RUNSPEC  - Basic meta data (dimensions &c).
%              GRID     - Grid specification (e.g., COORD/ZCORN and PERMX)
%              PROPS    - Fluid properties (e.g., PVT and relperm).
%              REGIONS  - Region partition.  May be empty.
%              SOLUTION - Initial conditions.
%              SCHEDULE - Time step and well controls (&c).
%
% SEE ALSO:
%   `readEclipseDeck`, `convertDeckUnits`, `initEclipseGrid`,
%   `initEclipseRock`, `initEclipseFluid`.

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


   deck.RUNSPEC  = initRunSpec;
   deck.GRID     = initGrid;
   deck.PROPS    = initProps;
   deck.REGIONS  = initRegions;
   deck.SOLUTION = initSolution;
   deck.SUMMARY  = initSummary;
   deck.SCHEDULE = initSchedule;

   sect  = reshape(fieldnames(deck), 1, []);
   empty = { {} };
   a     = [ sect ; repmat({ empty }, [ 1, numel(sect) ]) ];

   deck.UnhandledKeywords = struct(a{:});
end

%--------------------------------------------------------------------------

function runspec = initRunSpec()
   runspec.cartDims = repmat(-1, [1, 3]);
   runspec.START    = datenum('01 JAN 1983');  % E100 default value.
   runspec.DUALPORO = false;

   satopts = { 'DIRECT', 'HYSTER', 'IRREVERS' };
   runspec.SATOPTS = ...
      cell2struct(repmat({ false }, [numel(satopts), 1]), satopts, 1);
end

%--------------------------------------------------------------------------

function grid = initGrid()
   grid = struct();
end

%--------------------------------------------------------------------------

function props = initProps()
   props = struct();
end

%--------------------------------------------------------------------------

function regions = initRegions()
   regions = struct();
end

%--------------------------------------------------------------------------

function solution = initSolution()
   solution = struct();
end

%--------------------------------------------------------------------------

function summary = initSummary()
   summary = struct();
end

%--------------------------------------------------------------------------

function schedule = initSchedule()
   schedule.control = [];
   schedule.step    = struct('control', [], 'val', []);
end

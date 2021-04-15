function deck = readEclipseDeck(fn, varargin)
%Read Simplified ECLIPSE input deck.
%
% SYNOPSIS:
%   deck = readEclipseDeck(fn)
%
% PARAMETERS:
%   fn - String holding name of (readable) input deck specification.
%        Note: The specification is assumed to be a physical (seekable)
%        file on disk, not to be read in through (e.g.) a POSIX pipe.
%
% RETURNS:
%   deck - Output structure containing the known, though mostly
%          unprocessed, fields of the input deck specification.  This
%          structure contains the following fields, corresponding to
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
% NOTE:
%   This function does no error checking.  There is an implicit assumption
%   that the input deck is properly formed.
%
% SEE ALSO:
%   `convertDeckUnits`, `initEclipseGrid`, `initEclipseRock`, `initEclipseFluid`.

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


   [fid, msg] = fopen(fn, 'rt');
   if fid < 0, error([fn, ': ', msg]); end

   % Initialize input deck
   deck = initializeDeck;
   dirname = directory_name(fid);

   % Read input deck
   kw = getEclipseKeyword(fid);
   while ischar(kw)
      switch kw
         case 'RUNSPEC'
            deck = readRUNSPEC(fid, dirname, deck);

         case 'END'
            break;

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

          otherwise
            fclose(fid);
            error(msgid('Keyword:Unexpected'), ...
                  'Unexpected keyword ''%s'' in input deck.', kw);
      end

      kw = getEclipseKeyword(fid);
   end

   if isfield(deck, 'GRID') && isstruct(deck.GRID)
      if isfield(deck.GRID, 'ACTNUM')
         deck.GRID.ACTNUM = int32(reshape(deck.GRID.ACTNUM, [], 1));
      end

      % Convenience for other parts of MRST (e.g., processGRDECL).
      % Export DIMENS to GRID section.
      %
      deck.GRID.cartDims = deck.RUNSPEC.cartDims;
   end

   fclose(fid);

   if ~any(isfield(deck.RUNSPEC, {'METRIC', 'FIELD', 'LAB'}))
      % Default to METRIC conventions as per ECLIPSE specification.

      dispif(mrstVerbose, ...
            ['Unit system not specified in input deck. ', ...
             'Defaulting to ''METRIC'' conventions.\n']);

      deck.RUNSPEC.METRIC = true;
   end

   assert (sum(isfield(deck.RUNSPEC, {'METRIC', 'FIELD', 'LAB'})) == 1, ...
          ['Exactly one unit system must be specified in ', ...
           'input deck (defaulting to ''METRIC'' conventions).']);
end

%--------------------------------------------------------------------------

function dname = directory_name(fid)
   dname = fileparts(fopen(fid));
end

function [G, rock, fluid, ...
   state, wells, trans, dt, deck] = initEclipseModel(inputfile)
%Initialise basic MRST objects from ECLIPSE input file (*.DATA)
%
% SYNOPSIS:
%   [G, rock, fluid, ...
%    state, wells, trans, dt, deck] = initEclipseModel(inputfile)
%
% DESCRIPTION:
%   Function 'initEclipseModel' constructs the fundamental MRST simulation
%   objects (grid, rock properties, fluids &c) from a single ECLIPSE-style
%   input file (typically named '*.DATA').
%
% PARAMS:
%   inputfile - Name (string) of input file.  The input is assumed to be a
%               regular file on disk and not, say, a POSIX pipe.
%
% RETURNS:
%   G     - Grid structure as defined by function 'initEclipseGrid'.
%
%   rock  - Rock data structure as defined by function 'initEclipseRock'.
%
%   fluid - Fluid data structure as defined by function 'initEclipseFluid'.
%
%   state - Initial state object as defined by function 'initEclipseState'.
%
%   wells - Initial well object as defined by function 'processWells'. This
%           structure is derived from the initial well configuration.
%
%   trans - Background (absolute) one-sided ("half-face")
%           transmissibilities as defined by functions 'computeTrans' and
%           'computeTranMult'.
%
%   dt    - Report step sequence derived from 'TSTEP' and 'DATES'.
%
%   deck  - Raw input deck as defined by functions 'readEclipseDeck' and
%           'convertDeckUnits'.  Especially usefull for the dynamic
%           SCHEDULE information pertaining to switching well controls.

   require deckformat

   gravity reset
   gravity on

   deck  = convertDeckUnits(readEclipseDeck(inputfile));

   fprintf('----------------------------------------------------------\n');

   if isfield(deck.RUNSPEC, 'TITLE'),
      fprintf(' %s\n', deck.RUNSPEC.TITLE);
   else
      fprintf(' Untitled run\n');
   end

   fprintf('----------------------------------------------------------\n');

   G = computeGeometry(initEclipseGrid(deck));
   if ~ all(G.cells.volumes > 0),
      nc0 = G.cells.num;

      dispif(mrstVerbose, 'Removing cells with non-positive volume ... ');
      G = removeCells(G, ~ (G.cells.volumes > 0));
      dispif(mrstVerbose, '%d cells removed.\n', nc0 - G.cells.num);
   end

   rock  = compressRock(initEclipseRock(deck), G.cells.indexMap);
   fluid = initEclipseFluid(deck);

   state = initEclipseState(G, deck, fluid);
   if isfield(deck.SCHEDULE, 'control') && ...
         ~isempty(deck.SCHEDULE.control),

      wells = processWells(G, rock, deck.SCHEDULE.control(1), ...
                           'InnerProduct', 'ip_tpf');
      state.wellSol = initWellSol(wells, max(state.pressure)); % Galt

   else

      wells = [];

   end
   %state  = computeFacePressure(state, G, htrans, fluid);

   trans = computeTrans(G, rock);
   m     = computeTranMult(G, deck.GRID);
   if any(m),
      trans = trans.*m;
   end

   dt = deck.SCHEDULE.step.val;
end

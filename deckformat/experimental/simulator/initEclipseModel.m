function [G, rock, fluid, state, wells, ...
      htrans, trans, dt, deck] = initEclipseModel(inputfile,varargin)
%Initialise basic MRST objects from ECLIPSE input file (*.DATA)
%
% SYNOPSIS:
%   [G, rock, fluid, state, wells, ...
%    htrans, trans, dt, deck] = initEclipseModel(inputfile,'pn1','pv1',..)
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
% OPTIONAL PARAMETERS:
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  These
%            are passed on to the function 'processGRDECL', which is used
%            to construct corner-point grids. Of particular interest is
%            the parameter 'SplitDisconnected' which tells whether or not
%            to split disconnected grid components into separate
%            grids/reservoirs.
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
%   htrans -
%           Background (absolute) one-sided transmissibilities as defined
%           by functions 'computeTrans' and 'computeTranMult'.  Also
%           incorporates net-to-gross factors if available.
%
%   trans - Background (absolute) connection transmissibilities defined by
%           harmonic averaging of the one-sided transmissibilities (htrans)
%           while including fault-related transmissibility multipliers too.
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

   G = computeGeometry(initEclipseGrid(deck,varargin{:}));
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

   [htrans, trans] = transmissibility(G, rock, deck);

   dt = deck.SCHEDULE.step.val;
end

%--------------------------------------------------------------------------

function [htrans, trans] = transmissibility(G, rock, deck)
   htrans = computeTrans   (G, rock);
   htmult = computeTranMult(G, deck.GRID);

   if ~isempty(htmult),
      htrans = htrans .* htmult;
   end

   if isfield(rock, 'ntg') && (numel(rock.ntg) == G.cells.num),
      assert (all(isnumeric(rock.ntg)) && ...
              all(isfinite (rock.ntg)) && ...
              ~ any(rock.ntg < 0),        ...
             ['NTG must have one non-negative, finite value ', ...
              'for each active cell']);

      i         = false([max(G.cells.faces(:,2)), 1]);
      i(1 : 4)  = true;

      j         = i(G.cells.faces(:,2));
      cellno    = gridCellNo(G);
      htrans(j) = htrans(j) .* rock.ntg(cellno(j));
   end

   trans = 1 ./ accumarray(G.cells.faces(:,1), ...
                           1 ./ htrans,        ...
                           [G.faces.num, 1]);

   flt = processFaults(G, deck.GRID);
   if ~isempty(flt),
      trans = trans .* accumarray(vertcat(flt.faces), ...
                                  rldecode(vertcat(flt.mult), ...
                                           vertcat(flt.numf)), ...
                                  [G.faces.num, 1], @prod, 1);
   end
end

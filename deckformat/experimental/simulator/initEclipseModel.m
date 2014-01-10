function [G, rock, fluid, state, wells, trans, dt, deck] = initEclipseModel(fn)
   require deckformat

   gravity reset
   gravity on

   deck  = convertDeckUnits(readEclipseDeck(fn));

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

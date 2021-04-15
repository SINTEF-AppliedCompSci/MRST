function [G, rock_fracture, rock_matrix, fluid_matrix, fluid_fracture, state, wells, ...
      htransf, transf, htransm, transm, dt, deck, DP_info] = initEclipseModelDP(inputfile,varargin)

  
  
  
  %Initialise basic MRST objects from ECLIPSE input file (*.DATA)
%
% SYNOPSIS:
%   [G, rock_fracture, rock_matrix, fluid_matrix, fluid_fracture, state, wells, ...
%      htransf, transf, htransm, transm, dt, deck, DP_info] = initEclipseModelDP(inputfile,'pn1','pv1',..)
%
% DESCRIPTION:
%   Function 'initEclipseModelDP' constructs the fundamental MRST simulation
%   objects (grid, rock properties, fluids &c) from a single ECLIPSE-style
%   input file (typically named '*.DATA'). This reader has been modifed to
%   read in and setup dual porosity/permeability models.
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
%   rock_fracture  - 
%           Fracture Rock data structure as defined by function 'initEclipseRockDP'.
%
%   rock_matrix  - 
%            Matrix Rock data structure as defined by function 'initEclipseRockDP'.
%
%   fluid_matrix - 
%            Matrix Fluid data structure as defined by function 'initEclipseFluid'.
%
%   fluid_fracture - 
%           Fracture Fluid data structure as defined by function 'initEclipseFluid'.
%
%   wells - 
%           Initial well object as defined by function 'processWells'. This
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
%
%   DP_info - 
%           Dual porsity information sigma, matrix block sizes 
%
%   NOT Handled for DP models Yet       
%   state - 
%           Initial state object as defined by function 'initEclipseState'.


% SEE ALSO:
%   initEclipseGrid, initEclipseState, processWells.

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

   %require deckformat

   gravity reset
   gravity on

   deck  = convertDeckUnits(readEclipseDeck(inputfile));
   
   %ACTNUM conversion for DP Models
    deck.GRID.ACTNUM = reshape(deck.GRID.ACTNUM, [prod(deck.RUNSPEC.cartDims) 2]);
    FracACTNUM=deck.GRID.ACTNUM(:,2);
    deck.GRID.ACTNUM=deck.GRID.ACTNUM(:,1);


   fprintf('----------------------------------------------------------\n');

   if isfield(deck.RUNSPEC, 'TITLE'),
      fprintf(' %s\n', deck.RUNSPEC.TITLE);
   else
      fprintf(' Untitled run\n');
   end

   fprintf('----------------------------------------------------------\n');
   G = initEclipseGrid(deck,varargin{:});
   if numel(G) > 1
       warning(['Grid processing resulted in %d disconnected grids. '...
                'Picking the grid with the most cells.'], numel(G));
       G = G(1);
   end
   G = computeGeometry(G);
   if ~ all(G.cells.volumes > 0),
      nc0 = G.cells.num;

      dispif(mrstVerbose, 'Removing cells with non-positive volume ... ');
      G = removeCells(G, ~ (G.cells.volumes > 0));
      dispif(mrstVerbose, '%d cells removed.\n', nc0 - G.cells.num);
   end

   [rock_matrix,rock_fracture] = initEclipseDPRock(deck)
    rock_matrix = compressRock(rock_matrix, G.cells.indexMap)
    rock_fracture = compressRock(rock_fracture, G.cells.indexMap)

  
   
   [fluid_matrix, fluid_fracture] = initDeckADIFluidDP(deck, 'G', G)

   state=[];
   
    rock_matrix.satnum=deck.REGIONS.SATNUM(G.cells.indexMap,1);
    rock_fracture.satnum=deck.REGIONS.SATNUM(G.cells.indexMap,2);
   
   
   
   if isfield(deck.SCHEDULE, 'control') && ...
         ~isempty(deck.SCHEDULE.control),
                      
      wells = processWellsDP(G, rock_matrix, rock_fracture, deck.SCHEDULE.control(1), 'DualPoro', true);

   else

      wells = [];

   end
   

   %Create dual porosity info (not fully done yet)
        
        if isfield(deck.GRID, 'SIGMAV')
                    DP_info.sigma=deck.GRID.SIGMAV(G.cells.indexMap);
        end
        
        if isfield(deck.GRID, 'DZMTRXV')
                DP_info.DZmatrix=deck.GRID.DZMTRXV(G.cells.indexMap);
        end
        


   [htransf, transf] = transmissibility(G, rock_fracture, deck);
   [htransm, transm] = transmissibility(G, rock_matrix, deck);

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

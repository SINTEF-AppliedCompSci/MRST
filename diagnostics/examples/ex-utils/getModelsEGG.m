function [models, wells] = getModelsEGG(realizations)
% Helper function for the example script preprocessDiagnosticsEgg.m
%
% SYNOPSIS:
%   [models,wells] = getModelsEGG(realizations)
%
% DESCRIPTION:
%   This function returns models and wells for specified realisations of the Egg
%   ensemble in the correct format for input into DiagnosticsViewer(), the flow
%   diagnostics preprocessing GUI. 
%
% PARAMETERS:
%   realizations  - an array containing the desired realisations.
%
% RETURNS:
%     models - A structure of size numel(realizations) where models{i} is
%              an instance of GenericBlackOilModel containing grid, rock and
%              fluid properties for realizations(i).
%
%      wells - A structure of size numel(realizations) where wells{i} is
%              contains an mrst well structure for realizations(i).
%
% EXAMPLE:
%   [models, wells] = getModelsEGG([1:10]);
%
% SEE ALSO:
%   `preprocessDiagnosticsEgg`, `getModelsSAIGUP`
    
    [models, wells] = deal(cell(1, numel(realizations)));

    for k = 1:numel(realizations)
        
        deck = getDeckEGG('realization', realizations(k));
        
        G = initEclipseGrid(deck);
        G = computeGeometry(G);
        
        rock  = initEclipseRock(deck);
        rock  = compressRock(rock, G.cells.indexMap);
        
        fluid = initDeckADIFluid(deck);
        
        models{k} = selectModelFromDeck(G, rock, fluid, deck);
        
        schedule = convertDeckScheduleToMRST(models{k}, deck);
        wells{k}  = schedule.control(1).W;
    end
end


%{
  Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
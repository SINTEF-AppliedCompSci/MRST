classdef MpfaBlackOilModel < GenericBlackOilModel

    properties
        eta
        bcetazero
    end

    methods
        
        function model = MpfaBlackOilModel(G, rock, fluid, varargin)
            
            model = model@GenericBlackOilModel(G, rock, fluid, varargin{:});
            model = merge_options(model, varargin{:});
            
            % Process the grid for mechanical computation
            if ~ismember('createAugmentedGrid', model.G.type)
                model.G = createAugmentedGrid(model.G);
            end
        
            model.eta = 0;
            model.bcetazero = false;
            
            % Add mechanical operators  
            model.operators = setupMpfaAdOperators(model);

        end


        function model = setupStateFunctionGroupings(model, varargin) 
            
            model = setupStateFunctionGroupings@GenericBlackOilModel(model, varargin{:});
            fd = model.FlowDiscretization;
            fd = fd.setStateFunction('PermeabilityPotentialGradient', MpfaKgrad(model));
            fd = fd.setStateFunction('PhaseUpwindFlag', MpfaPhaseUpwindFlag(model));
            model.FlowDiscretization = fd;
            
        end


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

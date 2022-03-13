classdef FlowDiscretizationDG < FlowDiscretization
   
    properties
        ComponentTotalVelocity
        ComponentPhaseVelocity
        TotalVelocity
    end
    
    methods
        function fd = FlowDiscretizationDG(model)
            % Initialize
            fd = fd@FlowDiscretization(model);
            % Copy over properties from model.FlowDiscretization
            if ~isempty(model.FlowDiscretization)
                names = model.FlowDiscretization.getNamesOfStateFunctions();
                for i = 1:numel(names)
                    p = model.FlowDiscretization.getStateFunction(names{i});
                    fd = fd.setStateFunction(names{i}, p);
                end
            end
            % Set dg-specific state functions
            fd = fd.setStateFunction('TotalFlux', FixedTotalFluxDG(model));
            fd = fd.setStateFunction('GravityPotentialDifference', GravityPotentialDifferenceDG(model));
            fd = fd.setStateFunction('PhaseUpwindFlag', PhasePotentialUpwindFlagDG(model));
            fd = fd.setStateFunction('ComponentTotalVelocity', ComponentTotalVelocityDG(model));
            fd = fd.setStateFunction('ComponentPhaseVelocity', ComponentPhaseVelocityFractionalFlowDG(model));
            fd = fd.setStateFunction('TotalVelocity', FixedTotalVelocityDG(model));
            fd = fd.setFlowStateBuilder(FlowStateBuilderDG);
        end
        
        function [acc, v, vc, names, types] = componentConservationEquations(fd, model, state, state0, dt)
            % Compute discretized conservation equations in the interior of the domain.
            ncomp = model.getNumberOfComponents;
            [acc, types] = deal(cell(1, ncomp));
            names = model.getComponentNames();
            [types{:}] = deal('cell');
            % Get mass terms
            model = fd.expandRegions(model, 'cells');
            mass  = model.getProps(state.cellStateDG, 'ComponentTotalMass');
            mass0 = model.getProps(state0.cellStateDG, 'ComponentTotalMass');
            % Get boundary flux term
            flowStateFace = fd.buildFlowState(model, state, state0, dt, 'face');
            model     = fd.expandRegions(model, 'faces');
            v         = model.getProps(flowStateFace, 'ComponentTotalFlux');
            % Get cell flux term
            flowStateCell = fd.buildFlowState(model, state, state0, dt, 'cell');
            model = fd.expandRegions(model, 'cells');
            vc = model.getProps(flowStateCell, 'ComponentTotalVelocity');
            for c = 1:ncomp
                acc{c} = (mass{c} - mass0{c})./dt;
            end
        end
        
        function model = expandRegions(fd, model, type, elements)
            
            if nargin < 4
                elements = [];
            end
            switch type
                case 'cells'
                    cells = elements;
                    if isempty(cells)
                        [~, ~, cells] = model.operators.discretization.getCubature((1:model.G.cells.num)', 'cell');
                    end
                case 'faces'
                    faces = elements;
                    if isempty(faces)
                        [~, ~, ~, faces] = model.operators.discretization.getCubature(find(model.operators.internalConn), 'face');
                    end
                    cells = [model.G.faces.neighbors(faces,1); model.G.faces.neighbors(faces,2)];
            end
            
            pvtprops = model.PVTPropertyFunctions;
            pvtprops = expandPropsRegions(pvtprops, cells);
            model.PVTPropertyFunctions = pvtprops;
            
            fpprops = model.FlowPropertyFunctions;
            fpprops = expandPropsRegions(fpprops, cells);
            model.FlowPropertyFunctions = fpprops;
            
            fdprops = model.FlowDiscretization;
            fdprops = expandPropsRegions(fdprops, cells);
            model.FlowDiscretization = fdprops;
            
        end
        
        function flowState = buildFlowState(fd, model, state, state0, dt, type)
            if nargin < 6
                type = 'face';
            end
            flowState = fd.FlowStateBuilder.build(fd, model, state, state0, dt, type);
        end
        
    end
    
end

function props = expandPropsRegions(props, cells)
    names = props.getNamesOfStateFunctions();
    for j = 1:numel(names)
        if isprop(props, names{j})
            p = props.(names{j});
            p = p.subset(cells);
            props.(names{j}) = p;
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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

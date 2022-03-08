classdef AdvectiveHeatFlux < StateFunction & UpwindProperty
%State function for computing advective heat flux
    
    properties (Access = protected)
        upwind_name; % Name of state function where upwind flag comes from
    end
    
    methods
        %-----------------------------------------------------------------%
        function hf = AdvectiveHeatFlux(model, varargin)
            if nargin < 2
                upstr = UpwindFunctionWrapperDiscretization(model);
            end
            if nargin < 3
                upwind_name = 'PhaseUpwindFlag';
            end
            hf@StateFunction(model);
            hf@UpwindProperty(upstr)
            hf.upwind_name = upwind_name;
            % Dependenceies are phase flux, phase enthalpy and density
            hf = hf.dependsOn({upwind_name, 'ComponentPhaseFlux'});
            hf = hf.dependsOn({'PhaseEnthalpy'}, 'PVTPropertyFunctions');
            hf.label = 'H_{a,\alpha}';
        end
        
        %-----------------------------------------------------------------%
        function H = evaluateOnDomain(prop, model, state)
            % Get dependencies
            [v, flag] = prop.getEvaluatedDependencies(state, 'ComponentPhaseFlux', ...
                                                              prop.upwind_name   );
            h     = model.getProps(state, 'PhaseEnthalpy');
            nph   = model.getNumberOfPhases();
            ncomp = model.getNumberOfComponents();
            H     = cell(1, nph);
            for i = 1:nph
                vph = 0;
                for j = 1:ncomp
                    if isempty(v{j,i}), continue; end
                    vph = vph + v{j,i};
                end
                H{i} = prop.faceUpstream(model, state, flag{i}, h{i}).*vph;
            end
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
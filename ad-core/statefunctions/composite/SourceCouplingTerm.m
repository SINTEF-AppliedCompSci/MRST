classdef SourceCouplingTerm < CouplingTerm
   
    properties
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function prop = SourceCouplingTerm(model, modelA, modelB, varargin)
            
            % Optional input arguments
            opt = struct( ...
                'cellsA', [], ...
                'cellsB', []  ...
            );
            [opt, extra] = merge_options(opt, varargin{:});
            % Parent class constructor
            prop@CouplingTerm(model, extra{:});
            
            if isempty(prop.couplings)
                % Check that models share the same component names
                mA     = model.submodels.(modelA).validateModel();
                namesA = mA.getComponentNames();
                mB     = model.submodels.(modelB).validateModel();
                namesB = mB.getComponentNames();
                assert(all(ismember(namesA, namesB)), ...
                    ['Moldes do not share the same component names, '   , ...
                     'so I can''t figure out how to couple the models. ', ...
                     'Please specify coupling explicitly'               ]);
                % Specify coupling
                coupling = { ...
                    struct( ...
                        'model'    , modelA    , ... % Model name
                        'equations', {namesA}  , ... % Equations
                        'subset'   , opt.cellsA, ... % From cellsA
                        'sign'     , -1          ... % Sink
                    ), ...
                    struct( ...
                        'model'    , modelB    , ... % Model name
                        'equations', {namesB}  , ... % Equations
                        'subset'   , opt.cellsB, ... % To cellsB
                        'sign'     , 1           ... % Source
                    ), ...
                };
                % Set to state function
                prop.couplings = coupling;
                
            end
            % Set propeties
            prop.label            = 'q';
            prop.submodels.modelA = modelA;
            prop.submodels.modelB = modelB;
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function Q = evaluateOnDomain(prop, model, state)
            
            % Get model name and cell subset
            mname = prop.submodels.modelA;
            cells = prop.couplings{1}.subset;
            % Everything entering cellsA in modelA is extraced and injected
            % into cellsB in modelB.
            
            % Use simple formulation where we extract everything that
            % enters the cells. Not entirely physical, but at least
            % consistent with the physics of the model.
            nph = model.submodels.(mname).getNumberOfPhases();
            q   = model.submodels.(mname).getProps( ... 
                    state.(mname), 'ComponentTotalFlux');
            Q = q;
            for ph = 1:nph
                qph = model.submodels.(mname).operators.Div(q{ph});
                Q{ph} = qph(cells);
            end
            
        end
        %-----------------------------------------------------------------%
            
    end
    
end

%{
Copyright 2009-2024 SINTEF Digital, Mathematics & Cybernetics.

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

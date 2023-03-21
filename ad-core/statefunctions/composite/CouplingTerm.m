classdef CouplingTerm < StateFunction
% Template state function for implementing coupling terms between models in
% CompositeModels
    
    properties
        
        couplings             % Cell array of coupling structs
        submodels = struct(); % Struct of submodel names
        
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function prop = CouplingTerm(model, varargin)
        % Constructor
        
            prop = prop@StateFunction(model);
            prop = merge_options(prop, varargin{:});
            
        end
        %-----------------------------------------------------------------%
        
        %-----------------------------------------------------------------%
        function eqs = insertCoupling(prop, eqs, modelNames, names, q)
        % Insert coupling terms based on couplings structs
        
            for c = 1:numel(prop.couplings)
                % Get coupling
                cpl = prop.couplings{c};
                % Find model indices
                mix = find(strcmpi(modelNames, cpl.model));
                for m = mix
                    % Find equation indices
                    eqix = find(ismember(names{m}, cpl.equations));
                    for eq = eqix
                        % If we have a cell array of sources/sinks, we
                        % interpret them to be per equation
                        if iscell(q), qeq = q{eq}; else, qeq = q; end
                        % Insert in model equation, with correct sign
                        eqs{m}{eq}(cpl.subset) ...
                            = eqs{m}{eq}(cpl.subset) + qeq.*cpl.sign;
                    end
                end

            end
            
        end
        %-----------------------------------------------------------------%
        
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
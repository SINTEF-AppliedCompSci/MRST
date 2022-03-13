classdef RecoveryFactorQoI < BaseQoI
    
    properties
        phase           = 'oil';
        diagnosticsType = 'tof';
        diagnosticsArgs = {};
    end
    
    methods
        
        %-----------------------------------------------------------------%
        function qoi = RecoveryFactorQoI(varargin)
            qoi = qoi@BaseQoI('names', {'rf'});
            qoi = merge_options(qoi, varargin{:});
        end
        
        %-----------------------------------------------------------------%
        function qoi = validateQoI(qoi, problem, varargin)
            % Let base class do its thing
            qoi   = validateQoI@BaseQoI(qoi, problem, varargin{:});
            % Check that phase exists
            model = problem.SimulatorSetup.model;
            [~, phases] = model.getPhaseNames();
            assert(any(strcmpi(qoi.phase, phases)), ...
                   'Did not find %s in model phases', qoi.phase);
        end
        
        %-----------------------------------------------------------------%
        function u = computeQoI(qoi, problem)
            % Get state output handler
            [~, states] = getPackedSimulatorOutput(problem, 'readFromDisk', false);
            % Get model
            model = problem.SimulatorSetup.model.validateModel();
            % Get first state
            state = states{1};
            if isfield(state, 'flowDiagnostics')
                rf = qoi.computeRecoveryFactorFD(problem, state);
            else
                rf = qoi.computeRecoveryFactor(model, states);
            end
            u = struct('rf', rf);
        end
        
        %-----------------------------------------------------------------%
        function u = computeRecoveryFactor(qoi, model, states)
            % Compute recovery factor from full simulation
            % Compute initial volume
            [si, pvi] = model.getProps(states{1}, qoi.phase, 'PoreVolume');
            vi = sum(si.*pvi);
            % Compute final volume
            [sf, pvf] = model.getProps(states{numelData(states)}, qoi.phase, 'PoreVolume');
            vf = sum(sf.*pvf);
            % Compute recovery factor
            u = (vi-vf)./vi;
        end
        
        %-----------------------------------------------------------------%
        function u = computeRecoveryFactorFD(qoi, problem, state)
            % Compute recovery factor from representative flux field using
            % flow diagnostics
            tD = state.flowDiagnostics.tD;
            Ev = state.flowDiagnostics.Ev;
            pv = state.flowDiagnostics.pv;
            pvi = qoi.computePVI(problem, state, pv);
            % Find sweep efficiency at pvi
            [~, ix] = min(abs(tD - pvi));
            u = Ev(ix);
        end
        
        %-----------------------------------------------------------------%
        function model = adjustPoreVolumes(qoi, model)
            % TODO: implement pore volume adjustment based on BL speed
        end
        
        %-----------------------------------------------------------------%
        function pvi = computePVI(qoi, problem, state, pv)
            % Compute pvi for the schedule based on representative wellSol
            % Get well solution
            wellSol = state.wellSol;
            % Find injectors
            injector = [wellSol.sign] > 0;
            % Get injection fluxes at reservoir conditions
            flux = sum(sum(vertcat(wellSol(injector).flux), 1),2);
            % Total injection time
            time = sum(problem.SimulatorSetup.schedule.step.val);
            % Total injected volume
            vol = flux*time;
            % Translate to pore volumes
            pvi = vol./sum(pv);
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

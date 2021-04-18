classdef RecoveryFactorQoI < BaseQoI
    
    properties
        phase           = 'oil';
        diagnosticsType = 'tof';
        diagnosticsArgs = {};
    end
    
    methods
        
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
            statei = states{1};
            if isfield(statei, 'singlePressureSolve') && statei.singlePressureSolve
                u = qoi.computeRecoveryFactorFD(problem, model, states);
            else 
                u = qoi.computeRecoveryFactor(model, states);
            end
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
            u = {(vi-vf)./vi};
        end
        
        %-----------------------------------------------------------------%
        function u = computeRecoveryFactorFD(qoi, problem, model, states)
            % Compute recovery factor from representative flux field using
            % flow diagnostics
            require diagnostics
            % Get state
            state = states{1};
            % Get wells
            W = problem.SimulatorSetup.schedule.control(1).W;
            % Adjust pore volumes to account for multiphase flow (does
            % nothing at the moment)
            model = qoi.adjustPoreVolumes(model);
            % Get pore volume
            pv = model.getProp(state, 'PoreVolume');
            % Compute time of flight
            maxTOF = sum(problem.SimulatorSetup.schedule.step.val)*10;
            D      = computeTOFandTracer(state, model.G, model.rock, ...
                     'wells', W, 'maxTOF', maxTOF, qoi.diagnosticsArgs{:});
            % Use TOF or residence time distribution (more accurate)
            switch qoi.diagnosticsType
                case 'tof'
                    % Compute flow capacity and storage capacity from tof
                    [F, Phi] = computeFandPhi(pv, D.tof);
                case 'rtd'
                    WP  = computeWellPairs(state, model.G, model.rock, W, D);
                    rtd = computeRTD(state, model.G, pv, D, WP, W, ...
                                                     'showWaitbar', false);
                    [F, Phi] = computeFandPhiFromDist(rtd, 'sum', true);
            end
            % Compute sweep efficiency
            [Ev, tD] = computeSweep(F, Phi);
            % Estimate total injected pore volumes
            pvi = qoi.computePVI(problem, pv);
            % Find sweep efficiency at pvi
            [~, ix] = min(abs(tD - pvi));
            u = {Ev(ix)};
        end
        
        %-----------------------------------------------------------------%
        function model = adjustPoreVolumes(qoi, model)
            % TODO: implement pore volume adjustment based on BL speed
        end
        
        %-----------------------------------------------------------------%
        function pvi = computePVI(qoi, problem, pv)
            % Compute pvi for the schedule based on representative wellSol
            % Get well solution
            wellSol = problem.OutputHandlers.wellSols{1};
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

classdef State0Samples < BaseSamples
    % Class for holding stochastic samples that represent uncertain
    % initial conditions, such as initial saturation, pressure, and fluxes 
    % for an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize stochastic initial state properties. These  
    %   properties can either be pre-computed or generaxted on the fly.
    % 
    %   Supported data fields are:
    %    - initSw: initial water saturation per cell
    %    - initSo: initial oil saturation per cell 
    %    - initSg: initial gas saturation per cell
    %    - pressure: pressure per cell
    %    - flux: fluxes per face
    % 
    % SYNOPSIS
    %   samples = OperatorSamples('data', data);
    %   samples = OperatorSamples('generatorFn', generatorFn);
    %
    % OPTIONAL PARAMETERS
    %   'data' - precomputed sample data. Supported formats:
    %               * Cell array of data samples
    %               * Instance of ResultHandler class with
    %                 information about storage location and names. See
    %                 ResultHandler class for details
    %
    %   'generatorFn' - Function for generating a stochastic sample
    % 
    % NOTE:
    %   Either 'data' or 'generatorFn' must be provided.
    %   Each sample should consist of a struct with the fields also found
    %   in SimulatorSetup.model.operators
    %
    % SEE ALSO:
    %   `RockSamples`, `DeckSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the initial states to a
            % problem.
            %
            % SYNOPSIS:
            %   problem = sample.setSample(sampleData, problem)
            %
            % PARAMETERS:
            %   sampleData - The data for the specific sample realization.
            %
            %   problem - The problem which the sampleData will be applied
            %             to.
            % RETURNS:
            %   problem - A problem representing a single ensemble member
            
            % For now, only transmissibility and porevolume is supported
                     
            
            initSwIsField = isfield(sampleData, 'initSw');
            initSoIsField = isfield(sampleData, 'initSo');
            initSgIsField = isfield(sampleData, 'initSg');
            pressureIsField = isfield(sampleData, 'pressure');
            fluxIsField = isfield(sampleData, 'flux');
            
            assert(initSwIsField + initSoIsField + initSgIsField + pressureIsField + ...
                fluxIsField == numel(fields(sampleData)), ...
                'Invalid fields in sample data for State0Samples');
            
            % Initial saturation
            if (initSwIsField || initSoIsField || initSgIsField)
                numPhases = size(problem.SimulatorSetup.state0.s, 2);       
                sampleData = samples.fillInitialSaturation(sampleData, numPhases, ...
                                                           initSwIsField, initSoIsField, initSgIsField);
                if (numPhases == 3)
                    problem.SimulatorSetup.state0.s(:, :) = ...
                        [sampleData.initSw(:) sampleData.initSo(:) sampleData.initSg(:)];
                elseif (numPhases == 2)
                    problem.SimulatorSetup.state0.s(:, :) = ...
                         [sampleData.initSw(:) sampleData.initSo(:)];
                end
                
            end % if initial saturation
            
            % Pressure
            if pressureIsField
                problem.SimulatorSetup.state0.pressure(:) = sampleData.pressure;
            end
            
            % flux
            if fluxIsField
                problem.SimulatorSetup.state0.flux(:) = sampleData.flux;
            end
            
        end % function setSample
        
        
    end
    
    methods (Access = protected)
        function sizes = getSizes(samples)
            sizeInitSw = 0;
            sizeInitSo = 0;
            sizeInitSg = 0;
            sizePressure = 0;
            sizeFlux = 0;
            
            if isfield(samples.data{1}, 'initSw')
                sizeInitSw = numel(samples.data{1}.initSw);
            end
            if isfield(samples.data{1}, 'initSo')
                sizeInitSo = numel(samples.data{1}.initSo);
            end
            if isfield(samples.data{1}, 'initSg')
                sizeInitSg = numel(samples.data{1}.initSg);
            end
            if isfield(samples.data{1}, 'pressure')
                sizePressure = numel(samples.data{1}.pressure);
            end
            if isfield(samples.data{1}, 'flux')
                sizeFlux = numel(samples.data{1}.flux);
            end
            
            sizes = [sizeInitSw, sizeInitSo, sizeInitSg, sizePressure, sizeFlux];
        end
        
        function sampleData = fillInitialSaturation(samples, sampleData, numPhases, ...
                                                    initSwIsField, initSoIsField, initSgIsField)

            if (numPhases == 3)
                assert(initSwIsField + initSoIsField + initSgIsField > 1, ...
                    'Three phase initial conditions require that two of the phases are given in sampleData. Currently, only one phase is available.');
                if (initSwIsField && initSoIsField && initSgIsField)
                    assert(all(sampleData.initSw + sampleData.initSo + sampleData.initSg == 1), ...
                        'Three phase initial saturation does not sum to 1');
                    % Nothing to do, sampleData is fine.
                    
                else %(initSwIsField + initSoIsField + initSgIsField == 2)
                    if (initSwIsField && initSoIsField)
                        sampleData.initSg = 1 - sampleData.initSw - sampleData.initSo;
                    elseif (initSwIsField && initSgIsField)
                        sampleData.initSo = 1 - sampleData.initSw - sampleData.initSg;
                    else %(initSoIsField && initSgIsField)
                        sampleData.initSw = 1 - sampleData.initSg - sampleData.initSo;
                    end
                end

            elseif (numPhases == 2)
                assert(~initSgIsField, 'Initial gas can not be used for two-phase problems');

                if (initSwIsField && initSoIsField)
                    assert(all(sampleData.initSw + sampleData.initSo == 1), ...
                        'Two phase initial saturation does not sum to 1 for two phase problem');
                     % Nothing to do, sampleData is fine.
               else % initSw or initSo
                    if initSwIsField
                        sampleData.initSo = 1 - sampleData.initSw;
                    else % initSoIsField
                        sampleData.initSw = 1 - sampleData.initSo;
                    end
                end

            end % if numphases
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

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
                
                if (numPhases == 3)
                    assert(initSwIsField + initSoIsField + initSgIsField > 1, ...
                        'Three phase initial conditions require that two of the phases are given in sampleData. Currently, only one phase is available.');
                    if (initSwIsField && initSoIsField && initSgIsField)
                        assert(all(sampleData.initSw + sampleData.initSo + sampleData.initSg == 1), ...
                            'Three phase initial saturation does not sum to 1');
                        problem.SimulatorSetup.state0.s(:, :) = ...
                            [sampleData.initSw(:) sampleData.initSo(:) sampleData.initSg(:)];
                        
                    else %(initSwIsField + initSoIsField + initSgIsField == 2)
                        if (initSwIsField && initSoIsField)
                            initSg = 1 - sampleData.initSw - sampleData.initSo;
                            problem.SimulatorSetup.state0.s(:, :) = ...
                                [sampleData.initSw(:) sampleData.initSo(:) initSg(:)];
                        elseif (initSwIsField && initSgIsField)
                            initSo = 1 - sampleData.initSw - sampleData.initSg;
                            problem.SimulatorSetup.state0.s(:, :) = ...
                                [sampleData.initSw(:) initSo(:) sampleData.initSg(:)];
                        else %(initSoIsField && initSgIsField)
                            initSw = 1 - sampleData.initSg - sampleData.initSo;
                            problem.SimulatorSetup.state0.s(:, :) = ...
                                [initSw(:) sampleData.initSo(:) sampleData.initSg(:)];
                        end
                    end
                    
                elseif (numPhases == 2)
                    assert(~initSgIsField, 'Initial gas can not be used for two-phase problems');
                    
                    if (initSwIsField && initSoIsField)
                        assert(all(sampleData.initSw + sampleData.initSo == 1), ...
                            'Two phase initial saturation does not sum to 1 for two phase problem');
                        problem.SimulatorSetup.state0.s(:, :) = ...
                            [sampleData.initSw(:) sampleData.initSo(:)];
                    else % initSw or initSo
                        if initSwIsField
                            initSo = 1 - sampleData.initSw;
                            problem.SimulatorSetup.state0.s(:, :) = ...
                                [sampleData.initSw(:) initSo(:)];
                        else % initSoIsField
                            initSw = 1 - sampleData.initSo;
                            problem.SimulatorSetup.state0.s(:, :) = ...
                                 [initSw(:) sampleData.initSo(:)];
                        end
                    end
                    
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
    
end

%{
#COPYRIGHT#
%}
function samples = createEggSamplesFromFlowDiagnostics( ...
    baseNetworkModel, eggRealizations, initializationType, ...
    varargin)
%Undocumented Utility Function

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

    opt = struct( ...
        'regenerateInitialEnsemble', false, ...
        'clearAllPackedSimulatorOutputs', false, ...
        'WItype', 'mean', ...
        'transMin', eps, ...
        'transMaxFactor', 1.5, ...,
        'TScale', 1e-9, ...
        'pvMin', 0.1, ...
        'pvMaxFactor', 1.5, ...
        'pvScale', 1e4, ...
        'maxWIFactor', 7, ...
        'minWIFactor', 0.01, ...
        'WIScale', 1e-11, ...
        'postprocStateNumber', 20, ...
        'fullEnsembleDirectory', 'eggModels' ...
        );         
    opt = merge_options(opt, varargin{:});
    
    initParamFolder = fullfile(mrstOutputDirectory(), ...
                         'eggEnsembleNetworkESMDACalibration_examples');
    initParamFile = fullfile(initParamFolder, ...
                         [initializationType, '_parameters.mat']);

    ensembleSize = numel(eggRealizations);
                     
    injectors = find([baseNetworkModel.schedule.control.W.sign] ==   1);
    producers = find([baseNetworkModel.schedule.control.W.sign] ==  -1);
    numInjectors = numel(injectors);
    numProducers = numel(producers);
    numConnections = numInjectors*numProducers;
    numWells = numel(baseNetworkModel.schedule.control.W);
    
    overrideExistingParamFile = true;
    if exist(initParamFile, 'file') && ~opt.regenerateInitialEnsemble
        % Read pre-computed initial ensemble
        preComputedEnsemble = load(initParamFile);
        transmissibilities = preComputedEnsemble.transmissibilities;
        poreVolumes = preComputedEnsemble.poreVolumes;
        wellProductionIndicesSum = preComputedEnsemble.wellProductionIndicesSum;
        wellProductionIndicesMean = preComputedEnsemble.wellProductionIndicesMean;

        % Check if we got sufficient data or if we have to overwrite the
        % existing data
        overrideExistingParamFile = (size(poreVolumes, 1) < ensembleSize);
        
    end
    
    if overrideExistingParamFile
        % Use flow diagnostic analysis to obtain initial parameters
        fprintf('Running flow diagnostics analysis. This may take some time.\n');

        % Start with values that correspond to "closed" connections only.
        transmissibilities = ones(numel(eggRealizations), numConnections).*opt.transMin;
        poreVolumes = ones(numel(eggRealizations), numConnections).*opt.pvMin;
        % The rows are organized as follows for the different connections:
        % 1 (inj) -  9 (prod)
        % 1 (inj) - 10 (prod)
        % 1 (inj) - 11 (prod)
        % 1 (inj) - 12 (prod)
        % 2 (inj) - 9 (prod)
        % 2 (inj) - 10 (prod)
        % ...
        % 8 (inj) - 11 (prod)
        % 8 (inj) - 12 (prod)

        wellProductionIndicesSum = zeros(numel(eggRealizations), numWells);
        wellProductionIndicesMean = zeros(numel(eggRealizations), numWells);

        for eggRealization = eggRealizations

            % Create the full Egg model realization
            ensembleCase = TestCase('egg_wo', 'realization', eggRealization);

            % We consider only the timesteps up until the step specified
            % in the variable 'postprocStateNumber'
            ensembleCase.schedule = simpleSchedule(ensembleCase.schedule.step.val(1:opt.postprocStateNumber), ...
                                                       'W', ensembleCase.schedule.control.W);
           % Pack the problem 
           ensembleProblem = ensembleCase.getPackedSimulationProblem('Directory', ...
                fullfile(opt.fullEnsembleDirectory, ['full_egg_', num2str(eggRealization)]));

            % Configure wells to only have one perforation and store sum of the
            % well production indices to use as prior.
            WtmpNetwork = ensembleCase.schedule.control.W;
            for i = 1:numel(WtmpNetwork)
                WtmpNetwork(i).cells = WtmpNetwork(i).cells(7);
                wellProductionIndicesSum(eggRealization, i) = sum(WtmpNetwork(i).WI);
                wellProductionIndicesMean(eggRealization, i) = mean(WtmpNetwork(i).WI);
            end

            % Network from flow diagnostic analysis
            switch initializationType
                case 'fd_preprocessor'
                    % Flow diagnostic analysis directly on the geological model
                    tmpNetwork = Network(WtmpNetwork, ensembleProblem.SimulatorSetup.model.G, ...
                                         'type', initializationType,         ...
                                         'problem', ensembleProblem,        ...
                                         'flow_filter',1*stb/day);
                case 'fd_postprocessor'
                    % Run the simulation 
                    if opt.clearAllPackedSimulatorOutputs
                        clearPackedSimulatorOutput(ensembleProblem, 'prompt',  false);
                    end
                    simulatePackedProblem(ensembleProblem);

                    % Then do postprocessing flow diagnostic analysis
                    tmpNetwork = Network(WtmpNetwork, ensembleProblem.SimulatorSetup.model.G, ...
                                         'type', initializationType,         ...
                                         'problem', ensembleProblem,        ...
                                         'state_number', opt.postprocStateNumber,     ...
                                         'flow_filter', 1*stb/day);
                otherwise
                    error('\nNetwork of type %s is not implemented\n', initializationType);    
            end

            % Map the values from the smaller network to the full network
            for connectionID = 1:height(tmpNetwork.network.Edges)
                connection = tmpNetwork.network.Edges(connectionID,:);
                inj  = connection.EndNodes(1);
                prod = connection.EndNodes(2);

                rowInFullNetwork = numProducers*(inj-1) + prod-numInjectors;

                transmissibilities(eggRealization, rowInFullNetwork) = connection.T;
                poreVolumes(eggRealization, rowInFullNetwork) = connection.pv/baseNetworkModel.options.cellsPerConnection;
            end

            fprintf('Flow diagnostic analysis %d%% complete\n', floor(eggRealization*100/ensembleSize));

        end

        if ~exist(initParamFolder, 'dir')
            mkdir(initParamFolder);
        end
        save(initParamFile, 'transmissibilities', 'poreVolumes', ...
            'wellProductionIndicesSum', 'wellProductionIndicesMean');
    end


    %% Create sample object from initial ensemble
    % We store the initial ensemble parameters in the appropriate Sample
    % objects designed for history matching through the ensemble module. For 
    % this, we use a class that are designed specifically for network models
    % for holding the relevant parameters. We also specify maximum and minimum
    % allowed values to avoid that simulations break down due to nonphysical
    % configurations.

    % The parameters transmissibility and pore volume belongs to the operators
    % of a model.
    operatorData = cell(ensembleSize, 1);
    for i = 1:ensembleSize
        operatorData{i}.pv = poreVolumes(i, :);
        operatorData{i}.T = transmissibilities(i, :);
    end
    transMax = opt.transMaxFactor*max(max(transmissibilities));
    pvMax = opt.transMaxFactor*max(max(poreVolumes));

    connectionIndices = struct('faces', {baseNetworkModel.options.networkModel.graph.Edges.faceIx}, ...
                               'cells', {baseNetworkModel.options.networkModel.graph.Edges.cellIx});
    operatorSamples = NetworkOperatorSamplesHM('data', operatorData, ...
        'connectionIndices', connectionIndices, ...
        ... % parameter scaling 
        'pvScale', opt.pvScale, 'TScale', opt.TScale, ...
        'minPvValue', opt.pvMin, 'maxPvValue', pvMax, ...
        'minTValue', opt.transMin, 'maxTValue', transMax);

    % Well production indices is considered a well parameter
    wellSampleData = cell(ensembleSize, 1);

    % We choose either the mean or the sum of the original well production
    % indices
    switch opt.WItype
        case 'mean'
            wellProductionIndices = wellProductionIndicesMean;
        case 'sum'
            wellProductionIndices = wellProductionIndicesSum;
        otherwise
                error('\nWItype %s is not implemented\n', opt.WItype);    
    end
    for i = 1:ensembleSize
        wellSampleData{i}.WI = wellProductionIndices(i, :);
    end
    minWI = opt.minWIFactor*min(min(wellProductionIndices(:,:)));
    maxWI = opt.maxWIFactor*max(max(wellProductionIndices(:,:)));

    wellSamples = WellSamplesHM('data', wellSampleData, ...
                                'WIScale', opt.WIScale, ...
                                'minWIValue', minWI, 'maxWIValue', maxWI);

    % Wrap the two sample objects in a single composite object
    samples = CompositeSamplesHM({operatorSamples, wellSamples});
end

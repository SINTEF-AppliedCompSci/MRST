function samples = priorSamplesForUpscaledEggModel(baseUpscaledModel, ...
                                                   referenceExample, ...
                                                   regenerateInitialEnsemble, ...
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

    options = struct('eggRealizations', [1:100]);
    [options, extra] = merge_options(options, varargin{:});

    ensembleSize = numel(options.eggRealizations);

    pvMin = min(baseUpscaledModel.model.operators.pv)/1000;
    transMin = eps;
    numInternalFaces = numel(baseUpscaledModel.model.operators.T);
    numCells = baseUpscaledModel.model.G.cells.num;
    numWells = numel(baseUpscaledModel.schedule.control(1).W);

    initParamFolder = fullfile(mrstOutputDirectory(), ...
                             'eggEnsembleCoarseESMDACalibration_examples2');
    initParamFile = fullfile(initParamFolder, ...
                             ['upscaled_parameters.mat']);

    % To save processing time, we store the initial parameters in a .mat file
    % so that we don't need to do the upscaling for each time we run this
    % example.
    if exist(initParamFile, 'file') && ~regenerateInitialEnsemble

        fprintf('Read pre-computed initial ensemble from file\n');
        % Read pre-computed initial ensemble
        preComputedEnsemble = load(initParamFile);
        transmissibilities = preComputedEnsemble.transmissibilities;
        poreVolumes = preComputedEnsemble.poreVolumes;
        pertPoreVolumes = preComputedEnsemble.pertPoreVolumes;
        %wellProductionIndices = preComputedEnsemble.wellProductionIndices;

        wellProductionIndices = preComputedEnsemble.WImean;

        % Read transmissibilities from the old prior    
        %oldParams = load('C:\Users\havardh\Documents\MATLAB\mrst-bitbucket\EggFullUpscaleEnsemble.mat');
        %transmissibilities = oldParams.Transmisibility';
        %poreVolumes = oldParams.PoreVolume';
        %wellProductionIndices = oldParams.WI_mean';

    else
        % Allocate matrices for temporal storage of ensemble parameters
        poreVolumes = ones(ensembleSize, numCells).*pvMin;
        pertPoreVolumes = ones(ensembleSize, numCells).*pvMin;
        wellProductionIndices = zeros(ensembleSize, numWells);
        transmissibilities = ones(ensembleSize, numInternalFaces).*transMin;

        WIsum = zeros(ensembleSize, numWells);
        WImean = zeros(ensembleSize, numWells);
        WIharmonic = zeros(ensembleSize, numWells);


        totalVolume = sum(baseUpscaledModel.model.operators.pv);
        fractionsOfTotalVolume = baseUpscaledModel.model.operators.pv/totalVolume;


        for eggRealization = options.eggRealizations
            fullRealization = TestCase('egg_wo', 'realization', eggRealization);

            % Use the same time steps as the reference model, but keep the well
            % configurations
            fullRealization.schedule = simpleSchedule(referenceExample.schedule.step.val, ...
                                                      'W', fullRealization.schedule.control.W);

            %fullRealization.schedule.control.W(1).WI

            coarseRealization = TestCase('upscaled_coarse_network', ...
                                            'partition', [6 6 1], ...
                                            'referenceCase', fullRealization, ...
                                            'plotCoarseModel', false);

            poreVolumes(eggRealization, :) = coarseRealization.model.operators.pv(:);
            transmissibilities(eggRealization, :) = coarseRealization.model.operators.T(:);
            wellProductionIndices(eggRealization, :) = [coarseRealization.schedule.control.W.WI]';

            % Sample random numbers for perturbing pore volumes      
            pvPerturbations = randn(numCells,1).*sqrt(coarseRealization.model.operators.pv);

            % Ensure that the perturbation does not alter the total available
            % reservoir volume
            pertVolume = sum(pvPerturbations);
            pvPerturbations = pvPerturbations - fractionsOfTotalVolume.*pertVolume;
            pertPoreVolumes(eggRealization, :) = pvPerturbations;

            % Store the different well index upscaling methods
            WIsum(eggRealization, :) = [coarseRealization.options.sum]';
            WImean(eggRealization, :) = [coarseRealization.options.mean]';
            WIharmonic(eggRealization, :) = [coarseRealization.options.harmonic]';

            fprintf('Upscaling realization %d%% complete\n', eggRealization);

        end     

        % Store parameters:
        if ~exist(initParamFolder, 'dir')
            mkdir(initParamFolder);
        end
        save(initParamFile, 'transmissibilities', 'wellProductionIndices', ...
             'poreVolumes', 'pertPoreVolumes', ...
             'WIsum', 'WImean', 'WIharmonic');
    end

    %% Perturb pore volumes 
    % Add the perturbations to the pore volumes at the proper scale.
    pertStd = 10;
    poreVolumes = poreVolumes + pertPoreVolumes*pertStd;

    %% Create sample object from initial ensemble
    % We store the initial ensemble parameters in the appropriate Sample
    % objects designed for history matching through the ensemble module. Both
    % pore volumes and transmissibilities belong to the operator sample class,
    % whereas well productivity index is found in well samples. We also specify 
    % maximum and minimum allowed values to avoid that simulations break down 
    % due to nonphysical configurations.

    % The parameters transmissibility and pore volume belongs to the operators
    % of a model.
    operatorData = cell(ensembleSize, 1);
    for i = 1:ensembleSize
        operatorData{i}.pv = poreVolumes(i, :)';
        operatorData{i}.T = transmissibilities(i, :)';
    end
    transMax = 10*max(max(transmissibilities));
    pvMax = 10*max(max(poreVolumes));
    operatorSamples = OperatorSamplesHM('data', operatorData, ...
        ... % parameter scaling 
        'pvScale', 1e4, 'TScale', 1e-9, ...
        'minPvValue', pvMin, 'maxPvValue', pvMax, ...
        'minTValue', transMin, 'maxTValue', transMax);

    % Well production indices is considered a well parameter
    wellSampleData = cell(ensembleSize, 1);

    for i = 1:ensembleSize
        wellSampleData{i}.WI = wellProductionIndices(i, :)*2;
    end
    minWI = 0.01*min(min(wellProductionIndices(:,:)));
    maxWI = 8*max(max(wellProductionIndices(:,:)));

    wellSamples = WellSamplesHM('data' ,wellSampleData, ...
                                'WIScale', 1e-11, ...
                                'minWIValue', minWI, 'maxWIValue', maxWI);

    % Wrap the two sample objects in a single composite object
    samples = CompositeSamplesHM({operatorSamples, wellSamples});



end


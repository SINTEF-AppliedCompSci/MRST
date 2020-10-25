classdef RockSamples < BaseSamples

    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Get model
            model  = problem.SimulatorSetup.model;
            rmodel = getReservoirModel(model);
            % Set rock properties from sample data
            rmodel = samples.mapRockProps(rmodel, sampleData);
            % Set up operators
            rmodel = rmodel.setupOperators();
            % Replace problem model
            model = setReservoirModel(model, rmodel);
            problem.SimulatorSetup.model = model;
        end
        
        %-----------------------------------------------------------------%
        function model = mapRockProps(samples, model, sampleData)
            [sampleData, type] = samples.validateSample(model, sampleData);
            switch type
                case 'standard'
                    % Rock sample has one value per grid cell - simply replace
                    % rock data with sample data
                    model.rock.perm = sampleData.perm;
                    model.rock.poro = sampleData.poro;
                case 'upscale'
                    % Permeability and porosity are defined on an
                    % underlying fine grid - upscale to model.G
                    % Averaging matrix
                    S = sparse(model.G.partition, 1:nc, 1./accumarray(partition, 1));
                    % Average log10(permeability)
                    model.rock.perm = zeros(nc, size(sampleData.perm,2));
                    for i = 1:size(sampleData.perm,2)
                        lperm = S*log10(sampleData.rock.perm(:,i));
                        model.rock.perm(:,i) = 10.^lperm;
                    end
                    % Average porosity
                    model.rock.poro = S*sampleData.poro;
                case 'sampleFromBox'
                    % Rock sample is given as array, assumed to cover the
                    % bounding box of model.G - sample values
                    % Sample permeability
                    model.rock.perm = zeros(model.G.cells.num, numel(sampleData.perm));
                    for i = 1:numel(sampleData.perm)
                        model.rock.perm(:,i) = sampleFromBox(model.G, sampleData.perm{i});
                    end
                    % Sample porosity
                    model.rock.poro = sampleFromBox(model.G, sampleData.poro);
            end
        end
        
        %-----------------------------------------------------------------%
        function [sampleData, type] = validateSample(samples, model, sampleData)
            % Check that sample has perm and poro field
            assert(isfield(sampleData, 'perm') && isfield(sampleData, 'poro'), ...
                'Rock sample must have a perm and poro field');
            if size(sampleData.poro, 2) == 1
                % We are given a vector of cell values
                nc = numel(sampleData.poro);
                if nc < model.G.cells.num
                    % Sample data defined on an underlying fine grid. Check
                    % that we have the coarse grid information
                    assert(isfield(model.G, 'partition'), ['Data sample size ', ...
                        'is greater than G.cells.num, but G is not a coarse grid']);
                    type = 'upscale';
                else
                    % Data is defined directly on the grid
                    type = 'standard';
                end
                % Permeability should be [nc, 1] or [nc, G.griddim]
                assert(all(size(sampleData.perm) == [nc, 1]) ||             ...
                       all(size(sampleData.perm) == [nc, model.G.griddim]), ...
                    'Inconsistent permeability dimensions');
                % Porosity should be [nc, 1]
                assert(numel(sampleData.poro) == model.G.cells.num, ...
                    'Inconsistent porosity dimensions');
            else
                % Sample data is defined on in a cube assumed to be the
                % bounding box of model.G
                if ~iscell(sampleData.perm)
                    sampleData.perm = {sampleData.perm};
                end
                assert((numel(sampleData.perm) == 1 || ...
                        numel(sampleData.perm) == model.G.griddim) && ...
                        all(cellfun(@(perm) ndims(perm) == model.G.griddim, sampleData.perm)), ...
                    'Inconsistent permeability dimensions');
                assert(ndims(sampleData.poro) == model.G.griddim, ...
                    'Inconsistent porosity dimensions');
                type = 'sampleFromBox';
            end
        end
        
        
        %-----------------------------------------------------------------%
        % Functions related to history matching
        %-----------------------------------------------------------------%
        function sampleVectors = getSampleVectors(samples)
            
            assert(~isempty(samples.data), ...
                'This function only work with pre-generated samples');
            
            assert(numel(samples.data{1}.poro) == numel(samples.data{1}.perm), ...
                'Inconsistent permeability and porosity sizes');
            
            numCells = numel(samples.data{1}.poro);
            sampleVectors = zeros(2*numCells, samples.num);
            
            for i = 1:samples.num
                if samples.transformSampleVectors
                    % Logarithmic transform of permeability 
                    perm = convertTo(samples.data{i}.perm(:), milli*darcy);
                    sampleVectors(1:numCells, i) = log(perm);
                else
                    sampleVectors(1:numCells, i) = samples.data{i}.perm(:);
                end
                sampleVectors(numCells+1:numCells*2, i) = samples.data{i}.poro(:);
            end
        end
        
        function samples = setSampleVectors(samples, newSampleVectors)
            
            assert(size(newSampleVectors, 2) == samples.num, ...
                'number of columns of new samples does not match ensemble size');
            numCells = numel(samples.data{1}.poro);
            assert(size(newSampleVectors, 1) == numCells*2, ...
                'number of rows of new sample does not match old sample size');
            
            for i = 1:samples.num
                if samples.transformSampleVectors
                    perm = exp(newSampleVectors(1:numCells, i));
                    samples.data{i}.perm(:) = convertFrom(perm, milli*darcy);
                else
                    samples.data{i}.perm(:) = newSampleVectors(1:numCells, i);
                end
                samples.data{i}.poro(:) = newSampleVectors(numCells+1:numCells*2, i);
            end
            
        end
           
        
    end
    
end
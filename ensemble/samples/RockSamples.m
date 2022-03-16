classdef RockSamples < BaseSamples
    % Class for holding stochastic samples that represent uncertain rock
    % properties for an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize stochastic rock properties. These rock 
    %   properties can either be pre-computed or generated on the fly.
    % 
    % SYNOPSIS
    %   samples = RockSamples('data', data);
    %   samples = RockSamples('generatorFn', generatorFn);
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
    %   Each sample should consist of a struct with the fields 'perm' and 
    %   'poro'.
    %
    % SEE ALSO:
    %   `WellSamples`, `DeckSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        updateWells = true % Update well indices based on updated rock 
        WIargs      = {}   % Extra arguments passed to `computeWellIndex`
    end
    
    methods
        %-----------------------------------------------------------------%
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample realization of the rock properties to a 
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
            if samples.updateWells
                % Update well indices
                schedule = problem.SimulatorSetup.schedule;
                schedule = samples.updateWellIndices(rmodel, schedule);
                problem.SimulatorSetup.schedule = schedule;
            end
        end
    end
    
    methods (Access = protected)
        %-----------------------------------------------------------------%
        function model = mapRockProps(samples, model, sampleData)
            % Utility function for mapping the sample data to the actuel
            % rock model.
            
            [sampleData, type] = samples.validateSample(model, sampleData);
            switch type
                case 'standard'
                    % Rock sample has one value per grid cell - simply replace
                    % rock data with sample data
                    model.rock.perm = sampleData.perm;
                    model.rock.poro = sampleData.poro;
                    if isfield(sampleData, 'ntg')
                        model.rock.ntg = sampleData.ntg;
                    end
                    if isfield(sampleData, 'multiplies')
                        model.rock.multiplies = sampleData.multiplies;
                    end
                case 'upscale'
                    % Permeability and porosity are defined on an
                    % underlying fine grid - upscale to model.G
                    % Averaging matrix
                    p  = model.G.partition;
                    n  = accumarray(p, 1);
                    S = sparse(p, 1:numel(p), 1./n(p));
                    % Average log10(permeability)
                    model.rock.perm = zeros(max(p), size(sampleData.perm,2));
                    for i = 1:size(sampleData.perm,2)
                        lperm = S*log10(sampleData.perm(:,i));
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
            % Utility function for checking that th sample data is a valid
            % rock property, and evaluate how it should be mapped to the
            % model from the base problem.
            
            % Check that sample has perm and poro field
            assert(isfield(sampleData, 'perm') && isfield(sampleData, 'poro'), ...
                'Rock sample must have a perm and poro field');
            if size(sampleData.poro, 2) == 1
                % We are given a vector of cell values
                nc = numel(sampleData.poro);
                assert(size(sampleData.perm, 1) == nc, ['Rock sample ', ...
                            'perm/poro must have the same number of rows']);
                if nc > model.G.cells.num
                    % Sample data defined on an underlying fine grid. Check
                    % that we have the coarse grid information
                    assert(isfield(model.G, 'partition'), ['Data sample size '    , ...
                        'is greater than G.cells.num, but G is not a coarse grid']);
                    assert(numel(model.G.partition) == nc,['Data sample size '    , ...
                        'does not match underlying fine grid size']               );
                    type = 'upscale';
                else
                    % Data is defined directly on the grid
                    if nc == 1
                        % Homogeneous poro/perm
                        sampleData.poro = repmat(sampleData.poro, model.G.cells.num, 1);
                        sampleData.perm = repmat(sampleData.perm, model.G.cells.num, 1);
                        nc              = model.G.cells.num;
                    end
                    type = 'standard';
                end
                % Permeability should be [nc, 1] or [nc, G.griddim]
                assert(all(size(sampleData.perm) == [nc, 1]) ||             ...
                       all(size(sampleData.perm) == [nc, model.G.griddim]), ...
                    'Inconsistent permeability dimensions');
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
        function schedule = updateWellIndices(samples, model, schedule)
            % Update well indices for all wells based on sample rock
            % properties.
            updateWI = @(W) computeWellIndex(model.G, model.rock, W.r, W.cells, ...
                                             'Dir', W.dir, samples.WIargs{:});
            for i = 1:numel(schedule.control)
                if isfield(schedule.control(i), 'W')
                    W = schedule.control(i).W;
                    for j = 1:numel(W)
                        W(j).WI = updateWI(W(j));
                    end
                    schedule.control(i).W = W;
                end
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
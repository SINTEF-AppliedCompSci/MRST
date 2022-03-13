classdef DeckSamples < BaseSamples
    % Class that holds an ensemble of decks used for initializing ensemble
    % members in an MRSTEnsemble.
    %
    % DESCRIPTION:
    %   This class is an extention of `BaseSamples` and is used within a 
    %   MRSTEnsemble to organize ensemble members initialized from decks.
    %   These decks can either be pre-computed or generated on the fly.
    % 
    % SYNOPSIS
    %   samples = DeckSamples('data', data);
    %   samples = DeckSamples('generatorFn', generatorFn);
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
    %   'gridFromDeck' - Flag that indicates whether the grid itself is
    %                    also included in the deck samples. Default: true.
    %
    %   'initArgs' - Arguments passed directly to initEclipseProblemAD.
    % 
    % NOTE:
    %   Either 'data' or 'generatorFn' must be provided.
    %   Each sample should consist of a struct with the fields that are
    %   found in a DECK struct.
    %
    % SEE ALSO:
    %   `WellSamples`, `RockSamples`, `BaseSamples`, `MRSTExample`, `BaseQoI`
    
    properties
        gridFromDeck = true % Construct MRST grid from deck
        initArgs     = {}   % Arguments passed directly to initEclipseProblemAD
    end
    
    methods
        function problem = setSample(samples, sampleData, problem)
            % Applies the sample deck to a problem
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
            
            if samples.gridFromDeck
                % Construct MRST grid from deck
                G = [];
            else
                % Use grid from problem
                G = problem.simulatorSetup.model.G;
            end
            % Get initial state, model and schedule
            [state0, model, schedule] = initEclipseProblemAD(sampleData, ...
                                                 'G', G, samples.initArgs{:});
            % update problem
            problem.SimulatorSetup.state0   = state0;
            problem.SimulatorSetup.model    = model;
            problem.SimulatorSetup.schedule = schedule;
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
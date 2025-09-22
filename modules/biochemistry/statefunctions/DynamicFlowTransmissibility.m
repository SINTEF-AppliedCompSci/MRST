classdef DynamicFlowTransmissibility < StateFunction
    % Dynamic transmissibility computation accounting for flow properties
    %
    % SYNOPSIS:
    %   dt = DynamicFlowTransmissibility(model, 'conductivity_name', 'Permeability')
    %
    % DESCRIPTION:
    %   Computes dynamic transmissibilities based on conductivity fields
    %   (e.g., permeability or mobility) using two-point flux approximation
    %   with harmonic averaging.
    %
    % REQUIRED PARAMETERS:
    %   model - Reservoir model with grid and operators
    %
    % OPTIONAL PARAMETERS:
    %   'conductivity_name' - Name of conductivity property (default: 'Permeability')
    %
    % RETURNS:
    %   Class instance for transmissibility calculation
    %
    % SEE ALSO:
    %   ReservoirModel, TwoPointFluxApproximation

    properties
        harmonicAvgOperator    % Function handle for harmonic averaging
        twoPointOperator       % Function handle for two-point approximation
        conductivity_name      % Name of conductivity property to use
    end

    methods
        function dt = DynamicFlowTransmissibility(model, conductivity_name)
            % Constructor for dynamic transmissibility calculator

            % Set default conductivity property if not provided
            if nargin < 2
                conductivity_name = 'Permeability';
            end

            % Initialize base class
            dt@StateFunction(model);

            % Set up operators
            dt.twoPointOperator    = getTwoPointOperator(model.G);
            dt.harmonicAvgOperator = getHarmonicAvgOperator(model.G);
            dt.conductivity_name   = conductivity_name;

            % Declare dependencies
            dt = dt.dependsOn(conductivity_name);

            dt.label = 'T'; % Transmissibility label
        end

        function T = evaluateOnDomain(prop, model, state, allFaces)
            % Compute transmissibilities for current state
            %
            % PARAMETERS:
            %   prop      - Property function instance
            %   model     - Reservoir model instance
            %   state     - State struct containing fields
            %   allFaces  - Boolean, return all faces if true (optional)
            %
            % RETURNS:
            %   T - Transmissibility values for connections

            % Get conductivity field (permeability or mobility)
            lambda = prop.getEvaluatedDependencies(state, prop.conductivity_name);

            % Handle both cell and array inputs
            if iscell(lambda)
                T = cell(numel(lambda), 1);
                for i = 1:numel(lambda)
                    if ~isempty(lambda{i})
                        T{i} = prop.getTransmissibility(lambda{i});
                    end
                end
            else
                T = prop.getTransmissibility(lambda);

                % Return only internal connections by default
                if nargin < 4 || ~allFaces
                    T = T(model.operators.internalConn);
                end
            end
        end

        function T = getTransmissibility(prop, lambda)
            % Compute transmissibility from conductivity field
            %
            % PARAMETERS:
            %   lambda - Conductivity field (permeability/mobility)
            %
            % RETURNS:
            %   T - Transmissibility values

            % Apply two-point flux approximation
            T = prop.twoPointOperator(lambda);

            % Apply harmonic averaging
            T = prop.harmonicAvgOperator(T);

            % Ensure positive transmissibility
            T = abs(T);
        end
    end
end

function tp = getTwoPointOperator(G)
% Create two-point flux approximation operator
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   tp - Function handle for two-point approximation

% Mappings from cells to faces
cells = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
faces = G.cells.faces(:,1);

% Vector from cell to face centroids
C = G.faces.centroids(faces,:) - G.cells.centroids(cells,:);

% Oriented normals
sgn = 2*(cells == G.faces.neighbors(faces, 1)) - 1;
N   = bsxfun(@times, sgn, G.faces.normals(faces, :));

% Compute two-point weights
cn = sum(C.*N, 2)./sum(C.*C, 2);

% Return function handle
tp = @(lambda) cn.*lambda(cells);
end

function ha = getHarmonicAvgOperator(G)
% Create harmonic averaging operator
%
% PARAMETERS:
%   G - Grid structure
%
% RETURNS:
%   ha - Function handle for harmonic averaging

faces = G.cells.faces(:,1);
M = sparse(faces, 1:numel(faces), 1, G.faces.num, numel(faces));
ha = @(T) 1./(M*(1./T));
end

%{
Copyright 2009-2025 SINTEF Digital, Mathematics & Cybernetics.

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
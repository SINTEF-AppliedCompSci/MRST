function [vW, vM, vO, vU, model] = ...
        getFluxAndPropsMICP(model, pW, m, o, u, gdz, poro)
% Function to update permeability and compute water and component fluxes.
%
% This function is modified from a file in the MATLAB Reservoir Simulation
% Toolbox (MRST), see
%   mrst/modules/ad-eor/utils/getFluxAndPropsWaterPolymer_BO.m
%
% We refer to that function for a complete commented version of the file.

%{
Partial copyright 2009-2021, SINTEF Digital, Mathematics & Cybernetics.
Partial copyright 2021-2026, NORCE Research AS, Computational Geosciences and
Modeling.

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file. If not, see <http://www.gnu.org/licenses/>.
%}
    model = updateRockAndOperators(model, poro);

    [vW, vM, vO, vU] = calculateFluxes( ...
        model, pW, m, o, u, gdz);
end

function model = updateRockAndOperators(model, poro)
    % Compute current permeability
    permeability = model.fluid.K(poro);
    permeability = max(permeability, model.fluid.kmin);

    % Update rock with the current permeability and porosity
    model.rock = makeRock(model.G, permeability, poro);

    operators = model.operators;

    if ~isfield(operators, 'micpHalfFaceToFace')
        faceIndices = model.G.cells.faces(:, 1);
        numberOfHalfFaces = numel(faceIndices);

        operators.micpHalfFaceToFace = sparse( ...
            faceIndices, ...
            1 : numberOfHalfFaces, ...
            1, ...
            model.G.faces.num, ...
            numberOfHalfFaces);

        operators.micpInternalFaces = ...
            all(model.G.faces.neighbors ~= 0, 2);
    end

    halfTransmissibility = ...
        computeHalfTransmissibility(model.G, model.rock);

    inverseHalfTransmissibility = ...
        1 ./ halfTransmissibility;

    inverseFaceTransmissibility = ...
        operators.micpHalfFaceToFace * ...
        inverseHalfTransmissibility;

    faceTransmissibility = ...
        1 ./ inverseFaceTransmissibility;

    operators.T_all = faceTransmissibility;
    operators.T = faceTransmissibility( ...
        operators.micpInternalFaces);
    operators.pv = model.G.cells.volumes .* poro;

    model.operators = operators;
end

function halfTransmissibility = ...
        computeHalfTransmissibility(grid, rock)

    if hasCornerPointGeometry(grid)
        cornerPointGeometry = grid.cells.cpgeometry;

        halfTransmissibility = computeTrans( ...
            grid, ...
            rock, ...
            'K_system', ...
            'loc_xyz', ...
            'cellCenters', ...
            cornerPointGeometry.centroids, ...
            'cellFaceCenters', ...
            cornerPointGeometry.facecentroids);
    else
        halfTransmissibility = computeTrans(grid, rock);
    end
end

function present = hasCornerPointGeometry(grid)
    if ~isfield(grid, 'nodes')
        present = false;
        return
    end

    physicalDimension = size(grid.nodes.coords, 2);
    numberOfCells = grid.cells.num;
    numberOfCellFaces = size(grid.cells.faces, 1);

    present = ...
        isfield(grid.cells, 'cpgeometry') && ...
        all(isfield( ...
        grid.cells.cpgeometry, ...
        {'centroids', 'facecentroids'}));

    if ~present
        return
    end

    cornerPointGeometry = grid.cells.cpgeometry;

    present = ...
        ismatrix(cornerPointGeometry.centroids) && ...
        ismatrix(cornerPointGeometry.facecentroids) && ...
        isequal( ...
        size(cornerPointGeometry.centroids), ...
        [numberOfCells, physicalDimension]) && ...
        isequal( ...
        size(cornerPointGeometry.facecentroids), ...
        [numberOfCellFaces, physicalDimension]);
end

function [vW, vM, vO, vU] = ...
        calculateFluxes(model, pW, m, o, u, gdz)

    fluid = model.fluid;
    operators = model.operators;
    waterMobility = 1 ./ fluid.muw;

    waterPotentialDifference = ...
        operators.Grad(pW) - fluid.rhoWS .* gdz;
    upstreamWater = value(waterPotentialDifference) <= 0;

    % Water
    vW = -waterMobility .* operators.T .* ...
        waterPotentialDifference;

    % Microorganisms
    vM = getComponentFlux( ...
        operators, upstreamWater, m, vW);

    % Oxygen
    vO = getComponentFlux( ...
        operators, upstreamWater, o, vW);

    % Urea
    vU = getComponentFlux( ...
        operators, upstreamWater, u, vW);
end

function componentFlux = getComponentFlux( ...
        operators, upstreamWater, concentration, waterFlux)

    [faceConcentration, ~] = ...
        operators.splitFaceCellValue( ...
        operators, upstreamWater, concentration);

    componentFlux = faceConcentration .* waterFlux;
end

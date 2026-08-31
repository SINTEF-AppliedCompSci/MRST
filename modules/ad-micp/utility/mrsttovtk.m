function mrsttovtk(G, states, name, formatn)
%  Write one or more data sets into the files name.pvd and name-n.vtu
%  for visualization in ParaView.
%
% SYNOPSIS:
%       mrsttovtk(G, states, name, formatn)
%
% PARAMETERS:
%   G       - Grid data structure.
%
%   states  - Structure or cell array of structures containing cell data.
%
%   name    - Name of the output files name.pvd and name-n.vtu.
%
%   formatn - Format for writing numerical values.
%
% NOTES:
%   The function `mrsttovtk` only writes cell values (e.g., pressures), but
%   currently not values on faces (e.g., fluxes). This function supports 3D
%   tetrahedral, polyhedral, and corner-point grids, 2D Cartesian, triangular,
%   and polygonal grids, and 1D grids.
%
% EXAMPLE:
%   Let us consider the grid from the MRST 2019 book (pp. 97, Fig. 3.31):
%   G = simpleGrdecl([20, 20, 5], @(x) .25*(x-.5), 'flat', true);
%   G = processGRDECL(G);
%
%   % 1) Writing the content in G.cells.indexMap into the files states.pvd
%   and states-00001.vtu for visualization in ParaView (since indexMap
%   is the only field with the same number of values as grid cells).
%   mrsttovtk(G, G.cells, 'states', '%.6g');
%
%   % 2) After making a structure with 100 values, we run the mrsttovtk
%   function to write every tenth value for visualization in ParaView.
%   for i = 1 : 100
%       states{i}.index = i^2*G.cells.indexMap;
%   end
%   mrsttovtk(G, {states{10:10:100}}, 'states', '%.6g');

%{
Copyright 2021-2026, NORCE Research AS, Computational Geosciences and
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
    valueFormat = formatn;
    spacedValueFormat = [formatn ' '];
    geometryText = cell(0, 1);

    geometryText{end + 1} = sprintf('\t\t\t<Points>\n');
    geometryText{end + 1} = sprintf([ ...
        '\t\t\t\t<DataArray type="Float32" Name="Coordinates" ' ...
        'NumberOfComponents="3" format="ascii">\n']);
    geometryText{end + 1} = sprintf('\t\t\t\t\t');

    coordinates = zeros(G.nodes.num, 3);
    coordinates(:, 1 : G.griddim) = ...
        G.nodes.coords(:, 1 : G.griddim);

    if G.griddim == 1
        vtkCellType = 3;
    elseif G.griddim == 2
        if 4*G.cells.num == size(G.cells.faces, 1)
            vtkCellType = 8;
        else
            vtkCellType = 7;
        end
    else
        vtkCellType = 42;
    end

    if G.nodes.num > 1
        coordinateFormat = repmat( ...
            spacedValueFormat, 1, 3*(G.nodes.num - 1));
        geometryText{end + 1} = sprintf( ...
            coordinateFormat, ...
            reshape(coordinates(1 : end - 1, :)', [], 1));
    end

    geometryText{end + 1} = sprintf( ...
        [spacedValueFormat spacedValueFormat valueFormat '\n'], ...
        coordinates(end, :));
    geometryText{end + 1} = sprintf('\t\t\t\t</DataArray>\n');
    geometryText{end + 1} = sprintf('\t\t\t</Points>\n');
    geometryText{end + 1} = sprintf('\t\t\t<Cells>\n');
    geometryText{end + 1} = sprintf([ ...
        '\t\t\t\t<DataArray type="UInt8" Name="types" ' ...
        'NumberOfComponents="1" format="ascii">\n']);
    geometryText{end + 1} = sprintf('\t\t\t\t\t');

    if G.cells.num > 1
        geometryText{end + 1} = sprintf( ...
            '%d ', repmat(vtkCellType, 1, G.cells.num - 1));
    end

    geometryText{end + 1} = sprintf('%d\n', vtkCellType);
    geometryText{end + 1} = sprintf('\t\t\t\t</DataArray>\n');

    cellConnectivity = cell(G.cells.num, 1);
    faceOffsetLengths = zeros(G.cells.num, 1);
    faceNodes = cell(G.faces.num, 1);

    for face = 1 : G.faces.num
        faceNodeRows = ...
            G.faces.nodePos(face) : G.faces.nodePos(face + 1) - 1;
        faceNodes{face} = G.faces.nodes(faceNodeRows, 1)' - 1;
    end

    if G.griddim == 3
        faceText = cell(G.faces.num, 1);

        for face = 1 : G.faces.num
            nodes = faceNodes{face};
            faceText{face} = sprintf( ...
                repmat(' %d', 1, numel(nodes) + 1), ...
                numel(nodes), nodes);
        end

        geometryText{end + 1} = sprintf([ ...
            '\t\t\t\t<DataArray type="Int32" Name="faces" ' ...
            'NumberOfComponents="1" format="ascii">\n']);
        geometryText{end + 1} = sprintf('\t\t\t\t\t');

        for cellIndex = 1 : G.cells.num
            cellFaceRows = G.cells.facePos(cellIndex) : ...
                G.cells.facePos(cellIndex + 1) - 1;
            cellFaces = G.cells.faces(cellFaceRows, 1)';

            if cellIndex == 1
                geometryText{end + 1} = sprintf( ...
                    '%d', numel(cellFaces));
            else
                geometryText{end + 1} = sprintf( ...
                    ' %d', numel(cellFaces));
            end

            cellNodes = -1;
            numberOfFaceEntries = 0;

            for face = cellFaces
                nodes = faceNodes{face};
                geometryText{end + 1} = faceText{face};
                cellNodes = unique([nodes cellNodes]);
                numberOfFaceEntries = ...
                    numberOfFaceEntries + numel(nodes);
            end

            cellConnectivity{cellIndex} = cellNodes(2 : end);
            faceOffsetLengths(cellIndex) = ...
                numberOfFaceEntries + numel(cellFaces) + 1;
        end

        geometryText{end + 1} = sprintf('\n');
        geometryText{end + 1} = sprintf('\t\t\t\t</DataArray>\n');
    else
        for cellIndex = 1 : G.cells.num
            cellFaceRows = G.cells.facePos(cellIndex) : ...
                G.cells.facePos(cellIndex + 1) - 1;
            cellFaces = G.cells.faces(cellFaceRows, 1)';
            connectivity = zeros(numel(cellFaces), 1);
            nodes = faceNodes{cellFaces(1)};
            cellNodes = unique([nodes -1], 'stable');
            newNodes = setdiff(cellNodes, -1);

            connectivity(1) = newNodes(1);
            connectivity(2) = newNodes(2);

            numberOfFaceEntries = numel(nodes);
            connectivityIndex = 3;

            for faceIndex = 2 : numel(cellFaces) - 1
                face = cellFaces(faceIndex);
                nodes = faceNodes{face};
                updatedCellNodes = unique( ...
                    [nodes cellNodes], 'stable');
                newNodes = setdiff(updatedCellNodes, cellNodes);
                connectivity(connectivityIndex) = newNodes;
                cellNodes = updatedCellNodes;
                numberOfFaceEntries = ...
                    numberOfFaceEntries + numel(nodes);
                connectivityIndex = connectivityIndex + 1;
            end

            cellConnectivity{cellIndex} = connectivity;
            faceOffsetLengths(cellIndex) = ...
                numberOfFaceEntries + numel(cellFaces) + 1;
        end
    end

    connectivityLengths = cellfun(@numel, cellConnectivity);
    vtkOffsets = cumsum(connectivityLengths);
    connectivityColumns = cell(size(cellConnectivity));

    for cellIndex = 1 : numel(cellConnectivity)
        connectivityColumns{cellIndex} = ...
            cellConnectivity{cellIndex}(:);
    end

    vtkConnectivity = vertcat(connectivityColumns{:});

    geometryText{end + 1} = sprintf([ ...
        '\t\t\t\t<DataArray type="Int32" Name="offsets" ' ...
        'NumberOfComponents="1" format="ascii">\n']);
    geometryText{end + 1} = sprintf('\t\t\t\t\t');

    if numel(vtkOffsets) > 1
        geometryText{end + 1} = sprintf( ...
            '%d ', vtkOffsets(1 : end - 1));
    end

    geometryText{end + 1} = sprintf('%d\n', vtkOffsets(end));
    geometryText{end + 1} = sprintf('\t\t\t\t</DataArray>\n');
    geometryText{end + 1} = sprintf([ ...
        '\t\t\t\t<DataArray type="Int32" Name="connectivity" ' ...
        'NumberOfComponents="1" format="ascii">\n']);
    geometryText{end + 1} = sprintf('\t\t\t\t\t');

    if numel(vtkConnectivity) > 1
        geometryText{end + 1} = sprintf( ...
            '%d ', vtkConnectivity(1 : end - 1));
    end

    geometryText{end + 1} = sprintf('%d\n', vtkConnectivity(end));
    geometryText{end + 1} = sprintf('\t\t\t\t</DataArray>\n');

    if G.griddim == 3
        vtkFaceOffsets = cumsum(faceOffsetLengths);

        geometryText{end + 1} = sprintf([ ...
            '\t\t\t\t<DataArray type="Int32" Name="faceoffsets" ' ...
            'NumberOfComponents="1" format="ascii">\n']);
        geometryText{end + 1} = sprintf('\t\t\t\t\t');

        if numel(vtkFaceOffsets) > 1
            geometryText{end + 1} = sprintf( ...
                '%d ', vtkFaceOffsets(1 : end - 1));
        end

        geometryText{end + 1} = sprintf( ...
            '%d\n', vtkFaceOffsets(end));
        geometryText{end + 1} = sprintf( ...
            '\t\t\t\t</DataArray>\n');
    end

    geometryText{end + 1} = sprintf('\t\t\t</Cells>\n');
    geometryText{end + 1} = sprintf('\t\t</Piece>\n');
    geometryText{end + 1} = sprintf('\t</UnstructuredGrid>\n');
    geometryText{end + 1} = sprintf('</VTKFile>\n');
    geometryText = [geometryText{:}];

    if iscell(states)
        stateList = states;
    else
        stateList = {states};
    end

    initialState = stateList{1};
    fields = fieldnames(initialState);
    keepField = false(numel(fields), 1);

    for fieldIndex = 1 : numel(fields)
        values = initialState.(fields{fieldIndex});
        keepField(fieldIndex) = ...
            max(size(values)) == G.cells.num;
    end

    cellDataFields = fields(keepField);
    cellValueFormat = repmat( ...
        spacedValueFormat, 1, G.cells.num - 1);

    filename = sprintf('%s.pvd', name);
    pvdFileID = fopen(filename, 'w');

    assert(pvdFileID ~= -1, ...
        sprintf('Could not open file "%s" for writing.', filename));

    fprintf(pvdFileID, '<?xml version="1.0"?>\n');
    fprintf(pvdFileID, '<VTKFile type="Collection"\n');
    fprintf(pvdFileID, '\t\t\t\t version="0.1"\n');
    fprintf(pvdFileID, ...
        '\t\t\t\t byte_order="LittleEndian">\n');
    fprintf(pvdFileID, '\t<Collection>\n');

    hasTime = isfield(initialState, 'time');

    for stateIndex = 1 : numel(stateList)
        state = stateList{stateIndex};
        filename = sprintf('%s-%.05d.vtu', name, stateIndex);
        vtuFileID = fopen(filename, 'w');

        assert(vtuFileID ~= -1, ...
            sprintf('Could not open file "%s" for writing.', filename));

        fprintf(vtuFileID, '<?xml version="1.0"?>\n');
        fprintf(vtuFileID, [ ...
            '<VTKFile type="UnstructuredGrid" version="0.1" ' ...
            'byte_order="LittleEndian">\n']);
        fprintf(vtuFileID, '\t<UnstructuredGrid>\n');
        fprintf(vtuFileID, [ ...
            '\t\t<Piece NumberOfCells="%d" ' ...
            'NumberOfPoints="%d">\n'], ...
            G.cells.num, G.nodes.num);
        fprintf(vtuFileID, ...
            '\t\t\t<CellData Scalars="Variables">\n');

        for fieldIndex = 1 : numel(cellDataFields)
            fieldName = cellDataFields{fieldIndex};
            values = state.(fieldName);
            numberOfComponents = min(size(values));

            if numberOfComponents > 1
                for component = 1 : numberOfComponents
                    arrayName = sprintf( ...
                        '%s-%d', fieldName, component);
                    writeCellDataArray( ...
                        vtuFileID, arrayName, values, component, ...
                        cellValueFormat, valueFormat);
                end
            else
                writeCellDataArray( ...
                    vtuFileID, fieldName, values, 1, ...
                    cellValueFormat, valueFormat);
            end
        end

        fprintf(vtuFileID, '\t\t\t</CellData>\n');
        fprintf(vtuFileID, '%s', geometryText);
        fclose(vtuFileID);

        if hasTime
            timestep = state.time;
        else
            timestep = stateIndex - 1;
        end

        fprintf(pvdFileID, ...
            '\t\t<DataSet timestep="%d" file="%s"/>\n', ...
            timestep, filename);
    end

    fprintf(pvdFileID, '\t</Collection>\n');
    fprintf(pvdFileID, '</VTKFile>\n');
    fclose(pvdFileID);
end

function writeCellDataArray(fileID, arrayName, values, component, ...
                            cellValueFormat, valueFormat)
    fprintf(fileID, [ ...
        '\t\t\t\t<DataArray type="Float32" Name="%s" ' ...
        'NumberOfComponents="1" format="ascii">\n'], arrayName);
    fprintf(fileID, '\t\t\t\t\t');
    fprintf(fileID, cellValueFormat, ...
        values(1 : end - 1, component));
    fprintf(fileID, [valueFormat '\n'], ...
        values(end, component));
    fprintf(fileID, '\t\t\t\t</DataArray>\n');
end

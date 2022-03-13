classdef GPSNet
    %GPSNet - Model class implementing GPSNet type of network model
    %   Detailed explanation goes here

    properties
        graph
        W
        W_in
        model
        state0
        type
        nc
    end

    methods
        function obj = GPSNet(modelTrue, network, W_in, varargin)
            % Creates a GPSNet
            %
            % SYNOPSIS:
            %     obj =  GPSNet(model, network, W);
            %     obj =  GPSNet(model, network, W, nc);
            %
            % DESCRIPTION:
            % The function creates a GPSNet model starting from a model class that
            % defines the flow model and a graph that describes the interwell network.
            %
            % REQUIRED PARAMETERS:
            %   model       - An instance of the GenericBlackOilModel model class that
            %                 provides relevant description of fluid properties etc
            %                 from the flow model to be matched.
            %
            %   network     - graph containing the network elements nodes and edges,
            %                 typically found in the 'network' subfield of a Network
            %                 model.
            %
            %   W           - Well structure corresponding for the full-physics model
            %
            %   nc          - Number of cells in each conection.
            %                 Optional parameter. Default value: 10
            %
            %
            % OPTIONAL PARAMETERS:
            %
            %   'p0'        - Constant initial pressure. Default: 100 bar
            %
            %   'S0'        - Constant initial data. Must have the same
            %                 number of components as the fluid model.
            %                 Default: S_o=1, all other phases zero.
            %
            %   'fluid'     - Struct containing custom fluid model to be
            %                 used instead of the fine-scale fluid model
            %
            %   'scaling'   - cell array containing arguments used for
            %                 permeability scaling of custom fluid model
            %
            %   'Verbose'   - Indicate whether extra output is to be printed, such as
            %                 reports on how the faces and cell indices have been
            %                 assigned to each connection, and so on.
            % RETURNS:
            %   obj.model   - Model consisting of 1D flow paths packed into a 2D grid,
            %                 which thus has zero transmisibility in the vertical
            %                 direction and a number of non-neighboring connection
            %
            %   obj.W       - Converted  W for the new grid/model
            %
            % SEE ALSO:
            % `graph`, `Network`

            if mod(nargin,2)
                nc = 10;
            else
                nc = varargin{1};
                varargin = varargin(2:end);
            end
            opt = struct('Verbose',  mrstVerbose(),   ...
                         'p0',       100*barsa,       ...
                         'S0',       [],              ...
                         'fluid',    modelTrue.fluid, ...
                         'scaling',  []);
            opt = merge_options(opt, varargin{:});

            % Get the network graph and extract properties
            graph = network.network;
            numEdges = numedges(graph);
            numNodes = numnodes(graph);
            nodes    = graph.Nodes;
            
            % We subgrid each flow path into ten uniform cells and map the
            % resulting network onto a rectangular Cartesian grid having
            % the same number of rows as the number of flow paths. The grid
            % is set to have an aspect ratio of [5 1 1] and a volum that
            % matches the bulk volume of the reservoir. The constant
            % petrophysical properties are dummy values primarily used to
            % compute an initial guess for matching parameters that have
            % not be set by the network type.
            L = nthroot(sum(modelTrue.operators.pv./modelTrue.rock.poro)*25,3);
            G = cartGrid([nc, 1, numEdges], [L, L/5 ,L/5]*meter^3);
            G = computeGeometry(G);

            % Fluid model with endpoint scaling of relative permeabilities.
            % Unless a fluid model is supplied as an optional input
            % parameter, the fluid model is copied from the model we seek
            % to match.
            poro = mean(modelTrue.rock.poro);
            perm = mean(modelTrue.rock.perm(:,1));
            model = GenericBlackOilModel(G, makeRock(G, perm, poro), ...
                opt.fluid,'gas',modelTrue.gas);
            if ~isempty(opt.scaling)
                model = imposeRelpermScaling(model, opt.scaling{:});
            else
                pts   = opt.fluid.krPts;
                if isfield(pts,'o')
                    model = imposeRelpermScaling(model, ...
                        'SWL', pts.w(1,1), 'SWCR',  pts.w(1,2), ...
                        'SWU', pts.w(1,3), 'SOWCR', pts.o(1,2), ...
                        'KRW', pts.w(1,4), 'KRO',   pts.o(1,4));
                elseif isfield(pts,'ow')
                    model = imposeRelpermScaling(model, ...
                        'SWL', pts.w(1,1),  'SWCR',  pts.w(1,2),  ...
                        'SWU', pts.w(1,3),  'SOWCR', pts.ow(1,2), ...
                        'KRW', pts.ow(1,4), 'KRO',   pts.ow(1,4));
                end
            end
            model.OutputStateFunctions = {};
 
            %   Detailed explanation goes here
            obj.model = model;
            N  = model.operators.N;
            T  = model.operators.T;
            pv = model.operators.pv;
            
            % Set all transmissibility to zero except in the x-direction.
            % In a Cartesian grid, the faces are stored first in
            % x-direction, then in y and z. The T-vector only contains
            % values for the internal faces and thus only the
            % (nc-1)*numEdges values will be nonzero
            T((nc-1)*numEdges+1 :end) = 0;
            nodes.faceIx =zeros(numNodes, 1);
            
            wCells = nodes.cells;
            nf =  nc-1;
            en = graph.Edges.EndNodes;
            for i =  1:numEdges
                
                % Saving cell indices
                dispif(opt.Verbose,'Cell indices (edge %i): %i-%i \n',i,en(i,1),en(i,2));
                cellL = 1+(i-1)*nc;
                cellR = nc+(i-1)*nc;
                graph.Edges.cellIx{i} = cellL+1 : cellR-1;
                
                % Saving internal faces indices
                dispif(opt.Verbose,'in-Faces indices (edge %i): %i-%i \n',i,en(i,1),en(i,2));
                
                nodeL = en(i,1);
                nodeR = en(i,2);
                
                % Check if the first node was already defined and thus is
                % already connected to one or more rows in the 2D grid that
                % defines edges. If not, we assign this node to the first
                % cell of the next available row in the grid. If the node
                % is defined, we introduce a non-neighboring connection
                % between the cell it is defined in and the second cell in
                % the next available cell row.
                if nodes.faceIx(nodeL)== 0
                    
                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n',...
                        cellL ,W_in(nodes.well(nodeL)).name);
                    
                    firstIx =  1+(i-1)*nf;
                    nodes.faceIx(nodeL) = firstIx;
                    wCells{nodes.well(nodeL)} = [wCells{nodes.well(nodeL)},cellL];
                    nodes.cells{nodeL}= [nodes.cells{nodeL},cellL];

                else
                    %Create new connectivity by adding artifical face
                    N = [N; nodes.cells{nodeL}(end), cellL+1];             %#ok<AGROW>
                    T = [T;T(1)];                                          %#ok<AGROW>
                    
                    firstIx = length(T); % This the index of a new face
                    % between cells nodes.cells(nodeL) and cellL+1
                end
                
                % Internal faces number
                internIx = (2+(i-1)*nf:(nc-2)+(i-1)*nf) ;
                
                % Same as for the first node: check if node is assigned to
                % a cell. If not, assign it to the last cell of the current
                % cell row. If it exists, introduce a non-neighboring
                % connection between the second last cell in the current
                % row and the cell in which the node is defined already.
                if nodes.faceIx(nodeR)==0
                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n', ...
                        cellR ,W_in(nodes.well(nodeR)).name);
                    lastIx  =  (nc-1)+(i-1)*nf;
                    nodes.faceIx(nodeR) = lastIx;
                    wCells{nodes.well(nodeR)} = [wCells{nodes.well(nodeR)},cellR];
                    nodes.cells{nodeR}= [nodes.cells{nodeR},cellR];
                else
                    %Create new connectivity by adding a new artifical face
                    N = [N; cellR-1, nodes.cells{nodeR}(end)];             %#ok<AGROW>
                    T = [T;T(1)];                                          %#ok<AGROW>
                    
                    lastIx = length(T); % This the index of a new face
                    % between cells cellR-1 and nodes.cells(nodeR)
                end
                dispif(opt.Verbose,'\n');
                graph.Edges.faceIx{i} = [firstIx,internIx,lastIx];
            end
            
            % Creating the new well structure.
            W = [];
            for iw=1:numel(W_in)
                W = addWell(W, model.G, model.rock,  wCells{iw}, ...
                    'Type',   W_in(iw).type,  ...
                    'Val',    W_in(iw).val,   ...
                    'Radius', W_in(iw).r(1:numel(wCells{iw})),...
                    'Name',   W_in(iw).name,  ...
                    'Comp_i', W_in(iw).compi, ...
                    'sign',   W_in(iw).sign,  ...
                    'lims',   W_in(iw).lims,  ...
                    'status', W_in(iw).status);
            end
            
            % Initialize the model
            if isempty(opt.S0)
                opt.S0 = zeros(1,model.getNumberOfPhases());
                opt.S0(model.getPhaseIndex('O')) = 1;
            else
                assert(numel(opt.S0)==model.getNumberOfPhases(),...
                    'Initial data does not match number of phases');
            end
            state0 = initState(model.G, W, opt.p0, opt.S0);
            
            % Assign data objects
            graph.Nodes = nodes;
            obj.graph   =  graph;
            obj.W       =  W;
            obj.W_in    =  W_in;
            obj.state0  = state0;
            obj.type    = network.type;
            obj.nc      = nc;
            
            obj.model.operators = ...
                setupOperatorsTPFA(obj.model.G,  obj.model.rock, ...
                'neighbors',N,'trans',T,'porv',pv);
            obj.model = obj.model.validateModel();
            obj.model.toleranceCNV = 1e-6;
        end
               
        function plotGrid(obj,data,varargin)
            %Plot the grid used in the GPSNet with or without data.
            %
            % SYNOPSIS:
            %       plotGrid([])
            %       plotGrid(data, 'pn1', pv1, ...)
            %
            % PARAMETERS:
            %   data - Scalar cell data with which to colour the grid.  One
            %          scalar, indexed colour value for each cell in the
            %          grid. If empty, the routine will plot each row in
            %          the grid with a unique color.
            %% Plotting the data driven model
            G = obj.model.G;
            [I,~,K] = gridLogicalIndices(G);
            
            if isempty(data)
                plotCellData(G,K);
                plotGrid(G,I==1 | I==max(I), 'FaceColor',[.9 .9 .9]);
                plotGrid(G,vertcat(obj.W.cells),'FaceColor',[.7 .7 .7]);
            else
                plotCellData(G,data, 'EdgeAlpha', 0.1)
                plotGrid(G,I==1 | I==max(I), 'FaceColor','none','linewidth',1);
                plotGrid(G,vertcat(obj.W.cells),'FaceColor','none','linewidth',2);
            end
            
            for i =  1: numel(obj.graph.Nodes.well)
                cellNo =obj.graph.Nodes.cells(i);
                cellNo = cellNo{1};
                if ~isempty(cellNo)
                    XData = G.cells.centroids(cellNo,1);
                    YData = G.cells.centroids(cellNo,2);
                    ZData = G.cells.centroids(cellNo,3);
                    
                    if (abs(XData-G.cells.centroids(1,1)) < ...
                            abs(XData-G.cells.centroids(end,1)))
                        text(XData-70,YData,ZData,...
                            obj.graph.Nodes.name{i});
                    else
                        text(XData+20,YData,ZData,...
                            obj.graph.Nodes.name{i});
                    end
                end
            end
            view(0,0), axis off
        end
        
        function [edge,subset] = getMapping(obj,type)
            %Get mapping between graph edges and cell/faces in the grid.
            %
            % SYNOPSIS:
            %   [edges,subset] = getMapping(type)
            %       plotGrid(data, 'pn1', pv1, ...)
            %
            % PARAMETERS:
            %   type - Mapping type, either 'cells' or 'faces'
            %
            % OUTPUT:
            %   edges  - list of (repeated) edge numbers
            %   subset - list of cell/face indices for each edge
            
            switch type
                case 'cells'
                    indices = obj.graph.Edges.cellIx;
                case 'faces'
                    indices = obj.graph.Edges.faceIx;
                otherwise
                    error('Unknown mapping type');
            end
            A = cellfun(@(C)numel(C),indices);
            
            sizeAllIndices = sum(A);
            
            edge = zeros(sizeAllIndices,1);
            subset  = zeros(sizeAllIndices,1);
            
            reel = 1;
            for i=1:numel(indices)
                n = numel(indices{i});
                subset(reel:n+reel-1) = indices{i};
                edge(reel:n+reel-1)= i;
                reel = reel+ n;
            end
        end
        
        function pvec = getScaledParameterVector(obj, setup, params, varargin)
            %Get parameters scaled to the unit interval
            %
            % SYNOPSIS:
            %   pvec = getScaledParameterVector(setup, params)
            %   pvec = getScaledParameterVector(setup, params, 'pn', pv)
            %
            % DESCRIPTION:
            %   The function does the same as getScaledParameterVector from
            %   tha optimization module. The two exceptsion are in the case
            %   of a network created from flow diagnostics, in which the
            %   pore volumes and transmissibilities are overwritten by
            %   quantities stored in the network graph, and for well
            %   transmissibilities, for which the user can specify a
            %   multiplier value.
            %
            % PARAMETERS:
            %   setup  - simulation setup structure
            %   params - cell array of ModelParameter instances
            %
            % OPTIONAL PARAMETERS:
            %   'connscale' - scales well connectivities. Either a scalar
            %              value or one value per connection in params
            %   'randSw' - random initial guess for the water saturation
            %
            % OUTPUT:
            %   pvec - array containing scaled parameters
            
            opt = struct('connscale',  ones(numel(params),1), ...
                         'randSw',     false);
            opt = merge_options(opt, varargin{:});

            if numel(opt.connscale)==1
                opt.connscale = opt.connscale*ones(numel(params),1);
            else
                assert(numel(opt.connscale)==numel(params), ...
                    'Number of scaling parameters not equal number of connections');
            end
            values = applyFunction(@(p)p.getParameterValue(setup), params);
            
            isdiagnost = any(strcmp(obj.type,{'fd_preprocessor','fd_postprocessor'}));
            for k=numel(params):-1:1
                switch params{k}.name
                    case 'porevolume'
                        if isdiagnost
                            values{k} = obj.graph.Edges.pv/obj.nc;
                        end
                    case 'transmissibility'
                        if isdiagnost
                            values{k} = obj.graph.Edges.T;
                        end
                    case 'conntrans'
                        values{k} = opt.connscale(k)*values{k};
                    case 'sw'
                        if opt.randSw
                            values{k} =  rand(size(values{k}));
                        end
                end
                u{k} = params{k}.scale(values{k});
            end
            pvec = vertcat(u{:});   
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

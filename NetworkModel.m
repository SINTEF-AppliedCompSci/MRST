classdef NetworkModel
    %NetworkModel - Model class implementing GPSNet type of network model
    %   Detailed explanation goes here

    properties
        graph
        W
        W_in
        model
        state0
    end

    methods
        function obj = NetworkModel(modelTrue, network, W_in, varargin)
            % Creates a NetworkModel
            %
            % SYNOPSIS:
            %     obj =  NetworkModel(model, network, W);
            %     obj =  NetworkModel(model, network, W, nc);
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
            opt = struct('Verbose',  mrstVerbose());
            opt = merge_options(opt, varargin{:});

            % Get the network graph and extract properties
            graph = network.network;
            numEdges = numedges(graph);
            
            % We subgrid each flow path into ten uniform cells and map the
            % resulting network onto a rectangular Cartesian grid having
            % the same number of rows as the number of flow paths. The grid
            % is set to have an aspect ratio of [5 1 1] and a volum that
            % matches the bulk volume of the reservoir. The constant
            % petrophysical properties are dummy values primarily used to
            % compute an initial guess for matching parameters that have
            % not be set by the network type.
            L = nthroot(sum(modelTrue.operators.pv./modelTrue.rock.poro)*25,3)  ;
            G = cartGrid([nc, 1, numEdges], [L, L/5 ,L/5]*meter^3);
            G = computeGeometry(G);

            % Fluid model with endpoint scaling of relative permeabilities.
            % The fluid model is copied from the model we seek to match.            
            model = GenericBlackOilModel(G, makeRock(G, 1*darcy, 0.2), ...
                modelTrue.fluid,'gas',modelTrue.gas);
            pts   = modelTrue.fluid.krPts;
            scaling = {'SWL',   pts.w(1), 'SWCR', pts.w(2), 'SWU', pts.w(3), ...
                'SOWCR', pts.o(2), 'KRW',  pts.w(4), 'KRO', pts.o(4)};
            model = imposeRelpermScaling(model, scaling{:});
            model.OutputStateFunctions = {};
 
            %   Detailed explanation goes here
            obj.model = model;
            N  = model.operators.N;
            T  = model.operators.T;
            pv = model.operators.pv;
            
            number_of_vertical_faces = (nc-1)*numEdges; %Internal faces
            T(number_of_vertical_faces+1 :end) = 0;
            graph.Nodes.Face_Indx =0*(1:numnodes(graph))';%TODO: check this, include general counts for nodes
            
            Well_cells = graph.Nodes.Well_cells; %TODO Initialize better
            nf =  nc -1;
            en = graph.Edges.EndNodes;
            for i_ed =  1:numEdges
                
                % Saving cell indices
                dispif(opt.Verbose,'Cell indices Edge %i, %i-%i \n',i_ed,en(i_ed,1),en(i_ed,2));
                cellL   = 1+(i_ed-1)*nc;
                cellR = nc+(i_ed-1)*nc;
                graph.Edges.Cell_Indices{i_ed} = cellL+1 : cellR-1;
                
                % Saving internal faces indices
                dispif(opt.Verbose,'in-Faces indices Edge %i, %i-%i \n',i_ed,...
                    en(i_ed,1),en(i_ed,2));
                
                nodeL  = en(i_ed,1);
                nodeR = en(i_ed,2);
                
                % Check if the first node was already defined and has its
                % face number index stored
                if graph.Nodes.Face_Indx(nodeL)== 0
                    
                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n',...
                        cellL ,W_in(graph.Nodes.Well(nodeL)).name);
                    
                    First_Index =  1+(i_ed-1)*nf;
                    graph.Nodes.Face_Indx(nodeL) = First_Index;
                    Well_cells{graph.Nodes.Well(nodeL)} = [Well_cells{graph.Nodes.Well(nodeL)},cellL];
                    graph.Nodes.Well_cells{nodeL}= [graph.Nodes.Well_cells{nodeL},cellL];
                else
                    %Create New connectivity
                    N = [N;...
                        graph.Nodes.Well_cells{nodeL}(end), cellL+1]; % New artificial face
                    T = [T;T(1)];
                    
                    First_Index = length(T); % This the index of a new face
                    % between cells
                    % graph.Nodes.Well_cells(Node_1)
                    % and  Cell_1+1
                end
                
                % Internal faces number
                Internal = (2+(i_ed-1)*nf:(nc-2)+(i_ed-1)*nf) ;
                
                % Check if the last node was already defined and has its
                % face number index stored
                if graph.Nodes.Face_Indx(nodeR) ==0
                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n', ...
                        cellR ,W_in(graph.Nodes.Well(nodeR)).name);
                    Last_Index  =  (nc-1)+(i_ed-1)*nf;
                    graph.Nodes.Face_Indx(nodeR) = Last_Index;
                    Well_cells{graph.Nodes.Well(nodeR)} = [Well_cells{graph.Nodes.Well(nodeR)},cellR];
                    graph.Nodes.Well_cells{nodeR}= [graph.Nodes.Well_cells{nodeR},cellR];
                else
                    %Create new conectivity
                    N = [N;...
                        cellR-1, graph.Nodes.Well_cells{nodeR}(end)]; % New artificial face
                    T = [T;T(1)];
                    
                    Last_Index = length(T); % This the index of a new face
                    % between cells
                    % graph.Nodes.Well_cells(Node_1)
                    % and  Cell_1+1
                end
                dispif(opt.Verbose,'\n');
                graph.Edges.Face_Indices{i_ed} = [First_Index,Internal,Last_Index];
            end
            
            % Creating the new well structure.
            W = [];
            for iw = 1 : numel(W_in)
                W = addWell(W, model.G, model.rock,  Well_cells{iw}, ...
                    'Type', W_in(iw).type,...
                    'Val', W_in(iw).val, ...
                    'Radius',W_in(iw).r(1:numel(Well_cells{iw})),...
                    'Name',W_in(iw).name,...
                    'Comp_i',W_in(iw).compi,...
                    'sign', W_in(iw).sign);
            end
            
            % Initialize the model
            S0 = zeros(1,model.getNumberOfPhases());
            S0(model.getPhaseIndex('O')) = 1;
            state0 = initState(model.G, W, 100*barsa, S0);
            
            obj.graph =  graph;
            obj.W     =  W;
            obj.W_in  =  W_in;
            obj.state0 = state0;
            
            obj.model.operators = ...
                setupOperatorsTPFA(obj.model.G,  obj.model.rock, ...
                'neighbors',N,'trans',T,'porv',pv);
            obj.model = obj.model.validateModel();
            obj.model.toleranceCNV = 1e-6;
        end
        
        
        function plotNetwork(obj,G,data,varargin)
            plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
            
            plotWell(G,obj.W_in, 'color', 'k','fontsize',0);
            
            if isempty(data)
                linewidth = 5;
            else
                assert (numel(data) == numedges(obj.graph), ...
                    'The DATA should have one value for each grid cell in output.');
                linewidth = 10*data/max(data);
            end
            
            hold on, pg =  plot(obj.graph,...
                'XData',obj.graph.Nodes.XData,...
                'YData',obj.graph.Nodes.YData,...
                'ZData',obj.graph.Nodes.ZData,...
                'LineWidth',linewidth);
            labelnode(pg,1:numnodes(obj.graph),obj.graph.Nodes.Well_name);
            hold off;
            pg.NodeFontSize= 10;
            axis off ;
        end
        
        
        function plotGrid(obj,data,varargin)
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
            
            for i =  1: numel(obj.graph.Nodes.Well)
                cell_number =obj.graph.Nodes.Well_cells(i);
                cell_number = cell_number{1};
                if ~isempty(cell_number)
                    XData = G.cells.centroids(cell_number,1);
                    YData = G.cells.centroids(cell_number,2);
                    ZData = G.cells.centroids(cell_number,3);
                    
                    if (abs(XData-G.cells.centroids(1,1)) < ...
                            abs(XData-G.cells.centroids(end,1)))
                        text(XData-70,YData,ZData,...
                            obj.graph.Nodes.Well_name{i});
                    else
                        text(XData+20,YData,ZData,...
                            obj.graph.Nodes.Well_name{i});
                    end
                end
            end
            view(0,0), axis off
        end
        
        function [edge,subset] = getMapping(obj,type)
            
            switch type
                case 'cells'
                    indices = obj.graph.Edges.Cell_Indices;
                case 'faces'
                    indices = obj.graph.Edges.Face_Indices;
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
    end
end

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

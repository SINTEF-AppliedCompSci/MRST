classdef NetworkModel
    %NetworkModel - Model class implementing GPSNet type of network model
    %   Detailed explanation goes here

    properties
        Graph
        W
        W_in
        model
    end

    methods
        function obj = NetworkModel(model,nc,Graph,W_in,varargin)
            % Creates a NetworkModel
            %
            % SYNOPSIS:
            %     obj =  NetworkModel(model, nc, network, W);
            %
            %
            % DESCRIPTION:
            % The function creates a network model starting from a model class that
            % defines the flow model and graph that describes the interwell network.
            %
            % REQUIRED PARAMETERS:
            %   model       - An instance of the GenericBlackOilModel model class that
            %                 provides relevant description of fluid properties etc
            %                 from the flow model to be matched.
            %
            %   nc          - Number of cell on each conection.
            %                 TODO: Make it optional parameter
            %
            %   network     - Graph containing the network elements nodes and edges,
            %                 typically found in the 'network' subfield of a Network
            %                 model.
            %
            %  W            - Well structure corresponding with the full phisics model
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
            %   obj.indexs  - Structure containing the cell indices and faces indices
            %                 for each connection/edge in Graph. This are using as
            %                 lumping parameters inside ModelParameter Class
            %
            % SEE ALSO:
            % `graph`, `Network`

            opt = struct('Verbose',  mrstVerbose());
            opt = merge_options(opt, varargin{:});


            %TODO: get numedges and numlevels and numcros and define and model
            % using fluid and rock


            %   Detailed explanation goes here
            obj.model = model;
            N = model.operators.N;
            T = model.operators.T;
            pv= model.operators.pv;

            number_of_vertical_faces = (nc-1)*size(Graph.Edges.EndNodes,1); %Internal faces
            T(number_of_vertical_faces+1 :end) = 0;
            Graph.Nodes.Face_Indx =0*(1:numnodes(Graph))';%TODO: check this, include general counts for nodes

            Well_cells = Graph.Nodes.Well_cells; %TODO Initialize better
            nf =  nc -1;
            en = Graph.Edges.EndNodes;
            for i_ed =  1:numedges(Graph)

                %% Savinf cell indexs
                dispif(opt.Verbose,'Cell indices Edge %i, %i-%i \n',i_ed,en(i_ed,1),en(i_ed,2));
                Cell_1   = 1+(i_ed-1)*nc;
                Cell_nc = nc+(i_ed-1)*nc;
                Graph.Edges.Cell_Indices{i_ed} = Cell_1+1 : Cell_nc-1;
                %% Saving internal faces indices
                dispif(opt.Verbose,'in-Faces indices Edge %i, %i-%i \n',i_ed,...
                    en(i_ed,1),en(i_ed,2));

                Node_1  = en(i_ed,1);
                Node_nc = en(i_ed,2);
                % Check if the firs node was already defined and has its
                % face number index stored
                if Graph.Nodes.Face_Indx(Node_1)== 0

                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n',...
                        Cell_1 ,W_in(Graph.Nodes.Well(Node_1)).name);

                    First_Index =  1+(i_ed-1)*nf;
                    Graph.Nodes.Face_Indx(Node_1) = First_Index;
                    Well_cells{Graph.Nodes.Well(Node_1)} = [Well_cells{Graph.Nodes.Well(Node_1)},Cell_1];
                    Graph.Nodes.Well_cells{Node_1}= [Graph.Nodes.Well_cells{Node_1},Cell_1];
                else
                    %Create New connectivity
                    N = [N;...
                        Graph.Nodes.Well_cells{Node_1}(end), Cell_1+1]; % New artificial face
                    T = [T;T(1)];

                    First_Index = length(T); % This the index of a new face
                    % between cells
                    % Graph.Nodes.Well_cells(Node_1)
                    % and  Cell_1+1
                end

                % Internal faces number
                Internal    =  (2+(i_ed-1)*nf:(nc-2)+(i_ed-1)*nf) ;

                % Check if the last node was already defined and has its
                % face number index stored
                if Graph.Nodes.Face_Indx(Node_nc) ==0
                    dispif(opt.Verbose,2,'Cell number %i, will have Well %s \n', ...
                        Cell_nc ,W_in(Graph.Nodes.Well(Node_nc)).name);
                    Last_Index  =  (nc-1)+(i_ed-1)*nf;
                    Graph.Nodes.Face_Indx(Node_nc) = Last_Index;
                    Well_cells{Graph.Nodes.Well(Node_nc)} = [Well_cells{Graph.Nodes.Well(Node_nc)},Cell_nc];
                    Graph.Nodes.Well_cells{Node_nc}= [Graph.Nodes.Well_cells{Node_nc},Cell_nc];
                else
                    %Create new conectivity
                    N = [N;...
                        Cell_nc-1, Graph.Nodes.Well_cells{Node_nc}(end)]; % New artificial face
                    T = [T;T(1)];

                    Last_Index = length(T); % This the index of a new face
                    % between cells
                    % Graph.Nodes.Well_cells(Node_1)
                    % and  Cell_1+1
                end
                dispif(opt.Verbose,'\n');
                Graph.Edges.Face_Indices{i_ed} = [First_Index,Internal,Last_Index];
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

            obj.Graph =  Graph;
            obj.W     =  W;
            obj.W_in  =  W_in;

            obj.model.operators =  setupOperatorsTPFA( obj.model.G,  obj.model.rock,...
                'neighbors',N,...
                'trans',T,...
                'porv',pv);
        end


        function plotNetwork(obj,Grid,data,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            plotGrid(Grid, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);

            plotWell(Grid,obj.W_in, 'color', 'k','fontsize',0);

            if isempty(data)
                linewidth = 5;
            else
                assert (numel(data) == numedges(obj.Graph), ...
                    'The DATA should have one value for each grid cell in output.');
                linewidth = 10*data/max(data);
            end

            hold on, pg =  plot(obj.Graph,...
                'XData',obj.Graph.Nodes.XData,...
                'YData',obj.Graph.Nodes.YData,...
                'ZData',obj.Graph.Nodes.ZData,...
                'LineWidth',linewidth);
            labelnode(pg,1:numnodes(obj.Graph),obj.Graph.Nodes.Well_name);
            hold off;
            pg.NodeFontSize= 10;
            axis off ;
            title("Initial Transmisibility");
        end


        function plotGPSNetGrid(obj,data,varargin)
            %% Plotting the data driven model
            G = obj.model.G;
            [I,~,K] = gridLogicalIndices(G);

            if isempty(data)
                plotCellData(G,K);
                plotGrid(G,I==1 | I==max(I), 'FaceColor',[.9 .9 .9]);
                plotGrid(G,vertcat(obj.W.cells),'FaceColor',[.7 .7 .7]);
            else
                assert (numel(data) == G.cells.num, ...
                    'The DATA should have one value for each grid cell in output.');
                plotCellData(G,data, 'EdgeAlpha', 0.1)
                plotGrid(G,I==1 | I==max(I), 'FaceColor','none','linewidth',1);
                plotGrid(G,vertcat(obj.W.cells),'FaceColor','none','linewidth',2);
             end

            for i =  1: numel(obj.Graph.Nodes.Well)
                cell_number =obj.Graph.Nodes.Well_cells(i);
                cell_number = cell_number{1};
                if ~isempty(cell_number)
                    XData = G.cells.centroids(cell_number,1);
                    YData = G.cells.centroids(cell_number,2);
                    ZData = G.cells.centroids(cell_number,3);

                    if (abs(XData-G.cells.centroids(1,1)) < ...
                            abs(XData-G.cells.centroids(end,1)))
                        text(XData-70,YData,ZData,...
                            obj.Graph.Nodes.Well_name{i});
                    else
                        text(XData+20,YData,ZData,...
                            obj.Graph.Nodes.Well_name{i});
                    end
                end
            end
            view(0,0), axis off
        end
    end
end

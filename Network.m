classdef Network
    properties
        network
        G
        W
    end
    methods
        function obj = Network(W,G,varargin)
            % Creates a network between nodes defined by a set of wells
            %
            % SYNOPSIS:
            %     ntwrk =  Network(W, G, 'type', 'all_to_all');
            %
            %     ntwrk =  Network(W, G, 'type', 'injectors_to_producers',...
            %                            'injectors', 1:n_injectors,...
            %                            'producers', n_injectors + 1:n_producer);
            %
            %   Networks derived by flow diagnostics
            %
            %     ntwrk =  Network(W, G, 'type', 'fd_preprocessor',...
            %                            'problem', problem,...
            %                            'flow_filter', 1*stb/day);
            %
            %     ntwrk =  Network(W, G, 'type', 'fd_postprocessor',...
            %                            'problem', problem,...
            %                            'flow_filter', 1*stb/day,...
            %                            'state_number', 40);
            %
            % DESCRIPTION:
            %   Nodes are created from well cells given in the well structure W, and
            %   the edges of the graph are created depending on which type of network is
            %   required. These types are described briefly in the optimal parameter
            %   section.
            %
            % REQUIRED PARAMETERS:
            %   W  - Wells to be used as nodes in the network. Each perforated cell
            %        will be cosidered as a node.
            %
            %   G  - Grid structure with field cells.centroids defined
            %
            %
            % OPTIONAL PARAMETERS:
            %   'type' - Architecture/design of the network that specifies connections between nodes:
            %               'all_to_all':       Define edges to connect all
            %                                   nodes. Connections between nodes
            %                                   of the same well are not
            %                                   allowed.
            %
            %               'injectors_to_producers': Define edges only between
            %                                   nodes from injector well cells
            %                                   to producer well cells.
            %
            %               'fd_preprocessor':  Define edges only between
            %                                   connections discovered by the flow
            %                                   diagnostics preprosessor.
            %
            %               'fd_postprocessor': Define edges only between
            %                                   connections discovers by the flow
            %                                   diagnostics postprocessor.
            %
            %               'user_defined_edges': Defined edges by the user.
            %                                   This requires parameter 'edges'
            %
            %   'injectors'  - Indices for the injector wells inside the W structure
            %
            %   'producers'  - Indices for the producer wells inside the W structure
            %
            %   'edges'      - Array to indicate connections between well nodes.
            %
            %   'problem'    - A single packed problem from 'packSimulationProblem'
            %
            %   'flow_filter'- Flow value used to filter all conections that have a
            %                  lower flow allocation. This option is only used when
            %                  'type' is set 'fd_preprocessor' or 'fd_preprocessor'.
            %
            %   'state_number' - Indicates for time steps from which to read the states
            %                  passed on to the flow diagnostics solver. This option is
            %                  only used when 'type' is set 'fd_postprocessor'.
            %
            %   'Verbose'    - Indicates whether extra output is to be printed, such as
            %                  reports on how the faces and cell indices have been
            %                  assigned to each connection, and so on.
            % RETURNS:
            %   obj.network  - MATLAB Undirected Graph
            %
            % SEE ALSO:
            % `graph`, `NetworkModel`
            
            
            obj.G = G;
            obj.W = W;
            opt = struct('Verbose',      mrstVerbose(),...
                         'type',         'all_to_all', ...
                         'edges',        [], ...
                         'injectors',    [], ...
                         'producers',    [], ...
                         'problem',      [], ...
                         'flow_filter',  0,  ...
                         'state_number', []);
            opt = merge_options(opt, varargin{:});
            
            if ~isempty(opt.edges)
                opt.type = 'user_defined_edges';
            end
            
            if ~isempty(opt.injectors)&&~isempty(opt.producers)
                opt.type = 'injectors_to_producers';
            end
            
            if isempty(opt.state_number)&&~isempty(opt.problem)
                opt.state_number = numel(opt.problem.SimulatorSetup.schedule.step.val) ;
            end
            
            
            nW_cells  = arrayfun(@(x)numel(x.cells), W);
            num_nodes = sum(nW_cells);
            
            % Define everything for each node
            node = numel(vertcat(W.cells));
            for iw = numel(W):-1:1
                for iw_cell = numel(W(iw).cells):-1:1
                    Nodes(node,1)      = node;
                    Well(node,1)       = iw ;
                    SubWell(node,1)    = iw_cell;
                    Well_name{node,1}  = W(iw).name;
                    Well_cells{node,1} = [];
                    Type(node,1)       = 1;
                    cellNo             = W(iw).cells(iw_cell);
                    XData(node,1)      = obj.G.cells.centroids(cellNo,1);
                    YData(node,1)      = obj.G.cells.centroids(cellNo,2);
                    ZData(node,1)      = obj.G.cells.centroids(cellNo,3);
                    node =  node-1;
                end
            end
            
            % Create edges of the graph. Only rule is:
            % self-connections among wells are forbidden
            switch opt.type
                case 'all_to_all'
                    A = ones(num_nodes)- eye([num_nodes,num_nodes]);
                case 'user_defined_edges'
                    A = sparse([opt.edges(:,1); opt.edges(:,2)], ...
                               [opt.edges(:,2); opt.edges(:,1)], ...
                               ones(2*size(opt.edges,1),1));
                case 'injectors_to_producers'
                    A =  zeros(num_nodes);
                    inj  = any(Well==opt.injectors,2);
                    prod = any(Well==opt.producers,2);
                    A(inj,prod) = 1;
                    A(prod,inj) = 1;
                    %figure, spy(A);
                case {'fd_preprocessor','fd_postprocessor'}
                    % Compute flow diagnostics for the chosen state
                    require diagnostics
                    assert(all(nW_cells==1),['Flow diagnostics analysis to multiple connections',...
                        ' between wells is not yet supported.'])
                    if strcmp(opt.type,'fd_preprocessor')
                        state = [];
                        pressure_field = 'pressure';
                        ctrlNo = 1;
                    else
                        state = opt.problem.OutputHandlers.states{opt.state_number};
                        pressure_field = 'bhp';
                        ctrlNo = opt.problem.SimulatorSetup.schedule.step.control(opt.state_number);
                    end
                    [state, diagnostics] = computePressureAndDiagnostics(...
                        opt.problem.SimulatorSetup.model,...
                        'wells', opt.problem.SimulatorSetup.schedule.control(ctrlNo).W,...
                        'state',state, 'firstArrival', false);
                    
                    % Use diagnostics to define the graph
                    ix   = diagnostics.wellCommunication(:) > opt.flow_filter;
                    I    = diagnostics.WP.pairIx(ix,1);
                    P    = diagnostics.WP.pairIx(ix,2);
                    iWno = diagnostics.D.inj(I)';
                    pWno = diagnostics.D.prod(P)';
                    flux = diagnostics.wellCommunication(ix);
                    dP   = vertcat(state.wellSol(iWno).(pressure_field)) - ...
                        vertcat(state.wellSol(pWno).(pressure_field));
                    T    = flux./dP;
                    pv   = diagnostics.WP.vols(ix)';
                    A    = sparse([iWno; pWno],[pWno; iWno],ones(2*sum(ix(:)),1));
                otherwise
                    error('\nType of network: %s is not implemented\n', opt.type);
            end
                        
            obj.network = graph(A);
            obj.network.Nodes.Nodes     = Nodes;
            obj.network.Nodes.Well      = Well;
            obj.network.Nodes.SubWell   = SubWell;
            obj.network.Nodes.Type      = Type;
            obj.network.Nodes.Well_name = Well_name;
            obj.network.Nodes.Well_cells = Well_cells;
            obj.network.Nodes.XData     = XData;
            obj.network.Nodes.YData     = YData;
            obj.network.Nodes.ZData     = ZData;
            if any(strcmp(opt.type,{'fd_preprocessor','fd_postprocessor'}))
                obj.network.Edges.Transmissibility = T;
                obj.network.Edges.PoreVolume = pv;
            end
        end
        
        function plotNetwork(obj,varargin)
            opt = struct('FaceColor', 'none', 'EdgeAlpha', 0.1);
            opt = merge_options(opt, varargin{:});
            
            if size(obj.network.Edges,2)==4
                TT = obj.network.Edges.Transmissibility;
                pv = obj.network.Edges.PoreVolume;
                lineWidth = [10*TT/max(TT), 10*pv/max(pv)];
                names = {'Transmissibility', 'Pore volume'};
                n = 2;
            else
                lineWidth = 2*ones(numedges(obj.network),1);
                names = {''};
                n = 1;
            end
            for i=1:n
                subplot(1,n,i,'replace')
                plotGrid(obj.G, 'FaceColor',opt.FaceColor,...
                    'EdgeAlpha',opt.EdgeAlpha); view(2);
            
                hold on, pg =  plot(obj.network,...
                    'XData',obj.network.Nodes.XData,...
                    'YData',obj.network.Nodes.YData,...
                    'ZData',obj.network.Nodes.ZData,...
                    'LineWidth',lineWidth(:,i));
                labelnode(pg,obj.network.Nodes.Well,obj.network.Nodes.Well_name);
                title(names{i})
                hold off; axis off
            end
        end
        
    end
end
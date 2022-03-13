classdef Network
    properties
        network
        G
        W
        type
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
            
            
            nWcells  = arrayfun(@(x)numel(x.cells), W);
            numNodes = sum(nWcells);
            
            % Define everything for each node
            node = numel(vertcat(W.cells));
            for iw = numel(W):-1:1
                for i = numel(W(iw).cells):-1:1
                    nodes(node,1)   = node;
                    well(node,1)    = iw ;
                    subwell(node,1) = i;
                    name{node,1}    = W(iw).name;
                    cells{node,1}   = [];
                    type(node,1)    = 1;
                    cellNo          = W(iw).cells(i);
                    XData(node,1)   = obj.G.cells.centroids(cellNo,1);
                    YData(node,1)   = obj.G.cells.centroids(cellNo,2);
                    ZData(node,1)   = obj.G.cells.centroids(cellNo,3);
                    node =  node-1;
                end
            end
            
            % Create edges of the graph. Only rule is:
            % self-connections among wells are forbidden
            switch opt.type
                case 'all_to_all'
                    obj.network = graph(ones(numNodes) - eye([numNodes,numNodes]));
                    
                case 'user_defined_edges'
                    % Check that all nodes have at least one edge and that
                    % no edges involve non-existing nodes
                    a = accumarray(opt.edges(:),1);
                    assert(numel(a)<=numNodes,'Edges refer to non-existing node(s)');
                    assert((numel(a)==numNodes) && all(a>0),...
                        'Each node must have at least one edge');
                    obj.network = graph(opt.edges(:,1),opt.edges(:,2));
                    
                case 'injectors_to_producers'
                    ni = numel(opt.injectors);
                    np = numel(opt.producers);
                    obj.network = graph(...
                        rldecode(opt.injectors',np*ones(ni,1)), ...
                        repmat(opt.producers',ni, 1));
                    
                case {'fd_preprocessor','fd_postprocessor'}
                    % Compute flow diagnostics for the chosen state
                    require diagnostics
                    assert(all(nWcells==1),...
                        ['Flow diagnostics analysis to multiple ',...
                         'connections between wells is not yet supported.'])
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
                        'state', state, 'firstArrival', false);
                    
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
                    obj.network = graph(iWno, pWno);
                    
                otherwise
                    error('\nType of network: %s is not implemented\n', opt.type);
            end
            obj.type = opt.type;
                        
            obj.network.Nodes.nodes   = nodes;
            obj.network.Nodes.well    = well;
            obj.network.Nodes.subwell = subwell;
            obj.network.Nodes.type    = type;
            obj.network.Nodes.name    = name;
            obj.network.Nodes.cells   = cells;
            obj.network.Nodes.XData   = XData;
            obj.network.Nodes.YData   = YData;
            obj.network.Nodes.ZData   = ZData;
            if any(strcmp(opt.type,{'fd_preprocessor','fd_postprocessor'}))
                obj.network.Edges.T  = T;
                obj.network.Edges.pv = pv;
            end
        end
        
        function omap=plotNetwork(obj,varargin)
            %Plot the grid used in the GPSNet with or without data.
            %
            % SYNOPSIS:
            %          plotGrid()
            %   cmap = plotGrid(type, 'pn1', pv1, ...)
            %
            % PARAMETERS:
            %   type is a string that determines the layout of the plot
            %     'default' - plot in 3D space with uniform edge width
            %     'circle'  - plot the network with all nodes on a circle
            %     'transmissibility' - 3D space with edge width scaled
            %                  according to the transmissibilities
            %                  associated with each edge
            %     'porevolume' - 3D space with edge width scaled according
            %                  to pore volumes associated with each edge
            %   The last two layouts presume that the network has been computed
            %   using flow diagnostics
            %
            % OPTIONAL PARAMETERS:
            %  'Colors'   - use a unique color for each edge. The resulting
            %               colormap can be returned in the 'cmap'
            %  'onGrid'   - plot the underlying 3D grid
            %  'maxWidth' - the maximum line width for the edges
            %  'data'     - array with one element per edge, used to scale
            %               the line width of the edges
            %  'FaceColor' - for the 3D grid
            %  'EdgeAlpha' - for the 3D grid
            %
            if ~mod(nargin,2)
                plottype = varargin{1};
                varargin = varargin(2:end);
            else
                plottype = 'default';
            end
            opt = struct('Colors',    true,    ...
                         'onGrid',    true,    ...
                         'data',      [],      ...
                         'maxWidth',  6,      ...
                         'FaceColor', 'none',  ...
                         'EdgeAlpha', 0.05);
            opt = merge_options(opt, varargin{:});
            
            lineWidth = 2;
            ne = numedges(obj.network);
            if opt.Colors
                edgeCol = {'EdgeCData', 1:ne};
                cmap = max(tatarizeMap(ne+2)-.1,0); cmap=cmap(3:end,:);
                lineWidth = 3;
            else
                edgeCol = {};
                cmap = colormap;
            end
            if ~isempty(opt.data)
                assert(numel(opt.data)==ne, ...
                    'Size of data does not match number of edges');
                lineWidth = opt.maxWidth*opt.data/max(opt.data);
            end
            
            args = {'XData',obj.network.Nodes.XData,...
                    'YData',obj.network.Nodes.YData,...
                    'ZData',obj.network.Nodes.ZData };
            switch plottype
                case {'default', 'spacegraph'}
                    % default, do nothing special
                case 'circle'
                    args = {'Layout','circle'};
                    opt.onGrid = false;
                case 'transmissibility'
                    assert(size(obj.network.Edges,2)==4, ...
                        'Network contains no transmissibilities');
                    data = obj.network.Edges.T;
                    lineWidth = opt.maxWidth*data/max(data);
                case 'porevolume'
                    assert(size(obj.network.Edges,2)==4, ...
                        'Network contains no pore volumes');
                    data = obj.network.Edges.pv;
                    lineWidth = opt.maxWidth*data/max(data);
                otherwise
                    error('Plot type not defined');
            end
            
            pg = plot(obj.network, args{:}, ...
                'LineWidth', lineWidth', edgeCol{:});
            labelnode(pg,obj.network.Nodes.well,obj.network.Nodes.name);
            pg.NodeFontSize = 10;
            colormap(gca, cmap);
            if opt.onGrid
                plotGrid(obj.G, 'FaceColor',opt.FaceColor,...
                    'EdgeAlpha',opt.EdgeAlpha);
            end
            axis off, view(2)
            if nargout==1
                omap = cmap;
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

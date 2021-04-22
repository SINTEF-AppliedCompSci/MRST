classdef FracturedDomainManager
    %FRACTURED_DOMAIN_MANAGER
    
    properties
        
    end
    
    methods
        function obj = FracturedDomainManager()
            
        end
        
        function model = addFracturedDomain(obj, model, type, region, rock, fluid, varargin)
            
            opt = struct('transfer_model', [], 'connection_cell_list', []);
            opt = merge_options(opt, varargin{:});
            
            %% Initializing fractured domain
            if(~isfield(model.G,'FracturedDomains'))
                model.G.FracturedDomains = struct();
                model.G.FracturedDomains.number_virtual_cells = 0;
                model.G.FracturedDomains.number_virtual_connections = 0;
                model.G.FracturedDomains.domains = {};
                model.G.FracturedDomains.half_trans = computeTrans(model.G, model.rock);
            end
            
            %% Creating domain
            if strcmp(type,'multi_continuum')
                transfer_model = opt.transfer_model;
                connection_cell_list = opt.connection_cell_list;
                fdom = obj.multiContinuumDomain(model, region, connection_cell_list, rock, fluid, transfer_model);
                fdom.type = 'multi_continuum';
            elseif strcmp(type,'dfm')
                fdom = obj.dfmDomain(model, region, rock, fluid);
                fdom.type = 'dfm';
            end
            
            model.G.FracturedDomains.domains{end+1} = fdom;
            model.G.FracturedDomains.number_virtual_cells = model.G.FracturedDomains.number_virtual_cells + ...
                                                            fdom.ncells;
            model.G.FracturedDomains.number_virtual_connections = model.G.FracturedDomains.number_virtual_connections + ...
                                                            fdom.nconns;
                                                        
            
            %% Increasing the number of cells in the domain
            model.G.cells.num = model.G.cells.num + ...
                                fdom.ncells;                                                        
            
            %% Updating operators
            [newdiv,newgrad,newavg,newupw,T] = obj.updateOperators(model);
            
            %% Calculating transmissibilities for DFM domain
            if strcmp(type,'dfm')
                mat_conn_ids = (ismember(model.operators.N(:,1),fdom.connections_left(:,1))...
                            & ismember(model.operators.N(:,2),fdom.connections_right(:,2)));
                T(mat_conn_ids)=0.0;
                T(fdom.global_connection_ids) = calculateTransmissibilitiesDFM(model);
                
            end
            
            %% Updating operators
            model.operators.Div = newdiv;
            model.operators.Grad = newgrad;
            model.operators.faceAvg = newavg;
            model.operators.faceUpstr = newupw;
            model.operators.T = T;
            
            
            %% 'Assembling' fluid handles 
            model = obj.assembleFluid(model);
            
            %% State function groupings
            model.PVTPropertyFunctions = FracturedDomainPVTPropertyFunctions(model);
            model.FlowPropertyFunctions = FracturedDomainFlowPropertyFunctions(model);
            model.FlowDiscretization = FracturedDomainFlowDiscretization(model);
            
            %% TODO: there is a problem when outputFluxex = 1
            model.outputFluxes = 0;
                                                        
        end
        
        function model = assembleFluid(obj, model)
            % TODO: Treat gas phase pc and kr
            
            %% Base fluid
            base_fluid = model.fluid;
            has_pc = isfield(base_fluid, 'pcOW');
            if(has_pc)
                if(iscell(base_fluid.pcOW))
                    base_fluid.pcOW = {base_fluid.pcOW{1}};
                else
                    % dummy 
                    base_fluid.pcOW = {base_fluid.pcOW};
                end
            else
                base_fluid.pcOW = {@(sw)0*sw};
            end
            
            for i = 1:length(model.G.FracturedDomains.domains)
                f = model.G.FracturedDomains.domains{i}.fluid;
                has_pc = isfield(f, 'pcOW');
                if(has_pc)
                        base_fluid.pcOW{end+1} = f.pcOW;
                else
                    base_fluid.pcOW{end+1} = @(sw)0*sw;
                end
            end
            
            %% Now the relperms. It is easier, since we know they will always exist
            % for any consistent fluid model
            if(isfield(base_fluid, 'krW'))
                if ~iscell(base_fluid.krW)
                    base_fluid.krW = {base_fluid.krW};
                end
                % Now the current domain
                f = model.G.FracturedDomains.domains{end}.fluid;
                base_fluid.krW{end+1} = f.krW;
            end
            
            if(isfield(base_fluid, 'krO'))
                if ~iscell(base_fluid.krO)
                    base_fluid.krO = {base_fluid.krO};
                end
                % Now the current domain
                f = model.G.FracturedDomains.domains{end}.fluid;
                base_fluid.krO{end+1} = f.krO;
            end
            
            if(isfield(base_fluid, 'krG'))
                if ~iscell(base_fluid.krW)
                    base_fluid.krG = {base_fluid.krG};
                end
                % Now the current domain
                f = model.G.FracturedDomains.domains{end}.fluid;
                base_fluid.krG{end+1} = f.krG;
            end
            
            if(isfield(base_fluid, 'krOW'))
                if ~iscell(base_fluid.krOW)
                    base_fluid.krOW = {base_fluid.krOW};
                end
                % Now the current domain
                f = model.G.FracturedDomains.domains{end}.fluid;
                base_fluid.krOW{end+1} = f.krOW;
            end
            
            if(isfield(base_fluid, 'krOG'))
                if ~iscell(base_fluid.krOG)
                    base_fluid.krOG = {base_fluid.krOG};
                end
                % Now the current domain
                f = model.G.FracturedDomains.domains{end}.fluid;
                base_fluid.krOG{end+1} = f.krOG;
            end
            
            %% Assigning base fluid
            model.fluid = base_fluid;
            
            %% Defining regions
            sat = ones(size(model.G.cells.centroids,1),1);
            for i = 1:length(model.G.FracturedDomains.domains)
                nc = model.G.FracturedDomains.domains{i}.ncells;
                sat = [sat; (i+1)*ones(nc,1)];
            end
            model.rock.regions.saturation = sat;
            
        end
        
        function fdom = multiContinuumDomain(obj, model, region, connection_cell_list, rock, fluid, transfer_model)
            
            %% reshaping regions
            if(size(region,2)>1)
                region = reshape(region,size(region,2),1);
            end
            if isempty(connection_cell_list)
                connection_cell_list = region;
            elseif (size(connection_cell_list,2)>1)
                connection_cell_list = reshape(connection_cell_list,size(connection_cell_list,2),1);
            end
            
            
            %% Setting rock, fluid and transfer model
            fdom.rock = rock;
            fdom.fluid = fluid;
            fdom.transfer_model = transfer_model;
             
            %% number of cells
            nc = model.G.cells.num;
            
            %% Global number of connections
            nf = size(model.operators.N,1) + model.G.FracturedDomains.number_virtual_connections;
             
            %% Cells 
            fdom.ncells = length(region);
            fdom.region = region;
            fdom.connection_cell_list = connection_cell_list;
            fdom.virtual_cells = (nc+1 : nc+fdom.ncells)';
            fdom.connections = [connection_cell_list, fdom.virtual_cells];
            fdom.nconns = size(fdom.connections,1);
            fdom.global_connection_ids = (nf+1 : nf+fdom.nconns)';
        end
        
        function fdom = dfmDomain(obj, model, edges, rock, fluid)
            
            %% reshaping regions
            if(size(edges,2)>1)
                edges = reshape(edges,size(edges,2),1);
            end
            
            %% Setting rock, fluid and transfer model
            fdom.rock = rock;
            fdom.fluid = fluid;
             
            %% number of cells
            nc = model.G.cells.num;
            
            %% Global number of connections
            nf = size(model.operators.N,1) + model.G.FracturedDomains.number_virtual_connections;
             
            %% Cells 
            fdom.ncells = length(edges);
            fdom.region = edges;
            fdom.virtual_cells = (nc+1 : nc+fdom.ncells)';
            %matrix-fracture
            fdom.connections_left = [model.G.faces.neighbors(edges,1), fdom.virtual_cells];
            %fracture-matrix
            fdom.connections_right = [fdom.virtual_cells, model.G.faces.neighbors(edges,2)];
            %fracture-fracture
            [fdom.connections_ff,fracture_node_info] = obj.fractureFractureConnections(model.G,edges,fdom.virtual_cells);
            fdom.connections = [fdom.connections_left; fdom.connections_right; fdom.connections_ff.virtual_cells];

            fdom.nconns = size(fdom.connections,1);
            fdom.global_connection_ids = (nf+1 : nf+fdom.nconns)';
            fdom.fracture_node_info = fracture_node_info;
        end
        
        function [div, grad, avg, upw, T]= updateOperators(obj, model)
            N_ = obj.assemble_connections(model);
            nf = size(N_,1);
            nc = model.G.cells.num;
            nf_last = model.G.FracturedDomains.domains{end}.nconns;
            D_ = sparse([(1:nf)'; (1:nf)'], N_, ones(nf,1)*[-1 1], nf, nc);
            div = @(x) -D_'*x;
            grad = @(x) D_*x;
            avg   = @(x) 0.5 * (x(N_(:,1)) + x(N_(:,2)));
            upw   = @(flag, x) flag.*x(N_(:, 1)) + ~flag.*x(N_(:, 2));
            T = [model.operators.T; zeros(nf_last,1)];
        end
        
        function N = assemble_connections(obj, model)
            N = model.operators.N;
            for i = 1:length(model.G.FracturedDomains.domains)
                Ndom = model.G.FracturedDomains.domains{i}.connections;
                N = [N; Ndom];
            end
        end



        function [ffcons,fracture_node_info] = fractureFractureConnections(obj, G, edges,virtual_cells)
            %   This function computes the mesh connectivity of a fracture
            %   system. Fracture nodes that have empty connection info are
            %   fracture endpoints, nodes with two connections (fracture
            %   edges) are regular fracture nodes, nodes with more than two
            %   connections are fracture intersections
            edgesNodePos = G.faces.nodePos(edges);
            frac_nodes_pos = [G.faces.nodePos(edges);G.faces.nodePos(edges)+1];
            ffcons = struct();
            % fracture-fracture connections in terms of virtual cells ids
            ffcons.virtual_cells = [];
            % loacl ids of fracture-fracture connections
            ffcons.local_ids = [];
            % fracture_node_info is necessary to identify fracture
            % intersections for star-delta treatment of the corresponding
            % connections
            fracture_node_info = struct();
            % nodes ids that belong to fracture edges
            fracture_node_info.nodes = unique(G.faces.nodes(frac_nodes_pos));
            % fracture face ids (model.G.faces) connected to the node
            fracture_node_info.conns_edges = cell(length(fracture_node_info.nodes),1);
            % fracture virtual cell local ids connected to the node
            fracture_node_info.conns_cells = cell(length(fracture_node_info.nodes),1);
            % connection positions in ffcons associated with that node
            fracture_node_info.conns_pos = cell(length(fracture_node_info.nodes),1);
            connPos = 1;
            for i=1:size(edges,1)-1
                currentEdgeNodes = G.faces.nodes(edgesNodePos(i):edgesNodePos(i)+1);
                for j=i+1:size(edges,1)
                    nextEdgeNodes = G.faces.nodes(edgesNodePos(j):edgesNodePos(j)+1);
                    [~, ui] = unique([currentEdgeNodes(1); nextEdgeNodes]);
                    if size(ui,1)<size([currentEdgeNodes(1); nextEdgeNodes],1)
                        id_pos  = find(fracture_node_info.nodes==currentEdgeNodes(1));
                        fracture_node_info.conns_edges{id_pos} = ...
                            [fracture_node_info.conns_edges{id_pos} edges(i) edges(j)];
                        fracture_node_info.conns_cells{id_pos} = ...
                            [fracture_node_info.conns_cells{id_pos} i j];
                        fracture_node_info.conns_pos{id_pos} = ...
                            [fracture_node_info.conns_pos{id_pos} connPos];
                        ffcons.virtual_cells = [ffcons.virtual_cells; ...
                            [virtual_cells(i) virtual_cells(j)]];
                        ffcons.local_ids = [ffcons.local_ids; [i j]];
                        connPos = connPos + 1;
                    end
                    [~, ui] = unique([currentEdgeNodes(2); nextEdgeNodes]);
                    if size(ui,1)<size([currentEdgeNodes(2); nextEdgeNodes],1)
                        id_pos  = find(fracture_node_info.nodes==currentEdgeNodes(2));
                        fracture_node_info.conns_edges{id_pos} = ...
                            [fracture_node_info.conns_edges{id_pos} edges(i) edges(j)];
                        fracture_node_info.conns_cells{id_pos} = ...
                            [fracture_node_info.conns_cells{id_pos} i j];
                        fracture_node_info.conns_pos{id_pos} = ...
                            [fracture_node_info.conns_pos{id_pos} connPos];
                        ffcons.virtual_cells = [ffcons.virtual_cells; ...
                            [virtual_cells(i) virtual_cells(j)]];
                        ffcons.local_ids = [ffcons.local_ids; [i j]];
                        connPos = connPos + 1;
                    end
                end
            end
            
            %eliminate duplicates
            for i=1:length(fracture_node_info.nodes)
                fracture_node_info.conns_edges{i} = unique(fracture_node_info.conns_edges{i});
                fracture_node_info.conns_cells{i} = unique(fracture_node_info.conns_cells{i});
            end
        end
    end
   
end


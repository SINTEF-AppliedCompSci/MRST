classdef Network
   properties
      network
      G
      W
   end
   methods
      function obj = Network(W,G,varargin)
         obj.G = G;
         obj.W = W;
         opt = struct('Verbose',  mrstVerbose(),...
                      'type'     ,'all_to_all',...
                      'injectors',[],...
                      'producers',[],...
                      'problem',[],...
                      'flow_filter',0,...
                      'state_number',[]);                   
         opt = merge_options(opt, varargin{:});
         
         if ~isempty(opt.injectors)&&~isempty(opt.producers)
            opt.type = 'injectors_to_producers';            
         end
         
         if isempty(opt.state_number)&&~isempty(opt.problem)
            opt.state_number = numel(opt.problem.SimulatorSetup.schedule.step.val) ;
         end                  
         
         
         nW_cells  = arrayfun(@(x)numel(x.cells), W);         
         num_nodes = sum(nW_cells);
                  
         % Define everything for each node
         node = 1;
         for iw = 1 : numel(W)
             for iw_cell = 1 : numel(W(iw).cells)
                 Nodes(node,1)     = node;
                 Well(node,1)      = iw ;
                 SubWell(node,1)   = iw_cell;
                 Well_name{node,1} = W(iw).name;
                 Well_cell{node,1} = NaN;
                 Type(node,1)      = 1;
                 cell_number       = W(iw).cells(iw_cell);
                 XData(node,1)     = obj.G.cells.centroids(cell_number,1);
                 YData(node,1)     = obj.G.cells.centroids(cell_number,2);
                 ZData(node,1)     = obj.G.cells.centroids(cell_number,3);
                 node =  node+1;
             end
         end
         
         % Create edges of the graph. Only rule is: 
         % conection with the  same well are forbiden
         A =  zeros(num_nodes);
         switch opt.type
             case 'all_to_all'
                 for iA_row =  1:num_nodes
                     well_numer  = Well(iA_row);
                     A(iA_row,(Well~=well_numer)) = 1;
                 end 
             case 'injectors_to_producers'
                 A =  zeros(num_nodes);
                 for iA_row =  1:num_nodes
                     well_numer  = Well(iA_row);
                     
                      if any(opt.injectors==well_numer)
                         conections = all(Well~=opt.injectors,2);
                      elseif any(opt.producers==well_numer)
                         conections = all(Well~=opt.producers,2);
                      else
                         warning(['\n Well %s is not in either injector or producer list.\n',...
                             'This well will be ignored'], Well_name{iA_row});
                      end 
                     A(iA_row,conections) = 1;
                 end                 
             case {'fd_preprocessor','fd_postprocessor'}
                 if strcmp(opt.type,'fd_preprocessor')
                     state = [];
                     pressure_field = 'pressure';
                 else
                     [~, states, ~] = getPackedSimulatorOutput(opt.problem,...
                         'readWellSolsFromDisk',false,...
                         'readReportsFromDisk',false);
                     state = states{opt.state_number};
                     pressure_field = 'bhp';
                 end
                 % TODO what happens with more controls in the schedule
                 [state, diagnostics] = computePressureAndDiagnostics(opt.problem.SimulatorSetup.model,...
                     'wells', opt.problem.SimulatorSetup.schedule.control(1).W,...
                     'state',state);
                  %Determine Well pair indices
                    [IP_indices]=find(diagnostics.wellCommunication > opt.flow_filter);
                    [I,P]=find(diagnostics.wellCommunication > opt.flow_filter);
                    
                    % TODO: Its assume injector are first and then comes
                    % producers. This may cause some problem in case
                    % injector and producer does not follows this.
                    
                    P_indx = P +size(diagnostics.wellCommunication,1); 
                    for wp = 1:length(IP_indices)        
                            iwp =  IP_indices(wp);
                            fluxes(wp,1) = diagnostics.wellCommunication(I(wp),P(wp));
                            DP(wp,1) = state.wellSol(I(wp)).(pressure_field)-... % Injector pressure
                                       state.wellSol(P_indx(wp)).(pressure_field);          %
                            T(wp,1)  =   fluxes(wp,1)/DP(wp,1);
                            pv(wp,1) = diagnostics.WP.vols(iwp);
                            A(I(wp),P_indx(wp)) = 1;
                            A(P_indx(wp),I(wp)) = 1;
                    end          
             otherwise
                 error('\nType of network: %s is not implemented\n', opt.type);     
         end
         
         
         obj.network = graph(A);
         obj.network.Nodes.Nodes     = Nodes;
         obj.network.Nodes.Well      = Well;
         obj.network.Nodes.SubWell   = SubWell;
         obj.network.Nodes.Type      = Type;
         obj.network.Nodes.Well_name = Well_name;
         obj.network.Nodes.Well_cell = Well_cell;
         obj.network.Nodes.XData     = XData;
         obj.network.Nodes.YData     = YData;
         obj.network.Nodes.ZData     = ZData;
         if any(strcmp(opt.type,{'fd_preprocessor','fd_postprocessor'}))
             obj.network.Edges.Transmissibility = T;  
             obj.network.Edges.PoreVolume = pv;
         end
      end
      
      function f = plotNetwork(obj,varargin)
          opt = struct('FaceColor', 'none',...
                       'EdgeAlpha'   ,'0.1',...
                       'NetworkLineWidth',2*ones(numedges(obj.network),1));                   
          opt = merge_options(opt, varargin{:});
          
          f= plotGrid(obj.G, 'FaceColor',opt.FaceColor,...
                             'EdgeAlpha',opt.EdgeAlpha); view(2);
          
          hold on, pg =  plot(obj.network,...
                             'XData',obj.network.Nodes.XData,...
                             'YData',obj.network.Nodes.YData,...
                             'ZData',obj.network.Nodes.ZData,...
                             'LineWidth',opt.NetworkLineWidth);
                  labelnode(pg,obj.network.Nodes.Well,obj.network.Nodes.Well_name)
          hold off;

          
      end
      
   end
end
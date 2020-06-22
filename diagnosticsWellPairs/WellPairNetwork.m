classdef WellPairNetwork
    properties
        dataMRST
        q
        WP %Wellpair structure
        D  %diagnostic
        model  %MRST Model
        schedule
        wps_all
        wps
        wp_flow_filter =  100*stb/day;
        I_filter
        Graph
        I_Graph_index
    end
    methods
        
        function DD = WellPairNetwork(model,schedule,states,state,wellSols,varargin)
            
            
            
            DD.model    = model;
            DD.schedule = schedule;
            [DD.dataMRST,DD.WP,DD.D]...
                = computeWellPairsData(model,schedule,states,state,wellSols);
            
            % Taking information for each wellpair conection
           for wp_index = 1 : size(DD.WP.pairIx,1)
               DD.wps_all{wp_index} = WellPair(DD.dataMRST,DD.WP,DD.D,...
                                    wp_index,model.G,schedule.control.W);
           end
           
           DD = DD.filter_wps(DD.wp_flow_filter);    
        end 
        
        
        function  DD = filter_wps(DD,wp_flow_filter)
           DD.wp_flow_filter = wp_flow_filter;
           
           DD.I_filter = find(abs(DD.dataMRST.BFPD(end,:))> DD.wp_flow_filter);
           
           DD.wps = DD.wps_all(DD.I_filter);      
           
           for i =  1:  numel(DD.wps)
            edges(i,:)= [DD.wps{i}.WellSolsIx_inj , DD.wps{i}.WellSolsIx_prod];
           end            

            
            
            ne = size(edges,1);
    
            DD.Graph = graph(edges(:,1),edges(:,2),[]);   
            [~, DD.I_Graph_index] = sortrows(edges); % to be used for ordering the wellpair conections          
        end
        
        
        function f=plotWellPairVolume(DD,wp_index)
            f=DD.wps{wp_index}.plotWellPairVolume(DD.model.G,DD.model.rock.perm(:,1));
        end
                      
        function f = plotWellPairConnections(DD)
            
            for i =  1:numel(DD.wps)
                pv(i) = DD.wps{i}.volume;
                TT(i) = DD.wps{i}.Tr;
            end
            n_wells  = numel(DD.schedule.control.W);
            for i = 1:n_wells
                cell_number = DD.schedule.control.W(i).cells(1);
                XData(i) = DD.model.G.cells.centroids(cell_number,1);
                YData(i) = DD.model.G.cells.centroids(cell_number,2);
                ZData(i) = DD.model.G.cells.centroids(cell_number,3);
            end


            
            subplot(1,2,1);

                h2= plotGrid(DD.model.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(DD.Graph, 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*TT/max(TT))
                        labelnode(pg,[1:n_wells],{DD.schedule.control.W.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                title("Transmisibility")


                subplot(1,2,2);
                h2= plotGrid(DD.model.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg=plot(DD.Graph,'r','XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*pv/max(pv))
                         labelnode(pg,[1:n_wells],{DD.schedule.control.W.name})
                hold off;
                pg.NodeFontSize= 15;
                axis off ; 
                title("Pore volume")
        end
        
        function f = plotProducer(DD,prod_indx)    
               f = figure('Name',DD.WP.prod(prod_indx).name);
               inj_index = get_injector_index(DD,prod_indx);
               colorstring = 'kbgry';
               linewidth = 2;
               iter =1;
               
                error('plotProducer would be re-implemented')           

        end
        
        
%%
            
       function f = plotProducers(DD)    
           error('plotSaturation would be re-implemented')           
        end
        
        
       
        function f = plotWellPairsData(DD,varargin)    
             
             opt = struct('subplot',[4 5]);
             
             opt  = merge_options(opt, varargin{:});
             m = opt.subplot(1);
             n = opt.subplot(2);
             mn_subplot= min(n*m,numel(DD.wps));

             if mn_subplot<numel(DD.wps)
                 warning('Several well pairs  will not be display. To include more rows n or  collumns m in the figures add a subplot option:  DD.plotWellPairsData("subplot",[m,n]), ');
             end
             tt =  cumsum(DD.schedule.step.val);
             f = figure('Name','Flow rate');
                for i = 1 : mn_subplot
                    subplot(m,n,i)
                    plot(tt/day, -DD.wps{i}.data.BFPD*day/stb,'k')
                     title(DD.wps{i}.wellpair_name);
                end   
                
             f = figure('Name','Volume');
                for i = 1 : mn_subplot
                    subplot(m,n,i)
                    plot(tt/day, DD.dataMRST.VV(:,DD.I_filter(i)),'k')
                     title(DD.wps{i}.wellpair_name);
                end
                
             f = figure('Name','Average saturation');
                for i = 1 :mn_subplot
                    subplot(m,n,i)
                    plot(tt/day, DD.wps{i}.data.s_avg,'k')
                     title(DD.wps{i}.wellpair_name);
                end     

        end
            
        function inj_index = get_injector_index(DD,inj_producer)
            inj_index = []; iter = 0;
            for wp_index = 1 :  size(DD.wps,2)
                if (DD.wps{wp_index}.wellpairIx_prod == inj_producer)
                    inj_index(iter+1) = wp_index;
                    iter= iter+1;
                end
            end           
        end
        
    end
end
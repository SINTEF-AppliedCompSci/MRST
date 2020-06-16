classdef WellPair
    properties
       volume 
       length
       lambda
       fw
       Tr
       data
       domain
       c
       schedule
       u % Velocity of the faturation front
       dx
       num_1D_cells
       s % Saturation at wellpair
       q % flow at well pair
       qw
       WellSolsIx_inj % Injetor index in WellSols
       WellSolsIx_prod % Pro ducer index in WellSols       
       
       wellpairIx_inj  % Injetor index in WP
       wellpairIx_prod % Producer index in WP
       
       wellpair_name

    end  
    methods 
        function wp = WellPair(data,WP,D, wp_index ,G , W , varargin)
            wp.volume =  WP.vols(wp_index);
            
            wp.wellpairIx_inj  =  WP.pairIx(wp_index,1) ; % Injetor index in WP
            wp.wellpairIx_prod =  WP.pairIx(wp_index,2) ; % Producer index in WP
            
            wp.wellpair_name = WP.pairs{wp_index};
            
            % Intersection between sweep area from injector and drained area from producer
            wp.c  =  D.itracer(:,wp.wellpairIx_inj).*D.ptracer(:,wp.wellpairIx_prod);
            wp.domain = find(wp.c>0);
            wp.WellSolsIx_inj  =  D.inj(  WP.pairIx(wp_index,1) ); % Injetor index in WellSols
            wp.WellSolsIx_prod =  D.prod( WP.pairIx(wp_index,2) ); % Producer index in WellSols
                     
            cellnumber_inj  =  W(wp.WellSolsIx_inj).cells;  % Cell numbers
            cellnumber_prod =  W(wp.WellSolsIx_prod).cells; % Cell numbers

            position_inj    = G.cells.centroids(cellnumber_inj,:); % Positions                     
            position_prod   = G.cells.centroids(cellnumber_prod,:); % Positions

            wp.length       = norm(position_prod(end,:)-position_inj(end,:)); % TODO this line should be inprove to consider horizontal wells and vertical wells with diferent completions
            
            wp.num_1D_cells    = 30;
            
            wp.dx              = wp.length /wp.num_1D_cells;
            
            
            
            wp.lambda = @(s)interp1(data.s_avg(:,wp_index),...
                                    data.lambda(:,wp_index),...
                                    s,'nearest','extrap') ;
                                
%             wp.fw    = @(s)interp1(data.s_avg(:,wp_index),...
%                                    data.fw_avg(:,wp_index),...
%                                    s,'nearest','extrap') ;  
                               
            wp.Tr = data.Tr(wp_index);                
            
            wp.data.BFPD = data.BFPD(:,wp_index);
            wp.data.BWPD = data.BWPD(:,wp_index);
            wp.data.BOPD = data.BOPD(:,wp_index);
            
            wp.data.VV       = data.VV(:,wp_index);

            wp.data.DP       = data.DP(:,wp_index);
            wp.data.BHP_inf  = data.BHP_inf(:,wp_index);
            wp.data.BHP_prod = data.BHP_prod(:,wp_index);
            
            wp.data.s_avg    = data.s_avg(:,wp_index);
            wp.data.s_avg_0    = data.s_avg_0(wp_index);

            wp.data.fw_avg   = data.fw_avg(:,wp_index);
            wp.data.lambda   = data.lambda(:,wp_index);
            
            wp.schedule      = data.schedule;
            

        end        
        
        
         function f = plotWellPairVolume(wp,G,data,varargin)
             f =  figure;
             plotGrid(G,'FaceColor','none', 'EdgeAlpha', 0.2);
             view(3), axis tight
             plotCellData(G, data,wp.domain);
             plotWell(G,wp.schedule.control.W([wp.WellSolsIx_inj,wp.WellSolsIx_prod]));
             
         end
        function f =plotWellPairData(wp,varargin)
            f = figure;
            
            subplot(1,3,1); plot(wp.data.s_avg,wp.data.fw_avg,'.')
                   ylabel('f_w');  xlabel('S_w');
            %legend(WP_pairs_names,'Orientation','vertical','Location','SouthEast')

            subplot(1,3,2); plot(cumsum(wp.schedule.step.val)/day,wp.data.DP/barsa,'.')
                   ylabel('Dp [bar]');  xlabel('time [days]');

            subplot(1,3,3);             
                   plot(cumsum(wp.schedule.step.val)/day, -wp.data.BFPD*day/stb,'.k',...
                        cumsum(wp.schedule.step.val)/day, -wp.data.BFPD.*wp.data.fw_avg*day/stb,'.b',...
                        cumsum(wp.schedule.step.val)/day, -wp.data.BFPD.*(1-wp.data.fw_avg)*day/stb,'.r');
                    ylabel('STB');  xlabel('time [days]');
            legend('Well pair BFPD','Well pair BWPF','Well pair BOPD','Orientation','vertical','Location','SouthEast')
        end
                                                  
    end
end      



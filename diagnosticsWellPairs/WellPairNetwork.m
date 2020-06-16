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
           
           DD.I_filter = find(abs(DD.dataMRST.BFPD(end,:))> DD.wp_flow_filter);
           DD.wps = DD.wps_all(DD.I_filter);           
        end 
        
        function f=plotWellPairVolume(DD,wp_index)
            f=DD.wps{wp_index}.plotWellPairVolume(DD.model.G,DD.model.rock.perm(:,1));
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
           error('To be re-implemented')           
        end
        
        
       
         %% Plot and compare saturation given by the MRST model using flow diagnostic and the DD model
        function f = plotSaturation(DD,wp_index)    
            
            
               f = figure('Name',DD.wps{wp_index}.wellpair_name);
                        
               error('plotSaturation would be re-implemented')           


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
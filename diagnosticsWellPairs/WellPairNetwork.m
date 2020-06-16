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
        

        
        function f =plotWellPairModel(DD,wp_index,varargin)
            f =plotWellPairModel(DD.wps{wp_index},varargin);
        end
        
        function f = plotDataDrivenModel(DD)
            f = figure;
            % plotting flow rate given by the data driven model
            plot(cumsum(DD.schedule.step.val)/day,DD.q*day/stb);
            %Legend for each producer
            legend(DD.WP.prod.name);
            ylabel('total flux [STB]');  xlabel('time [days]');
            
            hold on
            % plotting flow rate given by the MRST model
            plot(cumsum(DD.schedule.step.val)/day,(DD.dataMRST.q_prod)*day/stb,'.');
            
                        legend(DD.WP.prod.name);

            hold off            
        end     
        
        function f = plotProducer(DD,prod_indx)    
               f = figure('Name',DD.WP.prod(prod_indx).name);
               inj_index = get_injector_index(DD,prod_indx);
               colorstring = 'kbgry';
               linewidth = 2;
               iter =1;
               
               q_prod_model  = 0*DD.wps{inj_index(1)}.q';
               qw_prod_model = q_prod_model;
               qo_prod_model = q_prod_model;
               Q_prod_model  = 0;
               Qw_prod_model  = 0;
               Qo_prod_model  = 0;
               
               Q_prod_data  = 0;
               Qw_prod_data  = 0;
               Qo_prod_data  = 0;
               
%                q_prod_data  =   DD.dataMRST.q_prod(:,prod_indx);
%                qw_prod_data =  -DD.wps{inj_index(1)}.data.BWPD;
%                qo_prod_data =  -DD.wps{inj_index(1)}.data.BOPD;   
               
               
               %qw_prod_data =  -DD.wps{inj_index(1)}.data.BWPD;
               %qo_prod_data =  -DD.wps{inj_index(1)}.data.BOPD;
               
               for i= 1 :length(inj_index)
                   
                   
               wp_index = inj_index(i);
               lineproperties = {'Color', colorstring(i),'LineStyle',':'};

               
               S_model =  DD.wps{wp_index}.s';
               
               
               Fw_model = DD.wps{wp_index}.fw.val(S_model);
               Fo_model = 1-Fw_model;
               

                q_prod_model  =  DD.wps{wp_index}.q';
                qw_prod_model =  DD.wps{wp_index}.qw';
                qo_prod_model =  (DD.wps{wp_index}.q'-DD.wps{wp_index}.qw');
                
                fw = DD.wps{wp_index}.data.fw_avg;
                q_prod_data  =   -DD.wps{wp_index}.data.BFPD; 
                qw_prod_data =   -Fw_model.*DD.wps{wp_index}.data.BFPD;
                qo_prod_data =   -Fo_model.*DD.wps{wp_index}.data.BFPD;   
                

              subplot(2,3,1); line(cumsum(DD.schedule.step.val)/day, DD.wps{wp_index}.data.s_avg,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,S_model,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
                                    ylabel('S_avg');  xlabel('time [days]');                
                                    
%                subplot(1,2,2); line(cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.data.fw_avg,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line(cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.fw.val(S_model),'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                     ylabel('f_w');  xlabel('time [days]'); 
%                subplot(1,2,3); line( DD.dataMRST.s_avg(:,wp_index),DD.wps{wp_index}.data.fw_avg,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line( S_model,DD.wps{wp_index}.fw.val(S_model),'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                     ylabel('f_w');  xlabel('S_avg'); 
%                
%                subplot(2,3,4); line(cumsum(DD.schedule.step.val)/day, DD.wps{wp_index}.data.lambda,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line(cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.lambda(S_model),'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                     ylabel('\lambda');  xlabel('time [days]');
                                    
%                subplot(2,3,4); line(cumsum(DD.schedule.step.val)/day,-DD.wps{wp_index}.data.BFPD*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line(cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.q*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                     ylabel('BFPD');  xlabel('time [days]')
%                                     
%                subplot(2,3,2); line(cumsum(DD.schedule.step.val)/day,-Fw_model.*DD.wps{wp_index}.data.BFPD*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line(cumsum(DD.schedule.step.val)/day,Fw_model.*DD.wps{wp_index}.q'*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                   ylabel('BWPD');  xlabel('time [days]')
%                                   
%                subplot(2,3,3); line(cumsum(DD.schedule.step.val)/day,-Fo_model.*DD.wps{wp_index}.data.BFPD*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line(cumsum(DD.schedule.step.val)/day,Fo_model.*DD.wps{wp_index}.q'*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                   ylabel('BOPD');  xlabel('time [days]')                                  
%                
                                  
               subplot(2,3,4); line(cumsum(DD.schedule.step.val)/day,-DD.wps{wp_index}.data.BFPD*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,q_prod_model*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
                                    ylabel('BFPD');  xlabel('time [days]')
                                    
               subplot(2,3,2); line(cumsum(DD.schedule.step.val)/day,qw_prod_data*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,qw_prod_model*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
                                  ylabel('BWPD');  xlabel('time [days]')
                                  
               subplot(2,3,3); line(cumsum(DD.schedule.step.val)/day,qo_prod_data*day/stb,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,qo_prod_model*day/stb,'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
                                  ylabel('BOPD');  xlabel('time [days]')                                       
%                subplot(2,3,6); line( DD.dataMRST.s_avg(:,wp_index),DD.wps{wp_index}.data.lambda,'Color', colorstring(i),'LineStyle',':','LineWidth',linewidth)
%                                line( S_model,DD.wps{wp_index}.lambda(S_model),'Color', colorstring(i),'LineStyle','-','LineWidth',linewidth)
%                                     ylabel('\lambda');  xlabel('S_avg'); 
                Leyenda{iter}=['INJEC', num2str(DD.wps{wp_index}.wellpairIx_inj)];
                Leyenda{iter+1}='' ;
                iter =iter+2;
                
                Q_prod_model  = Q_prod_model  +q_prod_model;
                Qw_prod_model = Qw_prod_model +qw_prod_model;
                Qo_prod_model = Qo_prod_model +qo_prod_model;
                
                Q_prod_data  = Q_prod_data   +q_prod_data ;
                Qw_prod_data = Qw_prod_data  +qw_prod_data ;
                Qo_prod_data = Qo_prod_data  +qo_prod_data ;
                
               end
               subplot(2,3,1);legend(Leyenda);
               subplot(2,3,2);legend(Leyenda);
               subplot(2,3,3);legend(Leyenda);
               subplot(2,3,4);legend(Leyenda);
               
            % It does not compare with raw production data but with
            % the production data between each well pair conection obtaine
            % by flow diagnostic. It can be seen as cheating,
            % but  we are calibrating against production data from flow
            % diagnostic so we compare against that data.
               
               subplot(2,3,5); line(cumsum(DD.schedule.step.val)/day,Q_prod_data*day/stb,'Color', 'k','LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,Q_prod_model*day/stb,'Color', 'k','LineStyle','-','LineWidth',linewidth)
                               
                               line(cumsum(DD.schedule.step.val)/day,Qw_prod_data*day/stb,'Color', 'b','LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,Qw_prod_model*day/stb,'Color', 'b','LineStyle','-','LineWidth',linewidth)

                               line(cumsum(DD.schedule.step.val)/day,Qo_prod_data*day/stb,'Color', 'g','LineStyle',':','LineWidth',linewidth)
                               line(cumsum(DD.schedule.step.val)/day,Qo_prod_model*day/stb,'Color', 'g','LineStyle','-','LineWidth',linewidth)
 ylabel('BPD');  xlabel('time [days]')
 
 legend('BFPD MRST model','BFPD data-driven model','BWPD MRST model','BWPD data-driven model','BOPD MRST model','BOPD data-driven model');
 xlim([0 3600])    
% title('Data driven model with no calibration')
     title('After calibration')
                    %title('Before calibration')
                 
        end
        
        
        %% This fucntion compares the production results between MRST model and DD model at each producer
            % It does not compare with raw production data but with
            % the production data between each well pair conection obtaine
            % by flow diagnostic. It can be seen as cheating,
            % but  we are calibrating against production data from flow
            % diagnostic so we compare against that data.
            
       function f = plotProducers(DD)    
           
           
                for prod_indx =1 :4
                            inj_index = get_injector_index(DD,prod_indx);
                       
                       colorstring = 'kbgry';
                       linewidth = 2;

                       q_prod_model  = 0*DD.wps{inj_index(1)}.q';
                       qw_prod_model = q_prod_model;
                       qo_prod_model = q_prod_model;
                       
                       q_prod_data  =  q_prod_model;
                       qw_prod_data =  q_prod_model;
                       qo_prod_data =  q_prod_model;
                       
                       Q_prod_model  = 0;
                       Qw_prod_model  = 0;
                       Qo_prod_model  = 0;

                       Q_prod_data  = 0;
                       Qw_prod_data  = 0;
                       Qo_prod_data  = 0;
                       
                       q_prod_data  =   DD.dataMRST.q_prod(:,prod_indx);

                       for i= 1 :length(inj_index)
                       wp_index = inj_index(i);

                       S_model =  DD.wps{wp_index}.s';
                       Fw_model = DD.wps{wp_index}.fw.val(S_model);
                       
                       Fo_model = 1-Fw_model;

                       %q_prod_data
                        q_prod_model  =  DD.wps{wp_index}.q';
                        qw_prod_model =  DD.wps{wp_index}.qw';
                        qo_prod_model =  (DD.wps{wp_index}.q'-DD.wps{wp_index}.qw');

                        fw = DD.wps{wp_index}.data.fw_avg;
                        q_prod_data  =   -DD.wps{wp_index}.data.BFPD; 
                        qw_prod_data =   -Fw_model.*DD.wps{wp_index}.data.BFPD;
                        qo_prod_data =   -Fo_model.*DD.wps{wp_index}.data.BFPD;   

                Q_prod_model  = Q_prod_model  +q_prod_model;
                Qw_prod_model = Qw_prod_model +qw_prod_model;
                Qo_prod_model = Qo_prod_model +qo_prod_model;
                
                Q_prod_data  = Q_prod_data   +q_prod_data ;
                Qw_prod_data = Qw_prod_data  +qw_prod_data ;
                Qo_prod_data = Qo_prod_data  +qo_prod_data ;
                       end


                       subplot(2,2,prod_indx); line(cumsum(DD.schedule.step.val)/day,Q_prod_data*day/stb,'Color', 'k','LineStyle',':','LineWidth',linewidth)
                                       line(cumsum(DD.schedule.step.val)/day,Q_prod_model*day/stb,'Color', 'k','LineStyle','-','LineWidth',linewidth)

                                       line(cumsum(DD.schedule.step.val)/day,Qw_prod_data*day/stb,'Color', 'r','LineStyle',':','LineWidth',linewidth)
                                       line(cumsum(DD.schedule.step.val)/day,Qw_prod_model*day/stb,'Color', 'r','LineStyle','-','LineWidth',linewidth)

                                       line(cumsum(DD.schedule.step.val)/day,Qo_prod_data*day/stb,'Color', 'g','LineStyle',':','LineWidth',linewidth)
                                       line(cumsum(DD.schedule.step.val)/day,Qo_prod_model*day/stb,'Color', 'g','LineStyle','-','LineWidth',linewidth)
                        title(['PROD', num2str(prod_indx)])
                       legend('BFPD','','BWPD','','BOPD','');
                        
                                            ylabel('BFPD');  xlabel('time [days]')
                end
                                    
        end
        
        
       
         %% Plot and compare saturation given by the MRST model using flow diagnostic and the DD model
        function f = plotSaturation(DD,wp_index)    
            
            
               f = figure('Name',DD.wps{wp_index}.wellpair_name);
                        
               S_model =  DD.wps{wp_index}.s';
               subplot(2,3,1); plot(cumsum(DD.schedule.step.val)/day, DD.wps{wp_index}.data.s_avg,'.',...
                                    cumsum(DD.schedule.step.val)/day,S_model,'--')
                                    ylabel('S_avg');  xlabel('time [days]');
                                    
               subplot(2,3,2); plot(cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.data.fw_avg,'.',...
                                    cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.fw.val(S_model),'--')
                                    ylabel('f_w');  xlabel('time [days]'); 
               subplot(2,3,3); plot( DD.dataMRST.s_avg(:,wp_index),DD.wps{wp_index}.data.fw_avg,'.',...
                                    S_model,DD.wps{wp_index}.fw.val(S_model),'--')
                                    ylabel('f_w');  xlabel('S_avg'); 
               
               subplot(2,3,4); plot(cumsum(DD.schedule.step.val)/day, DD.wps{wp_index}.data.lambda,'.',...
                                    cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.lambda(S_model),'--')
                                    ylabel('\lambda');  xlabel('time [days]');
                                    
               subplot(2,3,5); plot(cumsum(DD.schedule.step.val)/day,-DD.wps{wp_index}.data.BFPD*day/stb,'.',...
                                    cumsum(DD.schedule.step.val)/day,DD.wps{wp_index}.q*day/stb,'--')
                                    ylabel('BFPD');  xlabel('time [days]')
               
               subplot(2,3,6); plot( DD.dataMRST.s_avg(:,wp_index),DD.wps{wp_index}.data.lambda,'.',...
                                    S_model,DD.wps{wp_index}.lambda(S_model),'--')
                                    ylabel('\lambda');  xlabel('S_avg'); 
                                    
               legend('flow-diagnostic','DD model');

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
function [data,WP,D] = computeWellPairsData(model,schedule,states,state,wellSols)

%TODO: check if the diagnostic have been calculated before
  

     W = schedule.control(1).W;
     

for step = 1 : size(states,1)
    
     % Diagnostics
     D = computeTOFandTracer(states{step}, model.G, model.rock, 'wells', W, ...
                            'computeWellTOFs', true, 'firstArrival', true);                    
     % Well pairs                    
     WP = computeWellPairs(states{step}, model.G, model.rock, W, D);
    
     
     for  each_prod=1:size(WP.prod,2)
        q_prod(step,each_prod)= -sum(sum(WP.prod(each_prod).alloc));
        
     end
     
     %Well Solution
     wellSols_step = wellSols(step,1);


    % fw_1 = WP_fractional_flow_teoretical(model.fluid,states{step});

      if step==1
          for well_pair  = 1:size(WP.pairIx,1)
                s_avg_0(well_pair) = WP_pair_saturation(well_pair,state,model,WP,D);  
          end
      end
     

     for well_pair  = 1:size(WP.pairIx,1)
         
        wellpairIx_inj  =  WP.pairIx(well_pair,1) ; % Injetor index in WP
        wellpairIx_prod =  WP.pairIx(well_pair,2) ; % Producer index in WP
     
        % Intersection between sweep area from injector and drained area from producer
        c  = D.itracer(:,wellpairIx_inj).*D.ptracer(:,wellpairIx_prod);
        
        cells = find(c>0);
         
        %Well pair pressuredrop    
        [DP(step,well_pair)  , p_inj(step,well_pair), p_prod(step,well_pair)]= WP_pair_pressure_drop(well_pair,wellSols_step,WP,D) ;
        %Wellpais Saturation
        [s_avg(step,well_pair),VV(step,well_pair)] = WP_pair_saturation(well_pair,states{step}(1),model,WP,D);  
        %Wellpair Flux
        [BFPD(step,well_pair),q_inj(step,well_pair)]     =  WP_pair_flux(well_pair,WP);  % [m^3 /s] 

       
    %    fw_1 = WP_fractional_flow_teoretical(model.fluid,states{step});
     %   fw_avg_teorico(step,well_pair) = mean(fw_1(find(c>0)));% [%]
        
       [fw_avg_data(step,well_pair),...
        BWPD(step,well_pair),...
        BOPD(step,well_pair)] =  WP_fractional_flow_flowdata(well_pair,wellSols_step,WP,D);% [%]
        
        if step>1
%            lambda(step,well_pair)                    = -BFPD(step,well_pair) / (T_0(well_pair)*DP(step,well_pair));
        %elseif step == 1
            T_0(well_pair)                      = -BFPD(step,well_pair) / (DP(1,well_pair) ); %*mob(0)/mob(1)
            lambda(1,well_pair) = 1.0 ;       
        end
        
        %[L(step,well_pair),Ev,tD]                 = WP_pair_Lorenz(well_pair,model,WP,D);
     end
 
end 


data.cells     =cells;

% Production data from each well pair
data.BFPD     = BFPD;
data.BWPD     = BWPD;
data.BOPD     = BOPD;


data.q_prod   = q_prod;

data.DP       = DP;
data.VV       = VV;
data.BHP_inf  = p_inj;
data.BHP_prod = p_prod;

% Volumetric data from each well pair domain
data.Tr       = T_0;
data.s_avg    = s_avg;
data.s_avg_0  = s_avg_0;

data.fw_avg   = fw_avg_data;
data.lambda   = lambda;

data.schedule = schedule;


end
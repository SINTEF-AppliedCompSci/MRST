function [data,WP,D] = computeWellPairsData(model,schedule,states,state,wellSols)

%TODO: check if the diagnostic have been calculated before
  

    
     
n_controls = size(schedule.control,1);
for step = 1 : size(states,1)
    
     if n_controls >= step
         W = schedule.control(step).W;
     else 
         W = schedule.control.W;
     end
     % Diagnostics
     D{step} = computeTOFandTracer(states{step}, model.G, model.rock, 'wells', W, ...
                            'computeWellTOFs', true, 'firstArrival', true);                    
     % Well pairs                    
     WP{step} = computeWellPairs(states{step}, model.G, model.rock, W, D{step});
    
     
     for  each_prod=1:size(WP{step}.prod,2)
        q_prod(step,each_prod)= -sum(sum(WP{step}.prod(each_prod).alloc));
        
     end
     
     %Well Solution
     wellSols_step = wellSols(step,1);


    % fw_1 = WP_fractional_flow_teoretical(model.fluid,states{step});

      if step==1
          for well_pair  = 1:size(WP{step}.pairIx,1)
                s_avg_0(well_pair) = WP_pair_saturation(well_pair,state,model,WP{step},D{step});  
          end
      end
     

     for well_pair  = 1:size(WP{step}.pairIx,1)
         
        wellpairIx_inj  =  WP{step}.pairIx(well_pair,1) ; % Injetor index in WP
        wellpairIx_prod =  WP{step}.pairIx(well_pair,2) ; % Producer index in WP
     
        % Intersection between sweep area from injector and drained area from producer
        c  = D{step}.itracer(:,wellpairIx_inj).*D{step}.ptracer(:,wellpairIx_prod);
        
        cells = find(c>0);
         
        %Well pair pressuredrop    
        [DP(step,well_pair)  , p_inj(step,well_pair), p_prod(step,well_pair)]= WP_pair_pressure_drop(well_pair,wellSols_step,WP{step},D{step}) ;
        %Wellpais Saturation
        [s_avg(step,well_pair),VV(step,well_pair)] = WP_pair_saturation(well_pair,states{step}(1),model,WP{step},D{step});  
        %Wellpair Flux
        [BFPD(step,well_pair),q_inj(step,well_pair)]     =  WP_pair_flux(well_pair,WP{step});  % [m^3 /s] 

       
    %    fw_1 = WP_fractional_flow_teoretical(model.fluid,states{step});
     %   fw_avg_teorico(step,well_pair) = mean(fw_1(find(c>0)));% [%]
       % fprintf("Step %i wellpair %i- Size %i\n",step,well_pair,size(WP{step}.pairIx,1))
       [fw_avg_data(step,well_pair),...
        BWPD(step,well_pair),...
        BOPD(step,well_pair)] =  WP_fractional_flow_flowdata(well_pair,wellSols_step,WP{step},D{step});% [%]
        
%         if step>48
%             lambda(step,well_pair)                    = -BFPD(step,well_pair) / (T_0(well_pair)*DP(step,well_pair));
%         else
            T_0(step,well_pair)                      = -BFPD(step,well_pair) / (DP(1,well_pair) ); %*mob(0)/mob(1)
%             lambda(1,well_pair) = 1.0 ;       
%         end
        
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
%data.s_avg    = s_avg;
%data.s_avg_0  = s_avg_0;

data.fw_avg   = fw_avg_data;
%data.lambda   = lambda;

data.schedule = schedule;


end
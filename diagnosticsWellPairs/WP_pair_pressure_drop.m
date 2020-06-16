%%
%This function calculate the pressure drop between a well pair Injector
%Producer

% well_pair is the wellpair index
%
% wellSols_step Contain the Well solution from the simulaation
%
% WP contains the intex of the well pair and 
%
% D contains the global well index necesary for accesing the well solution 
% bhp

function [DP,p_inj,p_prod] = WP_pair_pressure_drop(well_pair,wellSols_step,WP,D)     

     %We need WP and D to the index pair in the wellSols


     pairIx_inj  =  WP.pairIx(well_pair,1) ; % Injetor index in WP
     pairIx_prod =  WP.pairIx(well_pair,2) ; % Producer index in WP
     
     Ix_inj  = D.inj(  WP.pairIx(well_pair,1) ); % Injetor index in WellSols
     Ix_prod = D.prod( WP.pairIx(well_pair,2) ); % Pruducer index in WellSols

 DP = wellSols_step{1,1}(Ix_inj).bhp - wellSols_step{1,1}(Ix_prod).bhp; %Pascal  [Pa] =  [kg⋅m^−1 /  s^−2]
if nargout == 3
 p_inj = wellSols_step{1,1}(Ix_inj).bhp;  %Pascal  [Pa] =  [kg⋅m^−1 /  s^−2]
 p_prod = wellSols_step{1,1}(Ix_prod).bhp; %Pascal  [Pa] =  [kg⋅m^−1 /  s^−2]
end
end
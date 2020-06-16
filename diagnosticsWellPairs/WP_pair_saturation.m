%%
%This function calculate the average saturation between a well pair Injector
%Producer

% well_pair is the wellpair index
%
% state Contain the reservoir state
%
% model Contain the reservoir state
%
% WP contains the intex of the well pair and 
%
% D contains the global well index necesary for accesing the well solution 
% bhp

function [s,varargout] = WP_pair_saturation(well_pair,state,model,WP,D)     

     
  % average saturation in interaction region 1-1
    pv = model.operators.pv;  % Por volume
      
     pairIx_inj  =  WP.pairIx(well_pair,1) ; % Injetor index in WP
     pairIx_prod =  WP.pairIx(well_pair,2) ; % Producer index in WP
     
    % Intersection between sweep area from injector and drained area from producer
    c  = D.itracer(:,pairIx_inj).*D.ptracer(:,pairIx_prod);  
    % fluid volume involved
    v = sum(c.*pv);
        
    s = sum(c.*pv.*state.s(:,1))/ v;
    if nargout > 1
        varargout{1} = v; %Volumen
    end
end
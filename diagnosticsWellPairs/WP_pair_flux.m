%%
%This function calculate volumetric flux between a well pair Injector
%Producer


% well_pair is the wellpair index
%
% WP contains the intex of the well pair and 


function [f_prod,f_inj] = WP_pair_flux(well_pair,WP)     

     pairIx_inj  =  WP.pairIx(well_pair,1) ; % Injetor index in WP
     pairIx_prod =  WP.pairIx(well_pair,2) ; % Producer index in WP
     

% flux given by the injector pairIx_inj to producer pairIx_prod
 f_inj   =  sum(WP.inj(pairIx_inj).alloc(: ,  pairIx_prod)); % [m^3 /s]
 
% flux receive by the producer pairIx_prod from injector pairIx_inj
 f_prod  =  sum(WP.prod(pairIx_prod).alloc(: , pairIx_inj)); % [m^3 /s]


end
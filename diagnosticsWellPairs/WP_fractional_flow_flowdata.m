function [fw,qw,qo] = WP_fractional_flow_flowdata(well_pair,wellSols,WP,D)

    
     pairIx_inj  =  WP.pairIx(well_pair,1) ; % Injetor index in WP
     pairIx_prod =  WP.pairIx(well_pair,2) ; % Producer index in WP
     
     Ix_inj  = D.inj(  WP.pairIx(well_pair,1) ); % Injetor index in WellSols
     Ix_prod = D.prod( WP.pairIx(well_pair,2) ); % Pruducer index in WellSols


qw = wellSols{1}(Ix_prod).qWs;
qo = wellSols{1}(Ix_prod).qOs;

    if qw == 0  
        fw = 0;
    else    
    fw =  1.0 /...
            ( 1.0 + (qo/qw));
    end

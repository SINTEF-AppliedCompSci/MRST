function s = setupSimCompVe(Gt, rock2D, varargin)
% Wrapper for setupSimComp to produce a vertically-integrated version.
    
    % Compute vertially-integrated transmissibilities
    rock_tmp = rock2D;
    rock_tmp.perm = rock2D.perm .* Gt.cells.H;
    T = computeTrans(Gt,rock_tmp);
    cf = Gt.cells.faces(:,1); 
    nf = Gt.faces.num;
    T  = 1 ./ accumarray(cf, 1./T, [nf, 1]);
    T_all = T;
    
    % Computing vertically-integrated pore-volume
    pv = poreVolume(Gt,rock2D);
    pv = pv.*Gt.cells.H;    
    
    % Calling setupSimComp, substituting the pore-volume and transmissiblity
    % values it computes with our vertically-integrated quantities.
    s = setupSimComp(Gt, rock2D, 'porv', pv,'trans',T);
    s.T_all=T_all;
end
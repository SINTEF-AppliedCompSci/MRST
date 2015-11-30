function [ grdecl, G, Gt, rock2D ] = makeJohansenFullField(  )
% Load NPD's full-field model of Johansen

    try
       jdir = fullfile(mrstPath('co2lab'), 'data', 'johansen');
       sector = 'FULLFIELD_IMAXJMAX';
       sector = fullfile(jdir, sector);
       grdecl = readGRDECL([sector '.GRDECL']);
    catch me
       disp(' -> Download data from: http://www.sintef.no/Projectweb/MatMoRA/')
       disp(['    Putting data in ', jdir]);
       unzip('http://www.sintef.no/project/MatMoRA/Johansen/FULLFIELD_Eclipse.zip', jdir);
       grdecl = readGRDECL([sector '.GRDECL']);
    end

    % Extract the part that represents the Johansen formation
    grdecl = cutGrdecl(grdecl,[1 grdecl.cartDims(1); 1 grdecl.cartDims(2);  10 14]);
    G  = processgrid(grdecl);
    G = mcomputeGeometry(G);
    [Gt, G] = topSurfaceGrid(G); % @@
    
    %% rock properties
    avgperm = mean(grdecl.PERMX( grdecl.PERMX~=0 ));
    avgporo = mean(grdecl.PORO( grdecl.PORO~=0 ));
    
    % should be able to get full heterogeneous data @@
    rock2D.perm = avgperm*ones(Gt.cells.num,1);
    rock2D.poro = avgporo*ones(Gt.cells.num,1);

end


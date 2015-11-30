function [ grdecl, G, Gt, rock2D ] = makeJohansenSector()
% Construct NPD's sector model of Johansen


    mrstModule add co2lab;
    moduleCheck('libgeometry','opm_gridprocessing','coarsegrid','matlab_bgl');

    %% Load NPD data: sector model
    try
       jdir = fullfile(mrstPath('co2lab'), 'data', 'johansen');
       sector = 'NPD5';
       sector = fullfile(jdir, sector);
       grdecl = readGRDECL([sector '.grdecl']);
       
    catch me
       disp(' -> Download data from: http://www.sintef.no/Projectweb/MatMoRA/')
       disp(['    Putting data in ', jdir]);
       unzip('http://www.sintef.no/project/MatMoRA/Johansen/NPD5.zip', jdir);
       grdecl = readGRDECL([sector '.grdecl']);
    end

    % Load permeability and porosity
    K = reshape(load([sector, '_Permeability.txt'])', [], 1);
    p = reshape(load([sector, '_Porosity.txt'])',     [], 1);
    grdecl.PERMX = K;
    grdecl.PORO = p;
    
    % Extract the part that represents the Johansen formation
    % (also cuts out rock properties)
    grdecl = cutGrdecl(grdecl, ...
        [1 grdecl.cartDims(1); 1 grdecl.cartDims(2);  6 11]);
    G  = processgrid(grdecl);
    G  = mcomputeGeometry(G);
    [Gt, G] = topSurfaceGrid(G); % @@ G has changed
    
    %% form rock properties
    % 2D rock structure
    rock2D.perm = grdecl.PERMX(Gt.cells.indexMap).*milli*darcy;
    rock2D.poro = grdecl.PORO(Gt.cells.indexMap);


end


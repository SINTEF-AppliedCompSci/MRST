function G=meshHFsystem(physdim,frac_centroid_s, varargin)
    %This will generate structured grid by putting hydraulic fractures on the grid faces
    %Author Harun Rashid
    opt = struct('aperture', 0.02,...
                 'numStages', [],...
                 'fracSpacing', 10,...
                 'fracHalfLength', 100,...
                 'fracHeight', 60,...
                 'nxRefine',10,...
                 'numFracsPerStage',[],...
                 'clusterSpacing',100,...
                 'ny',11, ...
                 'nxRefine_small',2,...
                 'nz',11 );

    opt = merge_options(opt, varargin{:});
    i=1;
   
    facesXcoords=[];
    pt1=0;
    for k=1:opt.numStages
        if k <=opt.numStages
            ptX1 = linspace(pt1, frac_centroid_s(i), opt.nxRefine+1);
            facesXcoords=[facesXcoords, ptX1];
            pt1=facesXcoords(end);
        end
        for j=1:opt.numFracsPerStage-1
            ptX1 = linspace(pt1, frac_centroid_s(j+i), opt.nxRefine_small+1);
            facesXcoords=[facesXcoords, ptX1];
            pt1=facesXcoords(end);
        end
        i=i+opt.numFracsPerStage;
        if k ==opt.numStages
            ptX1 = linspace(pt1, physdim(1), opt.nxRefine+1);
            facesXcoords=[facesXcoords, ptX1];
            pt1=facesXcoords(end);
        end
    end

    facesXcoords=unique(facesXcoords);
    
        
    facesYcoords = linspace(0, physdim(2),opt.ny);
    
    facesZcoords = linspace(0,physdim(3),opt.nz);
    
    G = tensorGrid(facesXcoords, facesYcoords, facesZcoords);
        

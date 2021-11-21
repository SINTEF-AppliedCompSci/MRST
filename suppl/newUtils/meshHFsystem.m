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
    aa=frac_centroid_s(1)./(opt.fracSpacing/(opt.nxRefine-1));
    bb=(physdim(1)-frac_centroid_s(end))./(opt.fracSpacing/(opt.nxRefine-1));
    
    % This part will create grid from the left boundary to the first
    % fracture
    ptX = linspace(0, frac_centroid_s(1), aa);
    wellCellcentroidX = frac_centroid_s(1);
    facesXcoords = wellCellcentroidX + -fliplr(ptX);
    i=1;
    n=1;
    for k=1:opt.numStages
        for j=1:opt.numFracsPerStage-1
            ptX1 = linspace(0, frac_centroid_s(n+j)-frac_centroid_s(n+j-1), opt.nxRefine_small);
            wellCellcentroidX=frac_centroid_s(n+j-1);
            facesXcoords1 = wellCellcentroidX + ptX1;
            facesXcoords=[facesXcoords, facesXcoords1];
        end
        if k <=opt.numStages-1
            ptX1 = linspace(0, frac_centroid_s(i+opt.numFracsPerStage)-frac_centroid_s(i+opt.numFracsPerStage-1), opt.nxRefine);
            wellCellcentroidX=frac_centroid_s(i+opt.numFracsPerStage-1);
            facesXcoords2 = wellCellcentroidX + ptX1;
            facesXcoords=[facesXcoords, facesXcoords2];
            i=i+opt.numFracsPerStage;
            n=n+opt.numFracsPerStage;
        end
    end
    % This part will create grid from the last fracture to the domain
    % boundary
    ptX = linspace(0, physdim(1)-frac_centroid_s(end), bb);
    wellCellcentroidX = frac_centroid_s(end);
    facesXcoords3 = wellCellcentroidX + ptX;
    facesXcoords=[facesXcoords, facesXcoords3];
    facesXcoords=unique(facesXcoords);
    
        
    facesYcoords = linspace(0, physdim(2),opt.ny);
    
    facesZcoords = linspace(0,physdim(3),opt.nz);
    
    G = tensorGrid(facesXcoords, facesYcoords, facesZcoords);
        

function [G,frac_ids,well_ids]=explicitStencil3Fracs(physdim, varargin)
    %{

    %}
    opt = struct('aperture', 0.02,...
                 'numStages'        ,  1,                      ...
                 'fracSpacing'      , 100,                      ...
                 'fracHalfLength'   , 100,                      ...
                 'fracHeight', 60,                             ...
                 'nxRefine',16,...
                 'nxRefineSmall',8,...
                 'clusterSpacing',65*ft,...
                 'tipNY',5,...
                 'ny',11, ...
                 'nz',11, 'gridType','geomspacing');

    opt = merge_options(opt, varargin{:});
    
    switch(lower(opt.gridType))
        case 'geomspacing'
            % left and right of middle frac 
            ptX = geomspace(opt.aperture/2, opt.clusterSpacing/2, opt.nxRefineSmall,opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords1 = wellCellcentroidX + [-fliplr(ptX), ptX];
            % right of left frac
            wellCellcentroidX = opt.fracSpacing/2 - opt.clusterSpacing;
            facesXcoords = wellCellcentroidX + ptX(1:end-1);
            facesXcoords = [facesXcoords, facesXcoords1];
            % left of right frac
            wellCellcentroidX = opt.fracSpacing/2 + opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX - fliplr(ptX(1:end-1));
            facesXcoords = [facesXcoords, facesXcoords1];
            % large space to the right of the right frac
            ptX = geomspace(opt.aperture/2, opt.fracSpacing/2-opt.clusterSpacing, opt.nxRefine,opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2 + opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX + ptX;
            facesXcoords = [facesXcoords, facesXcoords1];
            % large space to the left of the left frac
            wellCellcentroidX = opt.fracSpacing/2 - opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX - fliplr(ptX);
            facesXcoords = [facesXcoords1, facesXcoords];
        case 'geomspacing3frac' 
           
            ob=(opt.physdim1-(opt.fracSpacing*3))/2;
            ptX = linspace(opt.aperture/2, ob, opt.nxRefine);
            wellCellcentroidX = (opt.physdim1/2)-(opt.fracSpacing/2)*3;
            facesXcoords1 = wellCellcentroidX + [-fliplr(ptX)];
            
            ptX = geomspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine,opt.aperture);
            wellCellcentroidX = (opt.physdim1/2)-opt.fracSpacing;
            facesXcoords2 = wellCellcentroidX + [-fliplr(ptX) ,ptX];
            facesXcoords3 = [facesXcoords1, facesXcoords2];
            
            ptX = geomspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine,opt.aperture);
            wellCellcentroidX = opt.physdim1/2;
            facesXcoords4 = wellCellcentroidX + [-fliplr(ptX) ,ptX];
            facesXcoords5 = [facesXcoords3, facesXcoords4];
            
            ptX = geomspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine,opt.aperture);
            wellCellcentroidX = opt.physdim1/2+opt.fracSpacing;
            facesXcoords6 = wellCellcentroidX + [-fliplr(ptX) ,ptX];
            facesXcoords7 = [facesXcoords5, facesXcoords6];
            
            
            ptX = linspace(opt.aperture/2, ob, opt.nxRefine);
            wellCellcentroidX = (opt.physdim1/2)+(opt.fracSpacing/2)*3;
            facesXcoords8 = wellCellcentroidX + ptX;
            facesXcoords = [facesXcoords7, facesXcoords8];
            
%             ptX = linspace(opt.aperture/2, ob, opt.nxRefine);
%             wellCellcentroidX = (opt.physdim1/2)-(opt.fracSpacing/2)*3;
%             facesXcoords5 = wellCellcentroidX + [-fliplr(ptX)];
%             facesXcoords = [facesXcoords4, facesXcoords5];
%             
%             facesXcoords = [facesXcoords4, facesXcoords5];

            %             
        case 'logspacing'
            ptX = logspace(log10(opt.aperture/2),log10(opt.fracSpacing/2),opt.nxRefine);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'linspacing'
            ptX = linspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'cartesian'
%             facesXcoords = linspace(0, opt.fracSpacing, 2*opt.nxRefine);

            % left and right of middle frac 
            ptX = geomspace(opt.aperture/2, opt.clusterSpacing/2, opt.nxRefineSmall,opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords1 = wellCellcentroidX + [-fliplr(ptX), ptX];
            % right of left frac
            ptX = geomspace(2*opt.aperture/2, opt.clusterSpacing/2, opt.nxRefineSmall,2*opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2 - opt.clusterSpacing;
            facesXcoords = wellCellcentroidX + ptX(1:end-1);
            facesXcoords = [facesXcoords, facesXcoords1];
            % left of right frac
            wellCellcentroidX = opt.fracSpacing/2 + opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX - fliplr(ptX(1:end-1));
            facesXcoords = [facesXcoords, facesXcoords1];
            % large space to the right of the right frac
            ptX = geomspace(2*opt.aperture/2, opt.fracSpacing/2-opt.clusterSpacing, opt.nxRefine,2*opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2 + opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX + ptX;
            facesXcoords = [facesXcoords, facesXcoords1];
            % large space to the left of the left frac
            wellCellcentroidX = opt.fracSpacing/2 - opt.clusterSpacing;
            facesXcoords1 = wellCellcentroidX - fliplr(ptX);
            facesXcoords = [facesXcoords1, facesXcoords];
        otherwise
            error('Unknown case')
    end

    frontTipY = physdim(2)/2 - opt.fracHalfLength;
    backTipY = physdim(2)/2 + opt.fracHalfLength;
    pt2FrontTip = linspace(0,frontTipY,opt.tipNY);
    ptFront2BackTip = linspace(frontTipY,backTipY,opt.ny);
    ptBackTip2boundary = linspace(backTipY,physdim(2),opt.tipNY);
    facesYcoords = [pt2FrontTip(1:end-1), ptFront2BackTip(1:end-1), ptBackTip2boundary];

    facesZcoords = linspace(0,physdim(3),opt.nz);

    G = tensorGrid(facesXcoords, facesYcoords, facesZcoords);
%     plotGrid(G)

    G2D = tensorGrid(facesXcoords, facesYcoords);
    plotGrid(G2D)

    %OMO edit:
    %code snippet below spits out grid data for CMG simulator
    DX=diff(facesXcoords)';
    DY=diff(facesYcoords)';
    DZ=diff(facesZcoords)';
%     Folder = cd; Folder = fullfile(Folder, '..');
%     fid = fopen(fullfile(Folder, '\mat_data\DI_DJ_DK.txt'), 'w'); 
%     fprintf(fid, 'DI IVAR     ** %d float\n', numel(DX));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DX(:),ft));fprintf(fid, '\n');
%     fprintf(fid, 'DJ JVAR     ** %d float\n', numel(DY));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DY(:),ft)); fprintf(fid, '\n');
%     fprintf(fid, 'DK KVAR     ** %d float\n', numel(DZ));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DZ(:),ft)); fprintf(fid, '\n');
%     fclose(fid);

    %% Find out the cell index for all fracture plane and wells
    switch(lower(opt.gridType))
        case 'geomspacing'
            Frac_I = [opt.nxRefine,opt.nxRefine+2*opt.nxRefineSmall-1,opt.nxRefine+4*opt.nxRefineSmall-2];
        case 'cartesian'
            Frac_I = [opt.nxRefine,opt.nxRefine+2*opt.nxRefineSmall-1,opt.nxRefine+4*opt.nxRefineSmall-2];
%             Frac_I = [opt.nxRefine-ceil((66*ft)/(physdim(1)/(2*opt.nxRefine+1))),opt.nxRefine,opt.nxRefine+ceil((66*ft)/(physdim(1)/(2*opt.nxRefine+1)))];
        otherwise
            Frac_I = opt.nxRefine;
    end
    Frac_J = (opt.tipNY):(opt.tipNY-1+opt.ny-1);
    Frac_K = 1:opt.nz-1;
    NumFracs= numel(Frac_I);
    frac_ids=[];
    well_ids=[];
    for fi=1:NumFracs
        [II,JJ,KK]=meshgrid(Frac_I(fi),Frac_J,Frac_K);
        FracCellIds = sub2ind(G.cartDims, II,JJ,KK);
        frac_ids=[frac_ids; FracCellIds(:)'];

        centerId=ceil(numel(frac_ids(fi,:))/2);
        well_ids=[well_ids; frac_ids(fi,centerId)]; %[i,j,k]=ind2sub(G.cartDims, well_ids(2))
    end
    %TO DO
    well_ids = well_ids(2);

end

function [G,frac_ids,well_ids]=explicitStencilSuppl(physdim, varargin)
    %{

    %}
    opt = struct('aperture', 0.02,...
                 'numStages'        ,  1,                      ...
                 'fracSpacing'      , 100,                      ...
                 'fracHalfLength'   , 100,                      ...
                 'fracHeight', 60,                             ...
                 'nxRefine',10,...
                 'tipNY',5,...
                 'physdim1',200,...
                 'ny',11, ...
                 'nz',11, 'gridType','geomspacing');

    opt = merge_options(opt, varargin{:});
    
    switch(lower(opt.gridType))
        case 'geomspacing'
            ptX = geomspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine,opt.aperture);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        
        case 'geomspacing3frac' %Only for three fracture and can be extended for multisatge fracture in explicit model(HR edit)
           
            ob=(opt.physdim1-(opt.fracSpacing*3))/2;
            ptX = geomspace(opt.aperture/2, ob, opt.nxRefine,opt.aperture);
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
            
            
            ptX = geomspace(opt.aperture/2, ob, opt.nxRefine,opt.aperture);
            wellCellcentroidX = (opt.physdim1/2)+(opt.fracSpacing/2)*3;
            facesXcoords8 = wellCellcentroidX + ptX;
            facesXcoords = [facesXcoords7, facesXcoords8];
        case 'logspacing'
            ptX = logspace(log10(opt.aperture/2),log10(opt.fracSpacing/2),opt.nxRefine);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'linspacing'
            ptX = linspace(opt.aperture/2, opt.fracSpacing/2, opt.nxRefine);
            wellCellcentroidX = opt.fracSpacing/2;
            facesXcoords = wellCellcentroidX + [-fliplr(ptX) ,ptX];
        case 'cartesian'
            facesXcoords = linspace(0, opt.fracSpacing, 2*opt.nxRefine);             
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
    
    %Hassan Edit:
    G = tensorGrid(facesXcoords, facesYcoords, facesZcoords);
    G.XGrid = facesXcoords; G.YGrid = facesYcoords; G.ZGrid = facesZcoords;
%     plotGrid(G)

    %OMO edit:
    %code snippet below spits out grid data for CMG simulator
%     DX=diff(facesXcoords)';
%     DY=diff(facesYcoords)';
%     DZ=diff(facesZcoords)';
%     fid = fopen('C:\Users\hrashi1\Documents\Research_Codes\mrst-2018b\modules\pEDFM\mat_data\DI_DJ_DK.txt', 'w'); 
%     fprintf(fid, 'DI IVAR     ** %d float\n', numel(DX));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DX(:),ft));fprintf(fid, '\n');
%     fprintf(fid, 'DJ JVAR     ** %d float\n', numel(DY));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DY(:),ft)); fprintf(fid, '\n');
%     fprintf(fid, 'DK KVAR     ** %d float\n', numel(DZ));
%     fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DZ(:),ft)); fprintf(fid, '\n');
%     fclose(fid);
% 
%     %% Find out the cell index for all fracture plane and wells
    %Frac_I = [opt.nxRefine-2,opt.nxRefine,opt.nxRefine+2];%for 3 frac
    Frac_I = [opt.nxRefine];
    Frac_J = (opt.tipNY):(opt.tipNY-1+opt.ny-1);
    Frac_K = 1:opt.nz-1;
    NumFracs= numel(Frac_I);
    frac_ids=[];
    well_ids=[];
    G.Frac_I = Frac_I;
    for fi=1:NumFracs
        [II,JJ,KK]=meshgrid(Frac_I(fi),Frac_J,Frac_K);
        FracCellIds = sub2ind(G.cartDims, II,JJ,KK);
        frac_ids=[frac_ids; FracCellIds(:)'];

        centerId=ceil(numel(frac_ids(fi,:))/2);
        well_ids=[well_ids; frac_ids(fi,centerId)];
    end
   % well_ids = well_ids(2);
    %% Find out the cell index for all fracture planeand wells: This case only valid for three fracture cell 
%     Frac_I = [2*opt.nxRefine,4*opt.nxRefine,6*opt.nxRefine]; %edit HR
%     Frac_J = (opt.tipNY):(opt.tipNY-1+opt.ny-1);
%     Frac_K = 1:opt.nz-1;
%     NumFracs= numel(Frac_I);
%     frac_ids=[];
%     well_ids=[];
%     for fi=1:NumFracs
%         [II,JJ,KK]=meshgrid(Frac_I(fi),Frac_J,Frac_K);
%         FracCellIds = sub2ind(G.cartDims, II,JJ,KK);
%         frac_ids=[frac_ids; FracCellIds(:)'];
% 
%         centerId=ceil(numel(frac_ids(fi,:))/2);
%         well_ids=[well_ids; frac_ids(fi,centerId)];
%     end
%%
% DX=diff(facesXcoords)';
% DY=diff(facesYcoords)';
% DZ=diff(facesZcoords)';
% fid = fopen('C:\Users\ahass16\Desktop\Research Work\mrst-2020a\modules\shale_chapter\examples\Diffusion_CMG_Validation\DI_DJ_DK_Coarse.txt', 'w'); 
% fprintf(fid, 'DI IVAR     ** %d float\n', numel(DX));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DX(:),ft));fprintf(fid, '\n');
% fprintf(fid, 'DJ JVAR     ** %d float\n', numel(DY));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DY(:),ft)); fprintf(fid, '\n');
% fprintf(fid, 'DK KVAR     ** %d float\n', numel(DZ));
% fprintf(fid, '%0.10f %0.10f %0.10f %0.10f %0.10f %0.10f\n', convertTo(DZ(:),ft)); fprintf(fid, '\n');
% fclose(fid);
end
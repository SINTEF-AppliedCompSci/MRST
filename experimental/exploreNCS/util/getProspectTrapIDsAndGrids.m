function [ Grids, pts ] = getProspectTrapIDsAndGrids( G_st, Gt_st, ta_st )
% location of point corresponding to "prospects" or "structural closures":

    % Snohvit
    snohvit_x      = 9.202e5;   pts.Xpt_gsf = snohvit_x;
    snohvit_y      = 7.988e6;   pts.Ypt_gsf = snohvit_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [snohvit_x, snohvit_y]);
    [~,i] = min(sum(dv.^2, 2));
    snohvit_cellIndex = i; % or Gt.cells.indexMap(i);


    % Albatross includes prospects E and F
    albatross_x = 9.402e5;  pts.Xpt_alba = albatross_x;
    albatross_y = 7.974e6;  pts.Ypt_alba = albatross_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [albatross_x, albatross_y]);
    [~,i] = min(sum(dv.^2, 2));
    albatross_cellIndex = i; % or Gt.cells.indexMap(i);



    % Askeladd
    askeladd_x = 9.077e5;   pts.Xpt_aske = askeladd_x;
    askeladd_y = 7.957e6;   pts.Ypt_aske = askeladd_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [askeladd_x, askeladd_y]);
    [~,i] = min(sum(dv.^2, 2));
    askeladd_cellIndex = i; % or Gt.cells.indexMap(i);


    % Other prospects
    prospectC_x = 9.577e5;  pts.Xpt_prosC = prospectC_x;
    prospectC_y = 8.012e6;  pts.Ypt_prosC = prospectC_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [prospectC_x, prospectC_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectC_cellIndex = i; % or Gt.cells.indexMap(i);

    prospectD_x = 9.502e5;  pts.Xpt_prosD = prospectD_x;
    prospectD_y = 8.007e6;  pts.Ypt_prosD = prospectD_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [prospectD_x, prospectD_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectD_cellIndex = i; % or Gt.cells.indexMap(i);

    prospectG_x = 9.217e5;  pts.Xpt_prosG = prospectG_x;
    prospectG_y = 7.942e6;  pts.Ypt_prosG = prospectG_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [prospectG_x, prospectG_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectG_cellIndex = i; % or Gt.cells.indexMap(i);

    prospectH_x = 8.992e5;  pts.Xpt_prosH = prospectH_x;
    prospectH_y = 7.916e6;  pts.Ypt_prosH = prospectH_y;
    dv = bsxfun(@minus, Gt_st.cells.centroids(:,1:2), [prospectH_x, prospectH_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectH_cellIndex = i; % or Gt.cells.indexMap(i);





    % IDs of prosects or structural closures:
    snohvitTrapID   = ta_st.traps(snohvit_cellIndex);
    prospectCTrapID = ta_st.traps(prospectC_cellIndex);
    prospectDTrapID = ta_st.traps(prospectD_cellIndex);
    %prospectETrapID = ta_st.traps(prospectE_cellIndex);
    %prospectFTrapID = ta_st.traps(prospectF_cellIndex);
    prospectGTrapID = ta_st.traps(prospectG_cellIndex);
    prospectHTrapID = ta_st.traps(prospectH_cellIndex);
    albatrossTrapID = ta_st.traps(albatross_cellIndex);
    askeladdTrapID  = ta_st.traps(askeladd_cellIndex);

    % cut away rest of grid except cells included in this structural trap
    G_gsf       = removeCells(G_st, ta_st.traps ~= snohvitTrapID);      G_gsf.name = 'Greater Snohvit';
    G_prospectC = removeCells(G_st, ta_st.traps ~= prospectCTrapID);    G_prospectC.name = 'Prospect C';
    G_prospectD = removeCells(G_st, ta_st.traps ~= prospectDTrapID);    G_prospectD.name = 'Prospect D';
    %G_prospectE = removeCells(G_st, ta_st.traps ~= prospectETrapID);
    %G_prospectF = removeCells(G_st, ta_st.traps ~= prospectFTrapID);
    G_prospectG = removeCells(G_st, ta_st.traps ~= prospectGTrapID);    G_prospectG.name = 'Prospect G';
    G_prospectH = removeCells(G_st, ta_st.traps ~= prospectHTrapID);    G_prospectH.name = 'Prospect H';
    [G_albatross, cellmap_albatross] = removeCells(G_st, ta_st.traps ~= albatrossTrapID);    G_albatross.name = 'Greater Albatross';
    G_askeladd  = removeCells(G_st, ta_st.traps ~= askeladdTrapID);     G_askeladd.name  = 'Greater Askeladd';


    % Get prospect E and F from trapping analysis performed on G_albatross
    Gt_albatross = topSurfaceGrid(G_albatross);
    ta_albatross = trapAnalysis(Gt_albatross, false); % trap info corresponds to Gt_albatross grid

    prospectE_x = 9.492e5;  pts.Xpt_prosE = prospectE_x;
    prospectE_y = 7.983e6;  pts.Ypt_prosE = prospectE_y;
    dv = bsxfun(@minus, Gt_albatross.cells.centroids(:,1:2), [prospectE_x, prospectE_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectE_cellIndex = i; % or Gt.cells.indexMap(i);

    prospectF_x = 9.527e5;  pts.Xpt_prosF = prospectF_x;
    prospectF_y = 7.974e6;  pts.Ypt_prosF = prospectF_y; 
    dv = bsxfun(@minus, Gt_albatross.cells.centroids(:,1:2), [prospectF_x, prospectF_y]);
    [~,i] = min(sum(dv.^2, 2));
    prospectF_cellIndex = i; % or Gt.cells.indexMap(i);

    prospectETrapID = ta_albatross.traps(prospectE_cellIndex); % trapID in Gt_albatross
    prospectFTrapID = ta_albatross.traps(prospectF_cellIndex); % trapID in Gt_albatross

    [G_prospectE, cellmapE] = removeCells(G_albatross, ta_albatross.traps ~= prospectETrapID);  G_prospectE.name = 'Prospect E';
    [G_prospectF, cellmapF] = removeCells(G_albatross, ta_albatross.traps ~= prospectFTrapID);  G_prospectF.name = 'Prospect F';

    
    % cellmap_albatross is cinx of Gt_st
    % cellmapE and cellmapF is cinx of G_albatross
    

    Grids.G_gsf                 = G_gsf;
    Grids.G_gsf.trapID          = snohvitTrapID;
    Grids.G_gsf.trapcells       = getcinx(Gt_st, ta_st, snohvitTrapID);
    % or use [G_gsf, cellmapGSF] = removeCells(...
    % where trapcells = cellmapGSF;
    
    Grids.G_albatross           = G_albatross;
    Grids.G_albatross.trapID    = albatrossTrapID;
    Grids.G_albatross.trapcells = getcinx(Gt_st, ta_st, albatrossTrapID);
    
    Grids.G_askeladd            = G_askeladd;
    Grids.G_askeladd.trapID     = askeladdTrapID;
    Grids.G_askeladd.trapcells  = getcinx(Gt_st, ta_st, askeladdTrapID);
    
    Grids.G_prospectC           = G_prospectC;
    Grids.G_prospectC.trapID    = prospectCTrapID;
    Grids.G_prospectC.trapcells = getcinx(Gt_st, ta_st, prospectCTrapID);
    
    Grids.G_prospectD           = G_prospectD;
    Grids.G_prospectD.trapID    = prospectDTrapID;
    Grids.G_prospectD.trapcells = getcinx(Gt_st, ta_st, prospectDTrapID);
    
    Grids.G_prospectE           = G_prospectE;
    Grids.G_prospectE.trapID    = prospectETrapID;
    Grids.G_prospectE.trapcells = cellmap_albatross(cellmapE);
    
    Grids.G_prospectF           = G_prospectF;
    Grids.G_prospectF.trapID    = prospectFTrapID;
    Grids.G_prospectF.trapcells = cellmap_albatross(cellmapF);
    
    Grids.G_prospectG           = G_prospectG;
    Grids.G_prospectG.trapID    = prospectGTrapID;
    Grids.G_prospectG.trapcells = getcinx(Gt_st, ta_st, prospectGTrapID);
    
    Grids.G_prospectH           = G_prospectH;
    Grids.G_prospectH.trapID    = prospectHTrapID;
    Grids.G_prospectH.trapcells = getcinx(Gt_st, ta_st, prospectHTrapID);
    
    


end

function trapcinx = getcinx(Gt, ta, trapID)
% get cell indexes of grid Gt that make up trap given by trapID

    trapcinx                      = zeros(Gt.cells.num,1);
    trapcinx(ta.traps == trapID)  = trapID;
    trapcinx                      = find(trapcinx);

end


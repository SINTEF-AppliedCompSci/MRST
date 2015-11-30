function [ Gt_faults, faultFaces, faultCells, s ] = implementFaults(Gt, rock2D, varargin)
% Detects potential fault lines based on: a) elevation gradients above a
% tolerance value, or b) specified coordinates
%
% SYNOPSIS
%   implementFaults(Gt, rock2D)
%   implementFaults(Gt, rock2D, 'gradTol',100)
%   implementFaults(Gt, [], 'faultCoords',[Xcoords(:), Ycoords(:)])
%
% DESCRIPTION
%   The face index belonging to fault line are returned in faultFaces.
%   Cells along fault lines can be removed from grid.
%   Also, faces belonging to fault lines can be assigned a new
%   transmissibility (either lower or higher than original, depending on
%   whether fault is sealing or conducting).

% SEE ALSO:
%    makeInternalBoundary.m
%    topSurfaceGrid(...,'add_faults',true)
%


moduleCheck('ad-core')

opt.gradTol     = 100; % elevation diff between two adjacent cells divided by cell spacing.
opt.faultCoords = [];   

opt = merge_options(opt, varargin{:});

if isempty(opt.faultCoords)
    % Faults are implemented by approach (a)
        
    % Use gradient operator, Grad(), to compute elevation gradients of
    % top-surface grid. Gradients are computed on the interfaces via a
    % first order finite difference approximation using the values of the
    % cells connected to the face. The transmissibility values for all
    % faultFaces are assigned as zero.
    s                       = setupOperatorsTPFA(Gt, rock2D);
    gradZ                   = s.Grad(Gt.cells.z);
    s.T( abs(gradZ) > opt.gradTol ) = 0;
    tmp                     = zeros(Gt.faces.num,1);
    tmp( s.internalConn )   = s.T; % @@
    s.T_all                 = tmp;
    
    %faultfaces=find(s.T_all==0); % @@ this sets even external bdry trans
    %to zero. Avoid this.
    %plotFaces(Gt,faultfaces)     % @@

    % All faces that have an elevation gradient above a pre-set tolerance
    % (or that were just assigned zero transmissibilities) are assigned as
    % faultFaces.
    finx                    = zeros(size(s.T));
    finx( s.T==0 )          = 1; % @@
    finx_all                = zeros(Gt.faces.num,1);
    finx_all( s.internalConn ) = finx;
    faultFaces              = find( finx_all );
    
    figure; set(gcf,'Position',[1 1 1500 500],'name',['Faults detected using grad(Gt.cells.z)>', num2str(opt.gradTol)])
    subplot(1,2,1)
    title('Detected fault faces')
    plotFaces(Gt, faultFaces,'r', 'EdgeColor','r','LineWidth',3)
    plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1)
    view(3)
    

    % All cells that are neighbors to the faultFaces are assigned as
    % faultCells. All faultCells are removed from the grid. After removing
    % cells, need to clear variable Gt.cells.sortedCellNodes
    faultCells  = [Gt.faces.neighbors( faultFaces, 1 ); ...
                   Gt.faces.neighbors( faultFaces, 2)];
    faultCells  = unique( faultCells );
    faultCells  = faultCells( find(faultCells~=0) );
    
    subplot(1,2,2)
    title('Cells on either side of fault faces')
    plotCellData(Gt, Gt.cells.z, faultCells, 'EdgeAlpha',0.1)
    plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1)
    view(3)
    
    % attempt to remove cells from 2D grid - didn't work properly @@
    %Gt_faults       = removeCells(Gt, faultCells);
    %Gt_faults.cells = rmfield(Gt_faults.cells, 'sortedCellNodes');
    
    % remove cells from 3D grid, so top surface grid is derived from the
    % cut 3D grid 
    %G_faults              = removeCells(Gt.parent, faultCells);
    %[Gt_faults, G_faults] = topSurfaceGrid(mcomputeGeometry(G_faults)); % @@ time-consuming
    Gt_faults = Gt; % @@
    
    % Visualization
    %figure;
    %plotCellData(Gt_faults, Gt_faults.cells.volumes)

    
    % perform trapAnalysis on update grid Gt
    %ta_faults = trapAnalysis(Gt_faults, false);
    
    % ---------------------------------------------------------------------

    % Better approach is to alter face transmissibilities? 
    
elseif ~isempty(opt.faultCoords)
    % Faults are implemented by approach (b)

end

end

% -------------------------------------------------------------------------

function getTranslatedFaultCoords()
% Fault info is read from text file (containing line coordinates X,Y), and
% then translated and scaled such that coordinates are some-what consistent
% with Sto.


    %load BarentsSeaFaultCoords.mat
    load HammerfestBasinFaultCoords.mat

    % load Sto formation
    [Gt, rock2D, petrodata] = getFormationTopGrid('Stofm',1);

    % place fault lines ontop of Sto grid
    % requires adjustments to fault coords:
    a = 4.5048e5;
    b = 7.8702e6;
    A = 8.8545e5;
    B = 7.904e6;
    c = 5.3919e5;
    d = 7.9080e6;
    C = 9.65e5;
    D = 7.958e6;

    % translation st bottom left corner of fault coords match with Sto
    % grid's bottom left corner
    dx = A - a;
    dy = B - b;

    X = X + dx;
    Y = Y + dy;

    % now scaling
    W = C - A;
    w = c - a;
    alpha = W/w;

    H = D - B;
    h = d - b;
    beta = H/h;

    % (xnew, ynew) = (xc + alpha(xold - xc), yc + beta(yold - yc)) (xc,yc)
    % is the reference coordinate where old and new coordinates match.
    xc = A;
    yc = B;

    X = xc + alpha*(X - xc);
    Y = yc + beta*(Y - yc);

    % plot
    figure;
    plotGrid(Gt, 'EdgeAlpha',0.1)
    hold on
    line(X,Y)
    axis equal tight


end


function [FS_ieaghg, FS_original, FS_inhouse] = createGridInterpolants(Gt_ieaghg, Gt_original, Gt_inhouse)
%% Prepare and return grid interpolants
%  (to facilitate comparison, since grids are of different resolution) A
%  scattered-point interpolant is useful when the domain of the formation
%  has an irregular shape.

% Output: F


% See also:
%   resTiltUtsira.m

Gt1 = Gt_ieaghg;
Gt2 = Gt_original;
Gt3 = Gt_inhouse;



    %% Interpolant for Sleipner IEAGHG model grid
    FS_ieaghg = TriScatteredInterp(Gt1.nodes.coords(:,1), ...
                            Gt1.nodes.coords(:,2), ...
                            Gt1.nodes.z);
    
    
    %% Interpolant for Sleipner ORIGINAL model grid
    FS_original = TriScatteredInterp(Gt2.nodes.coords(:,1), ...
                            Gt2.nodes.coords(:,2), ...
                            Gt2.nodes.z);
    
    
    %% Interpolant for Sleipner INHOUSE model grid
    FS_inhouse = TriScatteredInterp(Gt3.nodes.coords(:,1), ...
                            Gt3.nodes.coords(:,2), ...
                            Gt3.nodes.z);


    %% from resTiltUtsira.m
%     %% Interpolant for Sleipner grid
%     [nx, ny]  = deal(GtS.cartDims(1)+1, GtS.cartDims(2)+1);
% 
%     X         = reshape(GtS.nodes.coords(:,1), nx, ny); 
%     Y         = reshape(GtS.nodes.coords(:,2), nx, ny);
%     Z         = reshape(GtS.nodes.z          , nx, ny);
% 
%     FS        = TriScatteredInterp(X(:), Y(:), Z(:));
%     
%     %% Interpolant for Utsira grid
%     FU = TriScatteredInterp(GtU.nodes.coords(:,1), ...
%                             GtU.nodes.coords(:,2), ...
%                             GtU.nodes.z);
end
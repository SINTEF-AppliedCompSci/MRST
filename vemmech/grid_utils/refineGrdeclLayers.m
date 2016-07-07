function ngrdecl = refineGrdeclLayers(grdecl, layers, ref)
%
%
% SYNOPSIS:
%   function ngrdecl = refineGrdeclLayers(grdecl, layers, ref)
%
% DESCRIPTION: Refine grid in vertical direction
%
% PARAMETERS:
%   grdecl   - Grid in Eclipse format
%   layers   - layers to be refined
%   ref      - refinement factor
%
% RETURNS:
%   ngrdecl - refined grid in eclipse format.
%
% EXAMPLE:
%
% SEE ALSO:
%

    [xyz, zcorn] = grdeclXYZ(grdecl);
    zcorn1 = zcorn(:, :, 1:2*(layers(1)-1));
    zcorn2 = zcorn(:, :, 2*(layers(2)-1)+1:end);
    cgrdecl = cutGrdecl(grdecl, [1 grdecl.cartDims(1);1 grdecl.cartDims(2);layers]);
    cgrdecl = refineGrdecl(cgrdecl, [1 1 ref]);
    [cxyz, czcorn] = grdeclXYZ(cgrdecl);
    ol = (layers(2)-layers(1)+1);
    newl = ol*ref;
    newcartdim = [grdecl.cartDims(1:2), grdecl.cartDims(3)-ol+newl];
    zcorn = nan(newcartdim*2);
    zcorn(:, :, 1:2*(layers(1)-1)) = zcorn1;
    zcorn(:, :, end-size(zcorn2, 3)+1:end) = zcorn2;
    zcorn(:, :, 2*(layers(1)-1)+1:(end-size(zcorn2, 3)+2)) = czcorn;
    ngrdecl = struct('cartDims', newcartdim, 'COORD', grdecl.COORD, 'ZCORN', zcorn);

end

function Gt = flattenTraps(Gt, res)
% Given output res from trapAnalysis, produce a grid with flat top surfaces in trap areas.
%
% SYNOPSIS:
%   Gt = flattenTraps(Gt, res)
%
% PARAMETERS:
%   Gt  - Top surface grid
%
%   res - Output from trapAnalysis called on Gt.
%
% RETURNS:
%   Gt  - Gt where cells which are part of traps have had their z value
%   reduced to the level to the trap. Useful to make nice plots. This
%   function also adds the sortedCellNodes field to Gt if not present.
%
% SEE ALSO:
%   interactiveTrapping.

%{
#COPYRIGHT#
%}
    for i = 1:max(res.traps)
        Gt.cells.z(res.traps == i) = res.trap_z(i);
    end
    z_spill_loc = zeros(Gt.cells.num, 1);
    
    ind = res.traps == 0;
    z_spill_loc(ind) = 0;
    z_spill_loc(~ind) = res.trap_z(res.traps((~ind)));

    cells     = find(z_spill_loc > 0); 
    eIX       = Gt.cells.facePos;      
    nn        = double(diff([Gt.cells.facePos(cells), ...
                             Gt.cells.facePos(cells + 1)], [], 2));
                         
    Gt.cells.sortedCellNodes = getSortedCellNodes(Gt);


    cn        = double(Gt.cells.sortedCellNodes(mcolon(eIX(cells), eIX(cells + 1) - 1), 1));
    zz        = rldecode(z_spill_loc(z_spill_loc>0),nn);
    Gt.nodes.z(cn) = zz;

    Gt = computeGeometry(Gt);
end

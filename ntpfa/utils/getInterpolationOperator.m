function op = getInterpolationOperator(Pcell, Pface, Pnode, globalIndices, types)

%     Pcell = speye(G.cells.num);
    isCell = types == 1;
    isFace = types == 2;
    isNode = ~isCell & ~isFace;

    
    pcell = Pcell(globalIndices(isCell), :);
    pface = Pface(globalIndices(isFace), :);
    pnode = Pnode(globalIndices(isNode), :);
    
    op = @(pc) operator(pcell, pface, pnode, isCell, isFace, isNode, pc);
end

function val = operator(P_cell, P_face, P_node, isCell, isFace, isNode, v_cell)
    val = zeros(size(isCell));
    if isa(v_cell, 'ADI')
        val = double2ADI(val, v_cell);
    end
    if any(isCell)
        val(isCell) = P_cell*v_cell;
    end
    
    if any(isFace)
        val(isFace) = P_face*v_cell;
    end
    
    if any(isNode)
        val(isNode) = P_node*v_cell;
    end
end
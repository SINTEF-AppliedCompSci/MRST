function cnodes = activeCellNodes(Gt)
% Return a (4 x m)-sized matrix where m is the number of active cells in Gt.  Each
% column holds the indices of the 4 nodes that are corners of that cell.  
% The algorithm presupposes that each cell has exactly 4 corner nodes.
    
    cnode_full_mat = cellNodes(Gt);
    
    % only keep the information we want
    cnodes = reshape(cnode_full_mat(:, 3), 4, []);
       
end

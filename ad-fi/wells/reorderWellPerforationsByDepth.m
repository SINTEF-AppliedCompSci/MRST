function W = reorderWellPerforationsByDepth(W)
   for i = 1:numel(W)
      cells   = W(i).cells;
      WI      = W(i).WI;
      dZ      = W(i).dZ;
      cstatus = W(i).cstatus;
      
      index = (1:numel(cells))';
      new_index = sortrows(horzcat(dZ, index, cells, WI, cstatus));
      W(i).dZ      = new_index(:, 1);
      W(i).cells   = new_index(:, 3);
      W(i).WI      = new_index(:, 4);
      W(i).cstatus = new_index(:, 5);
   end
end
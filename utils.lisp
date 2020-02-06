(with-current-buffer "computeMPSATrans.m"
  (setq entries nil)
  (save-excursion
    (beginning-of-buffer)
    (while (re-search-forward (rx (or "crossTable" "SparseTensor()")) nil t)
      (push (buffer-substring-no-properties (line-beginning-position) (line-end-position)) entries)
      )
    ;; (beginning-of-buffer)
    ;; (while (re-search-forward "SparseTensor()" nil t)
      ;; (push (buffer-substring-no-properties (line-beginning-position) (line-end-position)) entries)
      ;; )
    )
  (setq entries (nreverse entries))
  )


(--map (insert it "\n") entries)

cellcoltbl = crossTable(celltbl, coltbl, {}); % ordering is cell - col
nodefacecoltbl = crossTable(nodefacetbl, coltbl, {});
cellnodefacetbl = crossTable(cellfacetbl, nodefacetbl, {'faces'});
cellnodecoltbl = crossTable(cellnodetbl, coltbl, {});
cellnodecolrowtbl = crossTable(cellnodecoltbl, rowtbl, {});
cellnodefacecoltbl = crossTable(cellnodefacetbl, coltbl, {});
cellnodefacecoltbl = crossTable(cellnodefacecoltbl, cellnodecoltbl, fds);
gradnodeface_T = SparseTensor();
gradcell_T = SparseTensor();
divnodeface_T = SparseTensor();
colrowtbl = crossTable(coltbl, rowtbl, {});
divcell_T = SparseTensor();
colrowtbl = crossTable(coltbl, rowtbl, {});
nodecolrowtbl = crossTable(nodetbl, colrowtbl, {});
trans_T = SparseTensor();
[~, indstruct] = crossTable(cellnodetbl, nodetbl, {'nodes'});
nodeaverage_T = SparseTensor();
celldispatch_T = SparseTensor();
cornercellnodecolrowtbl = crossTable(cornernodetbl, cellnodecolrowtbl, ...
cornerfix_T = SparseTensor();
colrowtbl = crossTable(coltbl, rowtbl, {});
col2row2tbl = crossTable(colrowtbl, colrowtbl, {}, 'crossextend', fds);
C_T = SparseTensor();
[cellnodecol2row2tbl, indstruct] = crossTable(cellnodetbl, col2row2tbl, {});
C_T = SparseTensor();
extnodefacetbl = crossTable(nodefacetbl, extfacetbl, {'faces'});
D_T = SparseTensor();
extnodefacetbl = crossTable(nodefacetbl, extfacetbl, {'faces'});
D_T = SparseTensor();
extnodefacetbl = crossTable(nodefacetbl, extfacetbl, {'faces'});
extnodefacecoltbl = crossTable(extnodefacetbl, coltbl, {});
    sourcetbl = crossTable(sourcetbl, coltbl, {});



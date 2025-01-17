function [tbls, mappings] = setupMpxaStandardTables(G, varargin)
% Setup the index arrays for both mpsa and mpfa

    opt = struct('useVirtual', true, ...
                 'inittbls', []);
    opt = merge_options(opt, varargin{:});
    useVirtual = opt.useVirtual;
    
    [tbls, mappings] = setupMpsaStandardTables(G, 'useVirtual', useVirtual, 'inittbls', opt.inittbls);

    inittbls.celltbl         = tbls.celltbl;
    inittbls.facetbl         = tbls.facetbl;
    inittbls.nodetbl         = tbls.nodetbl;
    inittbls.cellnodetbl     = tbls.cellnodetbl;
    inittbls.nodefacetbl     = tbls.nodefacetbl;
    inittbls.cellfacetbl     = tbls.cellfacetbl;
    inittbls.cellnodefacetbl = tbls.cellnodefacetbl;
    
    [tbls2, mappings2] = setupMpfaStandardTables(G, 'useVirtual', useVirtual, 'inittbls', inittbls);

    % merge the structures
    tbls     = mergeStructs(tbls, tbls2);
    mappings = mergeStructs(mappings, mappings2);
    
end

function st = mergeStructs(st, st2)

    fds = fieldnames(st2);

    for ifd = 1 : numel(fds)
        fd = fds{ifd};
        st.(fd) = st2.(fd);
    end
    
end


function [partition, cells] = makeGridPartition(model, cells, varargin)
    opt = struct('type'         , 'metis', ...
                 'computeTrans' , true   , ...
                 'minBlockSize' , []     , ...
                 'pdims'        , []     );
    opt = merge_options(opt, varargin{:});
    G  = model.G;
    if ~isempty(cells)
        Gr = removeCells(G, ~cells);
        Gr = computeGeometry(Gr);
    else
        Gr = G;
    end
    if isempty(opt.minBlockSize)
        opt.minBlockSize = G.cells.num/20;
    end
    nb = Gr.cells.num/opt.minBlockSize;
    switch opt.type
        case 'metis'
            rock = model.rock;
            rock.perm = rock.perm(cells,:);
%             T = model.operators.T;
            if opt.computeTrans
                T = getFaceTransmissibility(Gr, rock);
            else
                N = G.faces.neighbors + 1;
                c = [false; cells];
                keep = any(c(N), 2);
                T = model.operators.T_all(keep);
            end
            pg = partitionMETIS(Gr, T, nb);
        case 'cart'
            pdims = opt.pdims;
            if isempty(pdims)
                n = ceil(sqrt(nb));
                pdims = ones(1, G.griddim);
                pdims(1:2) = n;
            end
            if isfield(G, 'cartDims')
                pg = partitionUI(Gr, pdims);
            else
                pg = partitionCartGrid(pdims*50, pdims);
                pg = sampleFromBox(Gr, reshape(pg, pdims*50));
            end
    end
    pg = compressPartition(pg);
    pg = mergePartition(Gr, pg, opt.minBlockSize);
    partition = zeros(G.cells.num, 1);
    partition(cells) = pg;

end

%-----------------------------------------------------------------%
function partition = mergePartition(G, partition, minBlockSize)
    T = G.faces.areas;
    partition = mergeBlocksByConnections(G, partition, T, minBlockSize);
end
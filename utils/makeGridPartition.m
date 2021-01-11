function [pg, cells] = makeGridPartition(model, cells, varargin)
    opt = struct('type'        , 'metis', ...
                 'minBlockSize', []     , ...
                 'pdims'       , []     );
    opt = merge_options(opt, varargin{:});
    G  = model.G;
    Gr = removeCells(G, ~cells);
    Gr = computeGeometry(Gr);
    if isempty(opt.minBlockSize)
        opt.minBlockSize = Gr.cells.num/20;
    end
    nb = Gr.cells.num/opt.minBlockSize;
    switch opt.type
        case 'metis'
            rock = model.rock;
            rock.perm = rock.perm(cells,:);
%             T = model.operators.T;
            T = getFaceTransmissibility(Gr, rock);
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

end
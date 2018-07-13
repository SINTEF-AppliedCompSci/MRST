function [Gt, p] = perturb_topSurface(Gt, varargin)
% Add a perturbation to elevations of the top-surface

    opt.pert_type       = 'gaussian';
    opt.pert_interval   = [-15 15]; % meters
    opt.pert_clevel     = 1;        % coarsening level to apply to 
                                    % perturbations. Default is no
                                    % coarsening
    opt.kernal_size     = 3; % default, can pass in [3,3,3], [1,3,1], etc.
    opt.std             = 0.65;

    opt = merge_options(opt, varargin{:});

    if numel(opt.kernal_size)==1
       opt.kernal_size = [opt.kernal_size, opt.kernal_size, opt.kernal_size]; 
    end

    
    % -----------------------------
    % Generate a gaussian field
    % -----------------------------
    % for the surface change with dimensions [m,n] equal to Gt.cartDims+1
    % (node dimensions, not cell)

    surf_change_interval = opt.pert_interval;
    node_dims = Gt.cartDims+1;

    % Generate (a possibly coarsened) p field (then use interp to get p of
    % size Gt.cartDims+1)
    N = opt.pert_clevel;
    coarse_node_dims = ceil(Gt.cartDims./N)+1;
    if strcmpi(opt.pert_type,'gaussian')
        pc = gaussianField([coarse_node_dims], surf_change_interval, opt.kernal_size(1:2), opt.std);
    else
       error('Undefined perturbation type.') 
    end

    [xc,yc] = meshgrid(1:coarse_node_dims(2), 1:coarse_node_dims(1));
    r = 1;
    % r becomes smaller if we require a smaller step size to ensure enough
    % resolution for grid
    while 1
        dx = coarse_node_dims(1)/node_dims(1) * r;
        dy = coarse_node_dims(2)/node_dims(2) * r;
        [x,y] = meshgrid(1:dy:coarse_node_dims(2), 1:dx:coarse_node_dims(1));
        if size(x,1) >= node_dims(1) && size(y,2) >= node_dims(2)
            break
        else
            r = r*0.95;
        end
    end
    p = interp2(xc,yc,pc,x,y,'cubic');
    assert(all(size(p) >= node_dims))

    % -----------------------------
    % Add the gaussian field values to the node depths 
    % -----------------------------
    % This requires finding which index of field(i,j) corresponds to each
    % cell's node indexes

    mod_node_z = Gt.nodes.z; % to prevent over-writing
    for c = 1:Gt.cells.num
        % cell index is: c
        % ij index of this cell index:
        [i,j] = deal(Gt.cells.ij(c,1), Gt.cells.ij(c,2));
        % value of p at 4 surrounding nodes: order is counter-clock-wise, starting
        % at lower left corner: p(i,j); p(i+1,j); p(i+1,j+1); p(i,j+1)
        field_vals = [p(i,j); p(i+1,j); p(i+1,j+1); p(i,j+1)];
        % add these values of p to correct nodes in Gt.nodes.z:
        % cell nodes index:
        curr_cell_nodes = Gt.cells.sortedCellNodes((c*4-3):(c*4));
        mod_node_z( curr_cell_nodes ) = Gt.nodes.z( curr_cell_nodes ) + field_vals;
    end

    % -----------------------------
    % Apply change, recompute geometry
    % -----------------------------
    Gt.nodes.z = mod_node_z;
    Gt = computeGeometryVE_2D(Gt);


end
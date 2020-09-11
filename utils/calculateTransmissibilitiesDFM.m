function Trans = calculateTransmissibilitiesDFM(model) 

    fdom = model.G.FracturedDomains.domains{end};
    %% compute half transmissibilities for fracture virtual cells 
    half_trans = struct();
    half_trans.fracfrac = zeros(fdom.ncells,1);
    half_trans.fracmat = zeros(fdom.ncells,1);
    half_trans.matfrac_left = zeros(fdom.ncells,1);
    half_trans.matfrac_right = zeros(fdom.ncells,1);

    areas = model.G.faces.areas(fdom.region);

    %% Convenience variables
    k = fdom.rock.perm;
    a = fdom.rock.aperture;

    %% Setting frac-frac half trans
    half_trans.fracfrac = (a.*k)./(areas/2);

    %% Setting frac-mat half trans
    half_trans.fracmat = (2*areas.*k)./a;


    %% Setting mat-frac half trans
    matrix_cells_left = model.G.faces.neighbors(fdom.region,1);
    matrix_cells_right = model.G.faces.neighbors(fdom.region,2);
    half_trans.matfrac_left =  zeros(size(fdom.region));
    half_trans.matfrac_right =  zeros(size(fdom.region));
    

    for kk=1:size(fdom.region)
        facePos = model.G.cells.facePos(matrix_cells_left(kk)):model.G.cells.facePos(matrix_cells_left(kk)+1)-1;
        cfaces = model.G.cells.faces(facePos);
        Tid = facePos(find(cfaces==fdom.region(kk)));
        half_trans.matfrac_left(kk) = model.G.FracturedDomains.half_trans(Tid);
        facePos = model.G.cells.facePos(matrix_cells_right(kk)):model.G.cells.facePos(matrix_cells_right(kk)+1)-1;
        cfaces = model.G.cells.faces(facePos);
        Tid = facePos(find(cfaces==fdom.region(kk)));
        half_trans.matfrac_right(kk) = model.G.FracturedDomains.half_trans(Tid);
    end

    %% Assemble the frac-mat full transmissiblties            
    trans_left = half_trans.matfrac_left .* half_trans.fracmat./...
        (half_trans.matfrac_left + half_trans.fracmat);
    trans_right = half_trans.matfrac_right .* half_trans.fracmat./...
        (half_trans.matfrac_right + half_trans.fracmat);

    %% star-delta transformation for fracture intersections 
    trans_ff = zeros(size(fdom.connections_ff.virtual_cells,1),1);
    T_denom = zeros(size(fdom.fracture_node_info.nodes,1),1);
    for i=1:size(fdom.fracture_node_info.nodes,1)
        T_denom(i) =  sum(half_trans.fracfrac(fdom.fracture_node_info.conns_cells{i}));
    end
    for i=1:size(fdom.fracture_node_info.nodes,1)
        connsPos = fdom.fracture_node_info.conns_pos{i};
        if ~isempty(connsPos)
            for j=1:size(connsPos,2)
                local_fracture_cells = fdom.connections_ff.local_ids(connsPos(j),:);
                %star-delta transformation 
                T = half_trans.fracfrac(local_fracture_cells(1))*...
                    half_trans.fracfrac(local_fracture_cells(2))/T_denom(i);
                trans_ff(connsPos(j)) = T;
            end                    
        end
    end

    Trans = [trans_left;trans_right;trans_ff];
end

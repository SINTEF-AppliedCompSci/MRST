function T = getVerticallyAveragedTrans(Gt, rock)
% code taken from CO2VEBlackOilTypeModel()

    rock_tmp      = rock; 
    rock_tmp.perm = rock.perm .* Gt.cells.H; 
    T             = computeTrans(Gt, rock_tmp); 
    cf            = Gt.cells.faces(:, 1); 
    nf            = Gt.faces.num; 
    T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]);
    
end
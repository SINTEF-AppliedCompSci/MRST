function fill = findPercolationVolumes(traps, injectionAmount, trapvols)
    nt = numel(traps);
    fill = zeros(nt, 1);
    for i = nt:-1:1
        locvol = trapvols(traps(i));
        if injectionAmount > locvol
            fill(i) = 1;
            injectionAmount = injectionAmount - locvol;
        else
            fill(i) = injectionAmount/locvol;
            break
        end
    end
end
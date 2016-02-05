function poroU = upPoro(block)
% Upscale porosity of the block by pore volume averaging

poroU = sum(block.pv) / sum(block.G.cells.volumes);

end


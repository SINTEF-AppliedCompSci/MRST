function updata = upPoro(block, updata)
% Upscale porosity of the block by pore volume averaging

if nargin==1
    updata = [];
end

updata.poro = sum(block.pv) / sum(block.G.cells.volumes);

end


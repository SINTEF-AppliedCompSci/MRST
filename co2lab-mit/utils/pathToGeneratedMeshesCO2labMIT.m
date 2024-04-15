function pth = pathToGeneratedMeshesCO2labMIT()
    pth = fullfile(mrstOutputDirectory(), 'generated_meshes');
    if ~isfolder(pth)
        mkdir(pth)
    end
    return
end

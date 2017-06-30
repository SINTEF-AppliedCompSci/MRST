# ChemicalModel #

A model for solving equilibrium geochemistry.
            
## SYNOPSIS ##
    chem = ChemicalModel(varargin)
            
## DESCRIPTION ##
    A class of PhysicalModel which can construct and solve
    arbitrary aqueous geochemical system including surface chemistry
    all under the assumption of local chemical equilibrium. 

## REQUIRED PARAMETERS ##

### elementNames  ###      
    A cell array of strings containing all
    elements to be considered in the chemical system. 
    Elements do not have to correspond to actual element names,
    but it is reccomened that they do.

    ~~~~
    elementNames = {'H', 'O'};
    ~~~~

### speciesNames ###

    A cell array of strings of the chemical
    species to be considered in the chemical system. Each 
    entry in speciesNames must be a combination of
    entries found in elementNames, or surfaces. The species
    charge can be desginated with a +/- at the end of the
    name followed by the charge number (i.e. Ca+2, Cl-).

    ~~~~
    speciesNames = {'H+', 'OH-', 'H2O'};
    ~~~~

### reactions ###          
    A cell array listing the chemical
    reactions of the chemical system as strings, followed
    by the equilibrium constant of the reaction as a scalar
    quantity. Entries in reactions must be given as a
    string in the form 'reactants <-> products'. 
    Equuilibrium constant must be given in SI units.

    ~~~~
    reactions = {'H2O <-> H+ + OH-', 1e-14*mol/litre};
    ~~~~

## OPTIONAL PARAMETERS ##

### surfaces ###            
    A cell structure containing the names
    of surface functional groups, followed by information
    regarding the specific surface. The first entry of
    surfaces is the name of the surface functional group.
    Surface functional group names must begin with the ">"
    symbol (i.e. '>FeO'). These names shoudl be used in the same
    way an entries in elementNames for the construction of
    speciesNames. There are three different surface
    chemistry models implemented in this ChemicalModel: the
    langmuir model, the triple layer model, and the
    constant capacitance model. The diffuse layer and basic
    stern models can be constructed by an adjustment of the
    parameters of the triple layer model. A chemical system
    can have any number of surfaces with different surface
    complexation models.

    The Langmuir model is the most simplistic case.

    ~~~~
    geometry = [1/(nano*meter)^2 50*meter^2/gram 1000*grams/litre];

    surfaces = {'>FeO', {geometry, 'langmuir'}
    ~~~~

    In this example the surface functional group is
    '>FeO'. The second entry is a cell containing the
    surface parameters. The variable geometry must
    always be the first entry in this parameters cell.
    The entries of geometry must correspond to the
    surface site density, the specific surface area, and
    the slurry density in that specific order. Values
    must be given in SI units. The second entry of the
    parameters cell must be the type of surface model to
    use for the surface. In the case of 'langmuir' as above
    no other entries are needed. 

    The triple layer model requires additional work.

    ~~~~
    capacitances = [1 0.2]*Farad/meter^2;

    surfaces = {'>FeO', {geometry, 'tlm', capacitance,...
                                   '>FeO-', [-1 0 0]}
    ~~~~

    In this example geometry is defined as before, but now
    the '>FeO' is a triple layer model surface, as is
    indicated by 'tlm.' Additional parameters are need in
    this case. The capacitance density of each of the
    Helmholtz layers must be given after the model type.
    The basic stern model can be simulated by increasing
    the capacitance of the second layer (i.e. [1 1e3]) and
    the diffuse layer model can be simulated by increasing
    the capacitance of both layers (i.e. [1e3 1e3]). The
    capacitance is followed by the name of each surface
    species associated with '>FeO' and the species charge
    contribution to each plane. 

    The constant capacitance model can be similarly
    defined.
    
    ~~~~
    capacitances = 1*Farad/meter^2;

    surfaces = {'>FeO', {geometry, 'ccm', capacitance, '>FeO-', -1}
    ~~~~

    In this example geometry is defined as before, but now
    the '>FeO' is a constant capacitance surface, as is
    indicated by 'ccm.' Additional parameters are need in
    this case. The capacitance density of the
    Helmholtz layers must be given after the model type.
    The capacitance is followed by the name of each surface
    species associated with '>FeO' and the species charge
    contribution to the surface.       
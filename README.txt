====================================================================================================================
MRST's Poroelastic Dual-Continuum Module
====================================================================================================================
This is an MRST module for the modelling of multiscale poroelastic materials
within a dual-continuum setting. We have incoporated material anisotropy upto
orthotropy, which is reflected in the poroelastic constitutive equations, and
allowed for both materials to have intrinsic mechanical properties. 

The module itself was developed by Mark Ashworth of the Mutliscale / Multiphysics 
Group in the Institutite of GeoEnergy Engineering, Heriot-Watt University. In 
the module we leverage the Dual-Porosity module developed by Rafael March and 
Victoria Spooner, as well as the AD-mechanics module developed at Sintef. 
Any inquiries, sugestions or feedback should be sent to ma174@hw.ac.uk and/or 
f.doster@hw.ac.uk.

It is strongly suggested to look at the examples to see how the model is instatiated 
and used for the various levels of physical detail required (i.e. isotropy, anisotropy, etc). 

Description of folders and files:

- Base class models:
--- DualContinuumReservoirModel: Base flow class for all reservoirs that show dual-continuum behaviour.
--- DualContinuumMechanicModel: Base mechanics class for for all reservoirs that show dual-continuum 
                                behaviour. Constitutive coefficient models are instatiated here. 
--- DualContinuumWaterModel: Single phase water model that will be required by the fully coupled 
                             models.

- Models for fully coupled method:
--- DualContinuumMechFluidModel: Base class model to set up fully coupled mechanical-fluid simulations. 
                            This class is derived for each particular fluid model that is used.
--- DualContinuumMechMechanicModel: Mechanical part of the fully coupled model. Adds some functionalities 
                               required for the full coupling. 
--- DualContinuumWaterFluidModel: Fluid part of the fully coupled model. Derived from the DualContinuumWaterModel
                             and adds some functionalities required for the fully coupled simulation.
--- DualContinuumMechWaterModel: Single phase poromechanical dual-continuum model.

- utils:
--- addDerivedQuantitiesDC: Computes extra mechanical fields (such as strain, stress) from the primary
                            variable of the mechanical state (state.xd).
--- AnisoEnu2C: Anisotropic (orthotropic) stiffness tensor.
--- computeInitDisp: Computes the initial displacement.
--- computeInvC_dif: Computes the inverse of (C_m - C_f), where C denotes is
                     stiffness tensor notation. 
--- doubledot: Used to compute double inner product (up to orthotropic) 
--- equationsDPPoroMechanics: Assembles the residual for dual-continuum momentum balance.
--- equationsDPWaterMech: Assembles the residuals for the dual-continuum mass balance.

- transfer_models:
--- SimpleTransferFunction: Single phase transfer function as originally described in Warren and Root (1963).

- shape_factors:
--- Shape factor according to Lim and Aziz (1995). 
					  
- constitutive_coefficients:
--- ConstitutiveCoefficients: Base class for use in calculating constitutive coefficients belonging to the 
                              poromechanical dual-continuum constitutive model. 
--- isotropicVoidCoefficientModels: Dervied class which instatiates constitutive coefficients using models 
                                    from Khalili and Valliappan (1996), for which the underlying assumption 
                                    is that the high-perm, low storage phase is all void space.  
--- isotropicStiffCoefficientModels: Derived class which instantiates constitutive coefficients using models 
                                     in which the high-perm continuum also has an intrinsic stiffness (for
                                     example fractures exhibit intrinsic stiffness properties).
--- anisotropicStiffCoefficientModels: Derived class which instantiates constitutive coefficients using
                                       models in which the continua can show upto orthotropic material
                                       material behaviours.
NOTE: details on the isotropic and anisotropic stiff coefficient models will be released in a future
      Ashworth and Doster paper, so if confused just have search on the net for said work. At the 
      very least we will have tried to put this in arxiv.

- state_functions:
--- DC_PoroelasticPropertyFunctions: Derived state function grouping class for dual-continuum poroelasticity. 
--- EffectiveStress, PoroelasticFracturePoro, etc: State functions for poroelastic quantities.
					  								
- examples:
--- example_void_fractures: How to setup a simple 2D test case using the isotropicVoidCoefficientModels.
--- example_stiff_fractures: How to setup a simple 2D test case using the isotropicStiffCoefficientModels.
--- example_anisotropic_fractures: How to setup a simple 2D test case using the anisotropicStiffCoefficientModels.
--- DC_NorneExample: 3D dual-continuum-mech example on a geological grid.

References:
[1] Warren, J.E. and Root, P.J., 1963. The behavior of naturally fractured reservoirs.
[2] Lim, K.T. and Aziz, K., 1995. Matrix-fracture transfer shape factors for dual-porosity simulators. 
    Journal of Petroleum Science and Engineering, 13(3-4), pp.169-178.
[3] Khalili, N. and Valliappan, S., 1996. Unified theory of flow and deformation in double 
    porous media. European Journal of Mechanics, A/Solids, 15(2), pp.321-336.
[4] Ashworth, M. and Doster, F., 2019, March. An Open Source Numerical Framework for Dual-Continuum 
    Geomechanical Simulation. In SPE Reservoir Simulation Conference.
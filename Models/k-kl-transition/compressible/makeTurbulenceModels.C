#include "turbulentFluidThermoModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

namespace Foam
{

    typedef ThermalDiffusivity<CompressibleTurbulenceModel<fluidThermo> >
        fluidThermoCompressibleTurbulenceModel;

    typedef RASModel<EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
        RASfluidThermoCompressibleTurbulenceModel;
        
    typedef LESModel<EddyDiffusivity<fluidThermoCompressibleTurbulenceModel> >
        LESfluidThermoCompressibleTurbulenceModel;

}

///////////////////     RANS models      ///////////////////

#include "WAWDF.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    WAWDF
);

#include "RGamma.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    RGamma
);

#include "gammaReThetatSST.H"
makeTemplatedTurbulenceModel
(
    fluidThermoCompressibleTurbulenceModel,
    RAS,
    gammaReThetatSST
);


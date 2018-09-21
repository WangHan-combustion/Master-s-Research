#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "RASModel.H"
#include "LESModel.H"
#include "wallDist.H"

namespace Foam
{
    typedef IncompressibleTurbulenceModel<transportModel> 
        transportModelIncompressibleTurbulenceModel;
        
    typedef RASModel<transportModelIncompressibleTurbulenceModel>
        RAStransportModelIncompressibleTurbulenceModel;
    
    typedef LESModel<transportModelIncompressibleTurbulenceModel>
        LEStransportModelIncompressibleTurbulenceModel;
}

///////////////////     RANS models      ///////////////////
#include "WrayAgarwalWR2018.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwalWR2018
);

#include "WrayAgarwal2018.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2018
);

#include "WrayAgarwalWR.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwalWR
);

#include "SpalartAllmarasWR.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    SpalartAllmarasWR
);

#include "WrayAgarwal2017a.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    WrayAgarwal2017a
);

#include "GReWAmm.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    GReWAmm
);





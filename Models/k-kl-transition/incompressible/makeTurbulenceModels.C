#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "RASModel.H"
#include "LESModel.H"

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
/*
#include "kkloneOmega.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkloneOmega
);
*/
#include "kkloneTranV2.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkloneTranV2
);

#include "kkloneTranMenter.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkloneTranMenter
);

#include "kkloneTranV3.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkloneTranV3
);

#include "GReWA.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    GReWA
);

#include "kklone.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kklone
);
/*
#include "kkl.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkl
);

#include "kkloneM.H"
makeTemplatedTurbulenceModel
(
    transportModelIncompressibleTurbulenceModel,
    RAS,
    kkloneM
);







*/



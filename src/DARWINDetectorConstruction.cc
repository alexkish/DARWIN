#include <G4Material.hh>
#include <G4NistManager.hh>
#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include "G4Torus.hh"
#include <G4Orb.hh>
#include <G4Polyhedra.hh>
#include <G4Polycone.hh>
#include <G4Ellipsoid.hh>
#include <G4Trd.hh>
#include <G4Cons.hh>
#include <G4UnionSolid.hh>
#include <G4IntersectionSolid.hh>
#include <G4SubtractionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4SDManager.hh>
#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <globals.hh>

#include <vector>
#include <numeric>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cassert>

using std::vector;
using std::stringstream;
using std::max;

#include "DARWINLXeSensitiveDetector.hh"

#include "DARWINDetectorConstruction.hh"
#include "DARWINDetectorMessenger.hh"

map<G4String, G4double> DARWINDetectorConstruction::m_hGeometryParameters;

DARWINDetectorConstruction::DARWINDetectorConstruction()
{
	m_pRotationXPlus225 = new G4RotationMatrix();
	m_pRotationXPlus225->rotateX(22.5*deg);

	m_pRotationXMinus225 = new G4RotationMatrix();
	m_pRotationXMinus225->rotateX(-22.5*deg);

	m_pRotationXPlus45 = new G4RotationMatrix();
	m_pRotationXPlus45->rotateX(45.*deg);

	m_pRotationXMinus45 = new G4RotationMatrix();
	m_pRotationXMinus45->rotateX(-45.*deg);

	m_pRotationXPlus90 = new G4RotationMatrix();
	m_pRotationXPlus90->rotateX(90.*deg);

	m_pRotationXMinus90 = new G4RotationMatrix();
	m_pRotationXMinus90->rotateX(-90.*deg);

	m_pRotationX180 = new G4RotationMatrix();
	m_pRotationX180->rotateX(180.*deg);

	RotM90y  = new G4RotationMatrix() ;
	RotP90y  = new G4RotationMatrix() ;
	RotM90x  = new G4RotationMatrix() ;
	RotP270x  = new G4RotationMatrix() ;
	RotP90x  = new G4RotationMatrix() ;
	RotP90z  = new G4RotationMatrix() ;
	RotM90z  = new G4RotationMatrix() ;
	RotP0x   = new G4RotationMatrix() ;
	RotP180x = new G4RotationMatrix() ;

	RotP0x->rotateX   (0. *deg);
	RotP90x->rotateX (90. *deg);
	RotP270x->rotateX (270. *deg);
	RotM90x->rotateX(-90. *deg);
	RotP180x->rotateX(180. *deg);
	RotP90y->rotateY (90. *deg);
	RotM90y->rotateY(-90. *deg);
	RotP90z->rotateZ (90. *deg);
	RotM90z->rotateZ(-90. *deg);

	m_pDetectorMessenger = new DARWINDetectorMessenger(this);
}

DARWINDetectorConstruction::~DARWINDetectorConstruction()
{
		delete m_pDetectorMessenger;
}

G4VPhysicalVolume*
DARWINDetectorConstruction::Construct()
{
	DefineMaterials();
	G4cout<<"Constructing Materials"<<G4endl;

	DefineGeometryParameters();
	G4cout<<"Constructing Geometry"<<G4endl;

	ConstructLaboratory(); G4cout<<"Constructing Lab"<<G4endl;

	ConstructVeto(); G4cout<<"Constructing Veto"<<G4endl; // water buffer
	//ConstructLSMVeto(); // in alternative to the water tank veto

	ConstructCryostat(); G4cout<<"Constructing Cryostat"<<G4endl;

	ConstructXenon(); G4cout<<"Constructing Xenon"<<G4endl;

	ConstructTPC(); G4cout<<"Constructing TPC"<<G4endl;

	//ConstructSensitiveLXe(); G4cout<<"Constructing Sensitive LXe"<<G4endl;

	ConstructFieldCage(); G4cout<<"Constructing Field Cage"<<G4endl;

	ConstructQUPIDArrays(); G4cout<<"Constructing QUPID arrays"<<G4endl;

	//ConstructVetoPMTArrays(); G4cout<<"Constructing Veto PMT arrays"<<G4endl;

	//CheckOverlapping();

	//PrintGeometryInformation();

	return m_pLabPhysicalVolume;
}

void
DARWINDetectorConstruction::DefineMaterials()
{
	G4NistManager* pNistManager = G4NistManager::Instance();

	//================================== elements ===================================
	G4Element *Xe = new G4Element("Xenon",     "Xe", 54., 131.293*g/mole);
	G4Element *H  = new G4Element("Hydrogen",  "H",  1.,  1.0079*g/mole);
	G4Element *C  = new G4Element("Carbon",    "C",  6.,  12.011*g/mole);
	G4Element *N  = new G4Element("Nitrogen",  "N",  7.,  14.007*g/mole);
	G4Element *O  = new G4Element("Oxygen",    "O",  8.,  15.999*g/mole);
	G4Element *F  = new G4Element("Fluorine",  "F",  9.,  18.998*g/mole);
	G4Element *Al = new G4Element("Aluminium", "Al", 13., 26.982*g/mole);
	G4Element *Si = new G4Element("Silicon",   "Si", 14., 28.086*g/mole);
	G4Element *Cr = new G4Element("Chromium",  "Cr", 24., 51.996*g/mole);
	G4Element *Mn = new G4Element("Manganese", "Mn", 25., 54.938*g/mole);
	G4Element *Fe = new G4Element("Iron",      "Fe", 26., 55.85*g/mole);
	G4Element *Ni = new G4Element("Nickel",    "Ni", 28., 58.693*g/mole);
	G4Element *Cu = new G4Element("Copper",    "Cu", 29., 63.546*g/mole);
	//G4Element* Ti = new G4Element("Titanium",  "Ti", 22., 47.867*g/mole);
	G4Element *Ti = pNistManager->FindOrBuildElement("Ti");
//	G4Element* Mg = new G4Element("Magnesium", "Mg", 12., 24.3050*g/mole);
//	G4Element *B  = pNistManager->FindOrBuildElement("B");
//	G4Element *Gd = pNistManager->FindOrBuildElement("Gd");
	
	//================================== materials ================================== 

	//------------------------------------- air -------------------------------------
	pNistManager->FindOrBuildMaterial("G4_AIR");

	//----------------------------------- vacuum ------------------------------------
	G4Material *Vacuum = new G4Material("Vacuum", 1.e-20*g/cm3, 2, kStateGas);
	Vacuum->AddElement(N, 0.755);
	Vacuum->AddElement(O, 0.245);

	//------------------------------------ water ------------------------
	G4Material *Water = new G4Material("Water", 1.*g/cm3, 2, kStateLiquid);
	Water->AddElement(H, 2);
	Water->AddElement(O, 1);

	//-------------------------------- liquid xenon ---------------------------------
	//G4Material *LXe = new G4Material("LXe", 2.85*g/cm3, 1, kStateLiquid, 168.15*kelvin, 1.5*atmosphere);
	G4Material *LXe = new G4Material("LXe", 2.828*g/cm3, 1, kStateLiquid, 168.15*kelvin, 1.5*atmosphere);
	LXe->AddElement(Xe, 1);

	const G4int iNbEntries = 3;

	G4double pdLXePhotonMomentum[iNbEntries]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdLXeScintillation[iNbEntries]    = {0.1,     1.0,     0.1};
	G4double pdLXeRefractiveIndex[iNbEntries]  = {1.63,    1.61,    1.58};
	G4double pdLXeAbsorbtionLength[iNbEntries] = {100.*cm, 100.*cm, 100.*cm};
	G4double pdLXeScatteringLength[iNbEntries] = {30.*cm,  30.*cm,  30.*cm};

	G4MaterialPropertiesTable *pLXePropertiesTable = new G4MaterialPropertiesTable();
	
	pLXePropertiesTable->AddProperty("FASTCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("SLOWCOMPONENT", pdLXePhotonMomentum, pdLXeScintillation, iNbEntries);
	pLXePropertiesTable->AddProperty("RINDEX", pdLXePhotonMomentum, pdLXeRefractiveIndex, iNbEntries);
	pLXePropertiesTable->AddProperty("ABSLENGTH", pdLXePhotonMomentum, pdLXeAbsorbtionLength, iNbEntries);
	pLXePropertiesTable->AddProperty("RAYLEIGH", pdLXePhotonMomentum, pdLXeScatteringLength, iNbEntries);
	
	pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./keV);
	pLXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
	pLXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pLXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pLXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);
	
	LXe->SetMaterialPropertiesTable(pLXePropertiesTable);

	//-------------------------------- gaseous xenon --------------------------------
	G4Material *GXe = new G4Material("GXe", 0.005887*g/cm3, 1, kStateGas, 173.15*kelvin, 1.5*atmosphere);
	GXe->AddElement(Xe, 1);

	G4double pdGXePhotonMomentum[iNbEntries]   = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdGXeScintillation[iNbEntries]    = {0.1,     1.0,     0.1};
	G4double pdGXeRefractiveIndex[iNbEntries]  = {1.00,    1.00,    1.00};
	G4double pdGXeAbsorbtionLength[iNbEntries] = {100*m,   100*m,   100*m};
	G4double pdGXeScatteringLength[iNbEntries] = {100*m,   100*m,   100*m};

	G4MaterialPropertiesTable *pGXePropertiesTable = new G4MaterialPropertiesTable();

	pGXePropertiesTable->AddProperty("FASTCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("SLOWCOMPONENT", pdGXePhotonMomentum, pdGXeScintillation, iNbEntries);
	pGXePropertiesTable->AddProperty("RINDEX", pdGXePhotonMomentum, pdGXeRefractiveIndex, iNbEntries);
	pGXePropertiesTable->AddProperty("ABSLENGTH", pdGXePhotonMomentum, pdGXeAbsorbtionLength, iNbEntries);
	pGXePropertiesTable->AddProperty("RAYLEIGH", pdGXePhotonMomentum, pdGXeScatteringLength, iNbEntries);

	pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 0./(keV));
	pGXePropertiesTable->AddConstProperty("RESOLUTIONSCALE", 0);
	pGXePropertiesTable->AddConstProperty("FASTTIMECONSTANT", 3.*ns);
	pGXePropertiesTable->AddConstProperty("SLOWTIMECONSTANT", 27.*ns);
	pGXePropertiesTable->AddConstProperty("YIELDRATIO", 1.0);

	GXe->SetMaterialPropertiesTable(pGXePropertiesTable);	

	//----------------------------------- quartz ------------------------------------
	// ref: http://www.sciner.com/Opticsland/FS.htm
	G4Material *Quartz = new G4Material("Quartz", 2.201*g/cm3, 2, kStateSolid, 168.15*kelvin, 1.5*atmosphere);
	Quartz->AddElement(Si, 1);
	Quartz->AddElement(O, 2);

	const 	G4int iNbEntriesMatch = 5;
	G4double pdQuartzPhotonMomentum[iNbEntriesMatch]   = {2.*eV, 6.9*eV, 6.91*eV, 6.98*eV, 7.05*eV}; 
	G4double pdQuartzRefractiveIndex[iNbEntriesMatch]  = { 1.50,   1.50,    1.50,    1.56,    1.60};
	G4double pdQuartzAbsorbtionLength[iNbEntriesMatch] = {30*m,    30*m,    30*m,    30*m,    30*m};

	G4MaterialPropertiesTable *pQuartzPropertiesTable = new G4MaterialPropertiesTable();

	pQuartzPropertiesTable->AddProperty("RINDEX", pdQuartzPhotonMomentum, pdQuartzRefractiveIndex, iNbEntriesMatch);
	pQuartzPropertiesTable->AddProperty("ABSLENGTH", pdQuartzPhotonMomentum, pdQuartzAbsorbtionLength, iNbEntriesMatch);

	Quartz->SetMaterialPropertiesTable(pQuartzPropertiesTable);

	//------------------------------- stainless steel -------------------------------
	G4Material *SS304LSteel = new G4Material("SS304LSteel", 8.00*g/cm3, 5, kStateSolid);
	SS304LSteel->AddElement(Fe, 0.65);
	SS304LSteel->AddElement(Cr, 0.20);
	SS304LSteel->AddElement(Ni, 0.12);
	SS304LSteel->AddElement(Mn, 0.02);
	SS304LSteel->AddElement(Si, 0.01);

	//---------------------------------- titanium ----------------------------------
	G4Material *Titanium = new G4Material("Titanium", 4.506*g/cm3, 1, kStateSolid);
	Titanium->AddElement(Ti, 1); 

	//---------------------------- photocathode aluminium ---------------------------
	G4Material *PhotoCathodeAluminium = new G4Material("PhotoCathodeAluminium", 8.00*g/cm3, 1, kStateSolid);
	PhotoCathodeAluminium->AddElement(Al, 1);

	G4double pdPhotoCathodePhotonMomentum[]   = {2.*eV, 6.9*eV, 6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdPhotoCathodeRefractiveIndex[]  = {1.50,  1.50,   1.50,    1.56,    1.60};
	G4double pdPhotoCathodeAbsorbtionLength[] = {1.*nm, 1*nm,   1.*nm,   1.*nm,   1.*nm};
	G4MaterialPropertiesTable *pPhotoCathodePropertiesTable = new G4MaterialPropertiesTable();

	pPhotoCathodePropertiesTable->AddProperty("RINDEX", pdPhotoCathodePhotonMomentum, pdPhotoCathodeRefractiveIndex, iNbEntriesMatch);
	pPhotoCathodePropertiesTable->AddProperty("ABSLENGTH", pdPhotoCathodePhotonMomentum, pdPhotoCathodeAbsorbtionLength, iNbEntriesMatch);

	PhotoCathodeAluminium->SetMaterialPropertiesTable(pPhotoCathodePropertiesTable);

	//---------------------------- QUPID aluminium coating---------------------------
	G4Material *CoatingAluminium = new G4Material("CoatingAluminium", 2.7*g/cm3, 1, kStateSolid);
	CoatingAluminium->AddElement(Al, 1);

	G4double pdCoatingAluminiumPhotonMomentum[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdCoatingAluminiumReflectivity[iNbEntries]   = {0.15,    0.2,     0.15};/// check
	G4MaterialPropertiesTable *pCoatingAluminiumPropertiesTable = new G4MaterialPropertiesTable();

	pCoatingAluminiumPropertiesTable->AddProperty("REFLECTIVITY", pdCoatingAluminiumPhotonMomentum, pdCoatingAluminiumReflectivity, iNbEntries);
	CoatingAluminium->SetMaterialPropertiesTable(pCoatingAluminiumPropertiesTable);

	//----------------------------- grid mesh aluminium------------------------------
	G4Material *GridMeshAluminium = new G4Material("GridMeshAluminium", 0.174*g/cm3, 1, kStateSolid);
	GridMeshAluminium->AddElement(Al, 1);
	
	G4double pdGridMeshPhotonMomentum[] = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double *pdGridMeshRefractiveIndex = pdLXeRefractiveIndex;
	G4double dAbsortionLength = 2.10*mm; // exp(-GridMeshThickness/2.10) = 0.94
	G4double pdGridMeshAbsortionLength[] = {dAbsortionLength, dAbsortionLength, dAbsortionLength};
	
	G4MaterialPropertiesTable *pGridMeshPropertiesTable = new G4MaterialPropertiesTable();

	pGridMeshPropertiesTable->AddProperty("RINDEX", pdGridMeshPhotonMomentum, pdGridMeshRefractiveIndex, iNbEntries);
	pGridMeshPropertiesTable->AddProperty("ABSLENGTH", pdGridMeshPhotonMomentum, pdGridMeshAbsortionLength, iNbEntries);
	GridMeshAluminium->SetMaterialPropertiesTable(pGridMeshPropertiesTable);

	//------------------------------------ teflon -----------------------------------
	G4Material* Teflon = new G4Material("Teflon", 2.2*g/cm3, 2, kStateSolid);
	Teflon->AddElement(C, 0.240183);
	Teflon->AddElement(F, 0.759817);

	G4double pdTeflonPhotonMomentum[iNbEntries]  = {6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdTeflonRefractiveIndex[iNbEntries] = {1.63,    1.61,    1.58};
	G4double pdTeflonReflectivity[iNbEntries]    = {0.95,    0.95,    0.95};
	G4double pdTeflonSpecularLobe[iNbEntries]    = {0.01,    0.01,    0.01};
	G4double pdTeflonSpecularSpike[iNbEntries]   = {0.01,    0.01,    0.01};
	G4double pdTeflonBackscatter[iNbEntries]     = {0.01,    0.01,    0.01};
	G4double pdTeflonEfficiency[iNbEntries]      = {1.0,     1.0,     1.0};
	G4MaterialPropertiesTable *pTeflonPropertiesTable = new G4MaterialPropertiesTable();

	pTeflonPropertiesTable->AddProperty("RINDEX", pdTeflonPhotonMomentum, pdTeflonRefractiveIndex, iNbEntries);
	pTeflonPropertiesTable->AddProperty("REFLECTIVITY", pdTeflonPhotonMomentum, pdTeflonReflectivity, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARLOBECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularLobe, iNbEntries);
	pTeflonPropertiesTable->AddProperty("SPECULARSPIKECONSTANT", pdTeflonPhotonMomentum, pdTeflonSpecularSpike, iNbEntries);
	pTeflonPropertiesTable->AddProperty("BACKSCATTERCONSTANT", pdTeflonPhotonMomentum, pdTeflonBackscatter, iNbEntries);
	pTeflonPropertiesTable->AddProperty("EFFICIENCY", pdTeflonPhotonMomentum, pdTeflonEfficiency, iNbEntries);
	Teflon->SetMaterialPropertiesTable(pTeflonPropertiesTable);

	//------------------------------------ copper -----------------------------------
	G4Material *Copper = new G4Material("Copper", 8.92*g/cm3, 1);
	Copper->AddElement(Cu, 1); 

	G4double pdCopperPhotonMomentum[iNbEntries] = {1.91*eV, 6.98*eV, 7.05*eV};//{6.91*eV, 6.98*eV, 7.05*eV};
	G4double pdCopperReflectivity[iNbEntries]   = {0.15,    0.2,     0.15};
	G4MaterialPropertiesTable *pCopperPropertiesTable = new G4MaterialPropertiesTable();

	pCopperPropertiesTable->AddProperty("REFLECTIVITY", pdCopperPhotonMomentum, pdCopperReflectivity, iNbEntries);
	Copper->SetMaterialPropertiesTable(pCopperPropertiesTable);

	//------------------------------------ polyethilene -----------------------------------
	G4Material *  poly = new G4Material("poly", 0.95*g/cm3, 2);
	poly->AddElement(C, 1);
	poly->AddElement(H, 2);  

	// ----------- optical properties of water --------------------
  const G4int nEntries = 32;

  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
//
// Water
//	      
  G4double RefractiveIndex1[nEntries] =
            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
              1.346,  1.3465, 1.347,  1.3475, 1.348,
              1.3485, 1.3492, 1.35,   1.3505, 1.351,
              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
              1.3545, 1.355,  1.3555, 1.356,  1.3568,
              1.3572, 1.358,  1.3585, 1.359,  1.3595,
              1.36,   1.3608};

  G4double Absorption1[nEntries] =
           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
           17.500*m, 14.500*m };

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",       PhotonEnergy, RefractiveIndex1,nEntries);
  myMPT1->AddProperty("ABSLENGTH",    PhotonEnergy, Absorption1,     nEntries);
  
  Water->SetMaterialPropertiesTable(myMPT1);
}

void
DARWINDetectorConstruction::DefineGeometryParameters()
{
	//================================== Laboratory =================================
	m_hGeometryParameters["LabHeight"]			= 12.*m;
	m_hGeometryParameters["LabRadius"]			= 6.*m;
	//================================== Water tank =================================
	m_hGeometryParameters["WaterTankThickness"]		= 2.*mm;
	m_hGeometryParameters["WaterTankDomeOuterHeight"]	= 60.*cm;
	m_hGeometryParameters["WaterTankDomeInnerHeight"]	= GetGeometryParameter("WaterTankDomeOuterHeight") - GetGeometryParameter("WaterTankThickness");
	m_hGeometryParameters["WaterTankCylinderHeight"]	= 10.*m - GetGeometryParameter("WaterTankDomeOuterHeight");  // MS it was 3.062 to have 1m of water on top
	m_hGeometryParameters["WaterTankCylinderInnerHeight"]	= GetGeometryParameter("WaterTankCylinderHeight") - GetGeometryParameter("WaterTankThickness");
	m_hGeometryParameters["WaterTankOuterRadius"]		= 5.*m;
	m_hGeometryParameters["WaterTankInnerRadius"]		= GetGeometryParameter("WaterTankOuterRadius") - GetGeometryParameter("WaterTankThickness");
	//================================== Detector dimensions =================================
	m_hGeometryParameters["FiducialDriftLength"]	= 200.385*cm; // z_e
	m_hGeometryParameters["FiducialRadius"]			= 82.1*cm; // r_e
	m_hGeometryParameters["LinearFiducialCut"]		= 13.*cm; // GetGeometryParameter("LinearFiducialCut")

// 5 TONS Fiducial
//	z_e = 154.521 cm	r_e: 60.35 cm	LinearFiducialCut: 10 cm	Tot mass: 7.67369 tons	R/r: 0.050462 	Filling factor: 0.847438	Number of QUPIDs required per array: 333
//	z_e = 123.702 cm	r_e: 67.45 cm	LinearFiducialCut: 10 cm	Tot mass: 7.65835 tons	R/r: 0.045836 	Filling factor: 0.850362	Number of QUPIDs required per array: 405
//	z_e = 101.262 cm	r_e: 74.55 cm	LinearFiducialCut: 10 cm	Tot mass: 7.70159 tons	R/r: 0.041987 	Filling factor: 0.858943	Number of QUPIDs required per array: 488
//	z_e = 84.4167 cm	r_e: 81.65 cm	LinearFiducialCut: 10 cm	Tot mass: 7.79227 tons	R/r: 0.0387343	Filling factor: 0.861881	Number of QUPIDs required per array: 575
//	z_e = 71.4503 cm	r_e: 88.75 cm	LinearFiducialCut: 10 cm	Tot mass: 7.92298 tons	R/r: 0.0359494	Filling factor: 0.860341	Number of QUPIDs required per array: 666
//	z_e = 154.521 cm	r_e: 60.35 cm	LinearFiducialCut: 12 cm	Tot mass: 8.30223 tons	R/r: 0.049067 	Filling factor: 0.848882	Number of QUPIDs required per array: 353
//	z_e = 123.702 cm	r_e: 67.45 cm	LinearFiducialCut: 12 cm	Tot mass: 8.2833  tons	R/r: 0.0446822	Filling factor: 0.853096	Number of QUPIDs required per array: 428
//	z_e = 101.262 cm	r_e: 74.55 cm	LinearFiducialCut: 12 cm	Tot mass: 8.33646 tons	R/r: 0.0410168	Filling factor: 0.854504	Number of QUPIDs required per array: 508
//	z_e = 84.4167 cm	r_e: 81.65 cm	LinearFiducialCut: 12 cm	Tot mass: 8.44775 tons	R/r: 0.0379071	Filling factor: 0.859077	Number of QUPIDs required per array: 598
//	z_e = 71.4503 cm	r_e: 88.75 cm	LinearFiducialCut: 12 cm	Tot mass: 8.60789 tons	R/r: 0.0352357	Filling factor: 0.862761	Number of QUPIDs required per array: 695
//	z_e = 154.521 cm	r_e: 60.35 cm	LinearFiducialCut: 14 cm	Tot mass: 8.96402 tons	R/r: 0.0477471	Filling factor: 0.847788	Number of QUPIDs required per array: 372
//	z_e = 123.702 cm	r_e: 67.45 cm	LinearFiducialCut: 14 cm	Tot mass: 8.94134 tons	R/r: 0.043585 	Filling factor: 0.854317	Number of QUPIDs required per array: 450
//	z_e = 101.262 cm	r_e: 74.55 cm	LinearFiducialCut: 14 cm	Tot mass: 9.00485 tons	R/r: 0.0400903	Filling factor: 0.858174	Number of QUPIDs required per array: 534
//	z_e = 84.4167 cm	r_e: 81.65 cm	LinearFiducialCut: 14 cm	Tot mass: 9.13756 tons	R/r: 0.0371145	Filling factor: 0.859123	Number of QUPIDs required per array: 624
//	z_e = 71.4503 cm	r_e: 88.75 cm	LinearFiducialCut: 14 cm	Tot mass: 9.32822 tons	R/r: 0.0345499	Filling factor: 0.864702	Number of QUPIDs required per array: 725
//	z_e = 154.521 cm	r_e: 60.35 cm	LinearFiducialCut: 16 cm	Tot mass: 9.65993 tons	R/r: 0.0464964	Filling factor: 0.850034	Number of QUPIDs required per array: 394
//	z_e = 123.702 cm	r_e: 67.45 cm	LinearFiducialCut: 16 cm	Tot mass: 9.63333 tons	R/r: 0.0425404	Filling factor: 0.85378 	Number of QUPIDs required per array: 472
//	z_e = 101.262 cm	r_e: 74.55 cm	LinearFiducialCut: 16 cm	Tot mass: 9.7076  tons	R/r: 0.0392049	Filling factor: 0.860074	Number of QUPIDs required per array: 560
//	z_e = 84.4167 cm	r_e: 81.65 cm	LinearFiducialCut: 16 cm	Tot mass: 9.86255 tons	R/r: 0.0363543	Filling factor: 0.86185 	Number of QUPIDs required per array: 653
//	z_e = 71.4503 cm	r_e: 88.75 cm	LinearFiducialCut: 16 cm	Tot mass: 10.0848 tons	R/r: 0.0338902	Filling factor: 0.862686	Number of QUPIDs required per array: 752
//	z_e = 154.521 cm	r_e: 60.35 cm	LinearFiducialCut: 18 cm	Tot mass: 10.3908 tons	R/r: 0.0453095	Filling factor: 0.855313	Number of QUPIDs required per array: 417
//	z_e = 123.702 cm	r_e: 67.45 cm	LinearFiducialCut: 18 cm	Tot mass: 10.3601 tons	R/r: 0.0415448	Filling factor: 0.858518	Number of QUPIDs required per array: 500
//	z_e = 101.262 cm	r_e: 74.55 cm	LinearFiducialCut: 18 cm	Tot mass: 10.4456 tons	R/r: 0.0383576	Filling factor: 0.858488	Number of QUPIDs required per array: 584
//	z_e = 84.4167 cm	r_e: 81.65 cm	LinearFiducialCut: 18 cm	Tot mass: 10.6236 tons	R/r: 0.0356247	Filling factor: 0.861336	Number of QUPIDs required per array: 679
//	z_e = 71.4503 cm	r_e: 88.75 cm	LinearFiducialCut: 18 cm	Tot mass: 10.8786 tons	R/r: 0.0332553	Filling factor: 0.865797	Number of QUPIDs required per array: 784

// 10 TONS fiducial
//  z_e = 240.121 cm	r_e: 75		 cm	LinearFiducialCut: 10	cm	Tot mass: 16.6971	tons	R/r: 0.0417647	Filling factor: 0.860015	Number of QUPIDs required per array:  494
//  z_e = 200.385 cm	r_e: 82.1	 cm	LinearFiducialCut: 10	cm	Tot mass: 16.6085	tons	R/r: 0.0385451	Filling factor: 0.86048		Number of QUPIDs required per array:  580
//  z_e = 169.755 cm	r_e: 89.2	 cm	LinearFiducialCut: 10	cm	Tot mass: 16.59		tons	R/r: 0.0357863	Filling factor: 0.860896	Number of QUPIDs required per array:  673
//  z_e = 145.646 cm	r_e: 96.3	 cm	LinearFiducialCut: 10	cm	Tot mass: 16.6294	tons	R/r: 0.033396	Filling factor: 0.864722	Number of QUPIDs required per array:  776
//  z_e = 126.331 cm	r_e: 103.4	 cm	LinearFiducialCut: 10	cm	Tot mass: 16.7183	tons	R/r: 0.0313051	Filling factor: 0.866866	Number of QUPIDs required per array:  885
//  z_e = 240.121 cm	r_e: 75		 cm	LinearFiducialCut: 12	cm	Tot mass: 17.7611	tons	R/r: 0.0408046	Filling factor: 0.854988	Number of QUPIDs required per array:  514
//  z_e = 200.385 cm	r_e: 82.1	 cm	LinearFiducialCut: 12	cm	Tot mass: 17.6524	tons	R/r: 0.0377258	Filling factor: 0.859456	Number of QUPIDs required per array:  604
//  z_e = 169.755 cm	r_e: 89.2	 cm	LinearFiducialCut: 12	cm	Tot mass: 17.6296	tons	R/r: 0.0350791	Filling factor: 0.861848	Number of QUPIDs required per array:  701
//  z_e = 145.646 cm	r_e: 96.3	 cm	LinearFiducialCut: 12	cm	Tot mass: 17.6779	tons	R/r: 0.0327793	Filling factor: 0.863916	Number of QUPIDs required per array:  806
//  z_e = 126.331 cm	r_e: 103.4	 cm	LinearFiducialCut: 12	cm	Tot mass: 17.7865	tons	R/r: 0.0307626	Filling factor: 0.868016	Number of QUPIDs required per array:  918
//  z_e = 240.121 cm	r_e: 75		 cm	LinearFiducialCut: 14	cm	Tot mass: 18.8686	tons	R/r: 0.0398876	Filling factor: 0.857919	Number of QUPIDs required per array:  540
//  z_e = 200.385 cm	r_e: 82.1	 cm	LinearFiducialCut: 14	cm	Tot mass: 18.7389	tons	R/r: 0.0369407	Filling factor: 0.861737	Number of QUPIDs required per array:  632
//  z_e = 169.755 cm	r_e: 89.2	 cm	LinearFiducialCut: 14	cm	Tot mass: 18.7118	tons	R/r: 0.0343992	Filling factor: 0.865304	Number of QUPIDs required per array:  732
//  z_e = 145.646 cm	r_e: 96.3	 cm	LinearFiducialCut: 14	cm	Tot mass: 18.7692	tons	R/r: 0.032185	Filling factor: 0.866948	Number of QUPIDs required per array:  837
//  z_e = 126.331 cm	r_e: 103.4	 cm	LinearFiducialCut: 14	cm	Tot mass: 18.8982	tons	R/r: 0.0302385	Filling factor: 0.865312	Number of QUPIDs required per array:  947
//  z_e = 240.121 cm	r_e: 75		 cm	LinearFiducialCut: 16	cm	Tot mass: 20.0204	tons	R/r: 0.039011	Filling factor: 0.861783	Number of QUPIDs required per array:  567
//  z_e = 200.385 cm	r_e: 82.1	 cm	LinearFiducialCut: 16	cm	Tot mass: 19.869	tons	R/r: 0.0361876	Filling factor: 0.862057	Number of QUPIDs required per array:  659
//  z_e = 169.755 cm	r_e: 89.2	 cm	LinearFiducialCut: 16	cm	Tot mass: 19.8374	tons	R/r: 0.0337452	Filling factor: 0.863611	Number of QUPIDs required per array:  759
//  z_e = 145.646 cm	r_e: 96.3	 cm	LinearFiducialCut: 16	cm	Tot mass: 19.9042	tons	R/r: 0.0316118	Filling factor: 0.865864	Number of QUPIDs required per array:  867
//  z_e = 126.331 cm	r_e: 103.4	 cm	LinearFiducialCut: 16	cm	Tot mass: 20.0542	tons	R/r: 0.029732	Filling factor: 0.867723	Number of QUPIDs required per array:  982
//  z_e = 240.121 cm	r_e: 75		 cm	LinearFiducialCut: 18	cm	Tot mass: 21.2175	tons	R/r: 0.038172	Filling factor: 0.858084	Number of QUPIDs required per array:  589
//  z_e = 200.385 cm	r_e: 82.1	 cm	LinearFiducialCut: 18	cm	Tot mass: 21.0435	tons	R/r: 0.0354645	Filling factor: 0.863586	Number of QUPIDs required per array:  687
//  z_e = 169.755 cm	r_e: 89.2	 cm	LinearFiducialCut: 18	cm	Tot mass: 21.0072	tons	R/r: 0.0331157	Filling factor: 0.865067	Number of QUPIDs required per array:  790
//  z_e = 145.646 cm	r_e: 96.3	 cm	LinearFiducialCut: 18	cm	Tot mass: 21.0838	tons	R/r: 0.0310586	Filling factor: 0.865745	Number of QUPIDs required per array:  898
//  z_e = 126.331 cm	r_e: 103.4	 cm	LinearFiducialCut: 18	cm	Tot mass: 21.2554	tons	R/r: 0.0292422	Filling factor: 0.86837		Number of QUPIDs required per array:  1016

	//==================================== QUPIDs Dimensions =====================================
	m_hGeometryParameters["QUPIDWindowOuterRadius"]				= 37.	*mm;
	m_hGeometryParameters["QUPIDPhotocathodeOuterRadius"]		= 36.4	*mm;
	m_hGeometryParameters["QUPIDPhotocathodeInnerRadius"]		= 36.3	*mm;
	m_hGeometryParameters["QUPIDWindowOuterHeight"]				= 26.57	*mm;
	m_hGeometryParameters["QUPIDPhotocathodeOuterHeight"]		= 25.97	*mm;
	m_hGeometryParameters["QUPIDPhotocathodeInnerHeight"]		= 25.87	*mm;
	m_hGeometryParameters["QUPIDBodyOuterRadius"]				= 35.5	*mm;
	m_hGeometryParameters["QUPIDBodyInnerRadius"]				= 34.2	*mm;
	m_hGeometryParameters["QUPIDBodyHeight"]					= 44.93	*mm;
	m_hGeometryParameters["QUPIDBodyAluminiumCoatingHeight"]	= 15.	*mm;
	m_hGeometryParameters["QUPIDAluminiumCoatingThickness"]		= 0.1	*mm;
	m_hGeometryParameters["QUPIDAPDRadius"]						= 8.	*mm;
	m_hGeometryParameters["QUPIDAPDHeight"]						= 23.	*mm;
	m_hGeometryParameters["QUPIDBaseRadius"]					= 36.	*mm;
	m_hGeometryParameters["QUPIDBaseThickness"]					= 5.	*mm;
	m_hGeometryParameters["QUPIDsMinimumAllowedDistance"] 		= 0.5	*mm;
	m_hGeometryParameters["PhotoSensorsVoltageDividerSpace"] 	= 1.	*cm;

	m_hGeometryParameters["PhotoSensorsHeight"]					= GetGeometryParameter("PhotoSensorsVoltageDividerSpace")+
																  GetGeometryParameter("QUPIDWindowOuterHeight")+
																  GetGeometryParameter("QUPIDBodyHeight")+
																  GetGeometryParameter("QUPIDBaseThickness"); // PS = Photo Sensor
	m_hGeometryParameters["OuterCryostatThickness"]				= 1		*cm;
	m_hGeometryParameters["VacuumThickness"]					= 10	*cm;
	m_hGeometryParameters["InnerCryostatThickness"]				= 1		*cm;
	m_hGeometryParameters["OuterLXeThicknessTop"]				= 2.	*cm;
	m_hGeometryParameters["OuterLXeThickness"]					= 3.	*cm;
	m_hGeometryParameters["OuterLXeHeight"]						= 0.5	*cm;
	m_hGeometryParameters["PTFEThickness"]						= 1		*cm;
	m_hGeometryParameters["BellTopThickness"]					= 0.3	*cm;
	m_hGeometryParameters["BellWidth"]							= 0.3	*cm;
	m_hGeometryParameters["GridRingWidth"]						= GetGeometryParameter("PTFEThickness") - 0.2*cm;
	m_hGeometryParameters["VeryBottomGridToPS"]					= 1.	*cm;
	m_hGeometryParameters["CathodeToVeryBottomGrid"]			= 1.	*cm;
	m_hGeometryParameters["TopGridsHeight"]						= 0.6	*cm;
	m_hGeometryParameters["BottomGridsHeight"]					= 0.8	*cm;
	m_hGeometryParameters["PhotoSensorsToScreeningMesh"]		= 2.	*cm;
	m_hGeometryParameters["ScreeningMeshToAnode"]				= GetGeometryParameter("TopGridsHeight")+
																0.1		*cm;
	m_hGeometryParameters["AnodeToBelowLiquidMesh"]				= GetGeometryParameter("TopGridsHeight")+
																0.1		*cm;
	m_hGeometryParameters["CathodeToVeryBottomMesh"]			= GetGeometryParameter("BottomGridsHeight")+
																0.7		*cm;
	m_hGeometryParameters["VeryBottomMeshToPhotoSensors"]		= GetGeometryParameter("BottomGridsHeight")+
																1.0		*cm;
	m_hGeometryParameters["TopPipeOuterRadius"]					= 10.0	*cm;
	m_hGeometryParameters["TopPipeThickness"]					= 1.5	*cm;
	m_hGeometryParameters["GridMeshThickness"]					= 0.13	*mm;

	//================================== Outer cryostat =================================
	m_hGeometryParameters["OuterCryostatOuterRadius"] 		= GetGeometryParameter("OuterCryostatThickness")+
															GetGeometryParameter("VacuumThickness")+
															GetGeometryParameter("InnerCryostatThickness")+
															GetGeometryParameter("OuterLXeThickness")+
															GetGeometryParameter("PTFEThickness")+
															GetGeometryParameter("LinearFiducialCut")+
															GetGeometryParameter("FiducialRadius");
	m_hGeometryParameters["OuterCryostatHeight"] 			= GetGeometryParameter("VeryBottomMeshToPhotoSensors")+
															GetGeometryParameter("CathodeToVeryBottomMesh")+
															GetGeometryParameter("FiducialDriftLength")+
															2.*GetGeometryParameter("LinearFiducialCut")+
															GetGeometryParameter("AnodeToBelowLiquidMesh")+
															GetGeometryParameter("ScreeningMeshToAnode")+
															GetGeometryParameter("PhotoSensorsToScreeningMesh")+
															2.*GetGeometryParameter("PhotoSensorsHeight")+
															GetGeometryParameter("BellTopThickness")+
															2.*GetGeometryParameter("OuterLXeHeight");
	//================================== Cryostat Vacuum =================================
	m_hGeometryParameters["VacuumOuterRadius"]	= GetGeometryParameter("OuterCryostatOuterRadius")-
												GetGeometryParameter("OuterCryostatThickness");
	m_hGeometryParameters["VacuumHeight"]		= GetGeometryParameter("OuterCryostatHeight");
	//================================== Inner cryostat =================================
	m_hGeometryParameters["InnerCryostatOuterRadius"]	= GetGeometryParameter("VacuumOuterRadius")-
														GetGeometryParameter("VacuumThickness");
	m_hGeometryParameters["InnerCryostatHeight"]		= GetGeometryParameter("OuterCryostatHeight");
	//================================== Outer LXe =================================
	m_hGeometryParameters["OuterLXeOuterRadius"]	= GetGeometryParameter("InnerCryostatOuterRadius")-
													GetGeometryParameter("InnerCryostatThickness");
	m_hGeometryParameters["OuterLXeOuterHeight"]	= GetGeometryParameter("OuterCryostatHeight");
	//================================== Diving Bell =================================
	m_hGeometryParameters["BellHeight"]			= GetGeometryParameter("PhotoSensorsHeight")+
												GetGeometryParameter("PhotoSensorsToScreeningMesh")+
												GetGeometryParameter("ScreeningMeshToAnode")+
												0.5*GetGeometryParameter("AnodeToBeolwLiquidMesh");
												//GetGeometryParameter("OuterLXeHeight");
	m_hGeometryParameters["BellOuterRadius"]	= GetGeometryParameter("OuterLXeOuterRadius") -
												GetGeometryParameter("OuterLXeThicknessTop");;
	m_hGeometryParameters["BellInnerRadius"]	= GetGeometryParameter("BellOuterRadius")-
												GetGeometryParameter("BellWidth");
	//================================== TPC Dimensions =================================
	m_hGeometryParameters["TPCHeight"]		= GetGeometryParameter("OuterLXeOuterHeight")-
											GetGeometryParameter("BellHeight")-
											GetGeometryParameter("BellTopThickness")-
											2*GetGeometryParameter("OuterLXeHeight");
	m_hGeometryParameters["TPCOuterRadius"]	= GetGeometryParameter("OuterLXeOuterRadius")-
											GetGeometryParameter("OuterLXeThickness");
	m_hGeometryParameters["TPCInnerRadius"]	= GetGeometryParameter("LinearFiducialCut")+
											GetGeometryParameter("FiducialRadius");
	//================================== Photo Sensors =================================
	m_hGeometryParameters["PSArrayOuterRadius"]	= GetGeometryParameter("TPCInnerRadius");
	//================================== Grids Dimensions =================================
	m_hGeometryParameters["GridsOuterRadius"]	= GetGeometryParameter("TPCInnerRadius")+
												GetGeometryParameter("GridRingWidth");
	m_hGeometryParameters["GridsInnerRadius"]	= GetGeometryParameter("TPCInnerRadius") + 0.1*mm;
	//================================== Gas Inside the TPC =================================
	m_hGeometryParameters["ObservedGXeHeight"]		= GetGeometryParameter("PhotoSensorsToScreeningMesh")+
													GetGeometryParameter("ScreeningMeshToAnode")+
													0.5*GetGeometryParameter("AnodeToBelowLiquidMesh");
	m_hGeometryParameters["ObservedGXeOuterRadius"]	= GetGeometryParameter("FiducialRadius")+
													GetGeometryParameter("LinearFiducialCut");
	//============================== S2 Dead LXe Layer (between Cathode and Bottom PS Array) =================================
	m_hGeometryParameters["DeadLXeOuterRadius"]	= GetGeometryParameter("FiducialRadius")+
												GetGeometryParameter("LinearFiducialCut");
	m_hGeometryParameters["DeadLXeHeight"]		= GetGeometryParameter("VeryBottomMeshToPhotoSensors")+
												GetGeometryParameter("CathodeToVeryBottomMesh");
	//================================== Sensitive LXe =================================
	m_hGeometryParameters["SensitiveLXeOuterRadius"]	= GetGeometryParameter("FiducialRadius")+
														GetGeometryParameter("LinearFiducialCut");
	m_hGeometryParameters["SensitiveLXeHeight"]	= GetGeometryParameter("FiducialDriftLength")+
												2*GetGeometryParameter("LinearFiducialCut")+
												0.5*GetGeometryParameter("AnodeToBelowLiquidMesh");
	//================================= 10'' PMTs Dimensions ================================== (Hamamatsu R7081MOD-ASSY)
	m_hGeometryParameters["PMTWindowOuterRadius"]		= 126.5*mm;
	m_hGeometryParameters["PMTWindowOuterHalfZ"]		= 89.*mm;
	m_hGeometryParameters["PMTWindowTopZ"]		        = 85.*mm;
	m_hGeometryParameters["PMTPhotocathodeOuterRadius"]	= 125.*mm;
	m_hGeometryParameters["PMTPhotocathodeOuterHalfZ"]	= 87.5*mm;
	m_hGeometryParameters["PMTPhotocathodeTopZ"]		= -43.*mm;
	m_hGeometryParameters["PMTPhotocathodeInnerRadius"]	= 124.5*mm;
	m_hGeometryParameters["PMTPhotocathodeInnerHalfZ"]	= 87.*mm;
	m_hGeometryParameters["PMTBodyOuterRadius"]			= 51.*mm;
	m_hGeometryParameters["PMTBodyInnerRadius"]			= 50.*mm;
	m_hGeometryParameters["PMTBodyHeight"]				= 42.*mm;
	m_hGeometryParameters["PMTBaseOuterRadius"]			= 60.*mm;
	m_hGeometryParameters["PMTBaseInnerRadius"]			= 59.*mm;
	m_hGeometryParameters["PMTBaseHeight"]				= 62.*mm;
	m_hGeometryParameters["PMTBaseInteriorHeight"]		= 60.*mm;
	//==================================== Number of PMTs =====================================
	m_hGeometryParameters["NbPMTs"]						= 121+121+0+73; // Top + Bottom + LS + Water
	m_hGeometryParameters["NbTopPMTs"]					= 121;
	m_hGeometryParameters["NbBottomPMTs"]				= 121;
	m_hGeometryParameters["NbLSPMTs"]					= 0; //101  
	m_hGeometryParameters["NbLSTopPMTs"]				= 0; //13;
	m_hGeometryParameters["NbLSBottomPMTs"]				= 0;//13;
	m_hGeometryParameters["NbLSSidePMTs"]				= 0;//75;
	m_hGeometryParameters["NbLSSidePMTColumns"]			= 0;//15;
	m_hGeometryParameters["NbLSSidePMTRows"]			= 0;//5;
	m_hGeometryParameters["NbWaterPMTs"]				= 73;
	m_hGeometryParameters["NbWaterTopPMTs"]				= 0;
	m_hGeometryParameters["NbWaterBottomPMTs"]			= 25;
	m_hGeometryParameters["NbWaterSidePMTs"]			= 48;
	m_hGeometryParameters["NbWaterSidePMTColumns"]		= 12;
	m_hGeometryParameters["NbWaterSidePMTRows"]			= 4;
	//==================================== PMT positions =====================================
	m_hGeometryParameters["TopQUPIDWindowZ"]			= 58.2*cm;
	m_hGeometryParameters["BottomQUPIDWindowZ"]			= -54.2*cm;
	m_hGeometryParameters["QUPIDDistance"]				= 8.*cm;
	m_hGeometryParameters["LSTopPMTWindowZ"]			= 195.*cm;
	m_hGeometryParameters["LSBottomPMTWindowZ"]			= -195.*cm;
	m_hGeometryParameters["LSSidePMTWindowR"]			= 170.*cm;
	m_hGeometryParameters["LSTopPMTDistance"]			= 80.*cm;
	m_hGeometryParameters["LSBottomPMTDistance"]		= 80.*cm;
	m_hGeometryParameters["LSSidePMTRowDistance"]		= 75.*cm;
	m_hGeometryParameters["WaterTopPMTWindowZ"]			= 440.*cm;
	m_hGeometryParameters["WaterBottomPMTWindowZ"]		= -440.*cm;
	m_hGeometryParameters["WaterSidePMTWindowR"]		= 470.*cm;
	m_hGeometryParameters["WaterTopPMTDistance"]		= 160.*cm; //330
	m_hGeometryParameters["WaterBottomPMTDistance"]		= 160.*cm; //300
	m_hGeometryParameters["WaterSidePMTRowDistance"]	= 250.*cm;

// verifications
//	assert(GetGeometryParameter("OuterCryostatOffsetZ") + GetGeometryParameter("InnerCryostatOffsetZ") + GetGeometryParameter("LXeOffsetZ") == 0);// this way the tpc is centered at 0
}

G4double
DARWINDetectorConstruction::GetGeometryParameter(const char *szParameter)
{
	return m_hGeometryParameters[szParameter];
}

void
DARWINDetectorConstruction::ConstructLaboratory()
{
	//================================== Laboratory =================================
	const G4double dLabHalfZ = 0.5*GetGeometryParameter("LabHeight");
	const G4double dLabRadius = GetGeometryParameter("LabRadius");

	G4Material *Air = G4Material::GetMaterial("G4_AIR");

	G4Tubs *pLabTubs = new G4Tubs("LabTubs", 0.*cm, dLabRadius, dLabHalfZ, 0.*deg, 360.*deg);

	m_pLabLogicalVolume = new G4LogicalVolume(pLabTubs, Air, "LabVolume", 0, 0, 0);

	m_pLabPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(), m_pLabLogicalVolume, "Lab", 0, false, 0);

	m_pMotherLogicalVolume = m_pLabLogicalVolume;

	m_pLabLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
}

void
DARWINDetectorConstruction::ConstructVeto()
{
	//================================== Water tank =================================
	const G4double dWaterTankCylinderHalfZ = 0.5*GetGeometryParameter("WaterTankCylinderHeight");
	const G4double dWaterTankDomeHalfZ = GetGeometryParameter("WaterTankDomeOuterHeight");
	const G4double dWaterTankRadius = GetGeometryParameter("WaterTankOuterRadius");

	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");

	G4Tubs *pWaterTankTubs = new G4Tubs("WaterTankTubs", 0.*cm, dWaterTankRadius, dWaterTankCylinderHalfZ, 0.*deg, 360.*deg);

	G4Ellipsoid *pWaterTankEllipsoid = new G4Ellipsoid("WaterTankEllipsoid", dWaterTankRadius, dWaterTankRadius, dWaterTankDomeHalfZ, 0, dWaterTankDomeHalfZ);

	G4UnionSolid *pWaterTankUnionSolid = new G4UnionSolid("WaterTankUnionSolid", pWaterTankTubs, pWaterTankEllipsoid, 0, G4ThreeVector(0,0,dWaterTankCylinderHalfZ));

	m_pWaterTankLogicalVolume = new G4LogicalVolume(pWaterTankUnionSolid, SS304LSteel, "WaterTankVolume", 0, 0, 0);

	m_pWaterTankPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0, 0, 0),
		m_pWaterTankLogicalVolume, "WaterTank", m_pMotherLogicalVolume, false, 0);

	//================================== Water =================================
	const G4double dWaterCylinderHalfZ = 0.5*GetGeometryParameter("WaterTankCylinderInnerHeight");
	const G4double dWaterDomeHalfZ = GetGeometryParameter("WaterTankDomeInnerHeight");
	const G4double dWaterRadius = GetGeometryParameter("WaterTankInnerRadius");
	const G4double dWaterOffsetZ = 0.5* GetGeometryParameter("WaterTankThickness");

	G4Material *Water = G4Material::GetMaterial("Water");

	G4Tubs *pWaterTubs = new G4Tubs("WaterTubs", 0.*cm, dWaterRadius, dWaterCylinderHalfZ, 0.*deg, 360.*deg);

	G4Ellipsoid *pWaterEllipsoid = new G4Ellipsoid("WaterEllipsoid", dWaterRadius, dWaterRadius, dWaterDomeHalfZ, 0, dWaterDomeHalfZ);

	G4UnionSolid *pWaterUnionSolid = new G4UnionSolid("WaterUnionSolid", pWaterTubs, pWaterEllipsoid, 0, G4ThreeVector(0,0,dWaterCylinderHalfZ));

	m_pWaterLogicalVolume = new G4LogicalVolume(pWaterUnionSolid, Water, "WaterLogicalVolume", 0, 0, 0);

	m_pWaterPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0, 0, dWaterOffsetZ),
		m_pWaterLogicalVolume, "Water", m_pWaterTankLogicalVolume, false, 0);

	//================================== attributes =================================
	//G4Colour hWaterTankColor(0.500, 0.500, 0.500, 0.1);
	G4Colour hWaterTankColor(0.6, 0.6, 0.6, 0.005);
	//G3Colour hWaterTankColor(1.0,1.0,0.0);
	G4VisAttributes *pWaterTankVisAtt = new G4VisAttributes(hWaterTankColor);
	pWaterTankVisAtt->SetVisibility(false);
	m_pWaterTankLogicalVolume->SetVisAttributes(pWaterTankVisAtt);
	//m_pWaterTankLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hWaterColor(0.0, 0.7, 0.7, 0.01);
	G4VisAttributes *pWaterVisAtt = new G4VisAttributes(hWaterColor);
	pWaterVisAtt->SetVisibility(false);
	m_pWaterLogicalVolume->SetVisAttributes(pWaterVisAtt);
	//m_pWaterLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

}

void
DARWINDetectorConstruction::ConstructCryostat()
{
	//G4Material *Titanium = G4Material::GetMaterial("Titanium");
	G4Material *Copper = G4Material::GetMaterial("Copper");
	G4Material *Vacuum = G4Material::GetMaterial("Vacuum");

	const G4double dOuterCryostatOuterRadius = GetGeometryParameter("OuterCryostatOuterRadius");
	const G4double dOuterCryostatHeight = GetGeometryParameter("OuterCryostatHeight");
	const G4double dOuterCryostatDomeOuterRadius = 2*GetGeometryParameter("OuterCryostatOuterRadius");
	//================================== Outer cryostat =================================
	// ******** ROUNDED PARTS ********//
	G4double dOuterCryostatDomeInnerRadius = 0*cm;
	G4Sphere* OuterCryostatDomeTopSphere = new G4Sphere("OuterCryostatDomeTopSphere",		// const G4String& pName,
							dOuterCryostatDomeInnerRadius, dOuterCryostatDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4Sphere* OuterCryostatDomeBottomSphere = new G4Sphere("OuterCryostatDomeBottomSphere",		// const G4String& pName,
							dOuterCryostatDomeInnerRadius, dOuterCryostatDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4double OuterCryostatDomeXOffset = 0.0*cm;
	G4double OuterCryostatDomeYOffset = 0.0*cm;
	G4double OuterCryostatDomeBottomZOffset = -1*dOuterCryostatHeight*0.5001 + sqrt(3)*dOuterCryostatOuterRadius;
	G4double OuterCryostatDomeTopZOffset = -1*(dOuterCryostatOuterRadius*sqrt(3) - dOuterCryostatHeight*0.5001);
//	G4cout<<"Bottom OC dome offset: "<<OuterCryostatDomeBottomZOffset<<G4endl;
//	G4cout<<"Top OC dome offset: "<<OuterCryostatDomeTopZOffset<<G4endl;
//	G4cout<<"OC height: "<<dOuterCryostatHeight<<G4endl;
//	G4cout<<"OC OR: "<<dOuterCryostatOuterRadius<<G4endl;

	G4ThreeVector OuterCryostatDomeBottomZTranslation(OuterCryostatDomeXOffset,OuterCryostatDomeYOffset,OuterCryostatDomeBottomZOffset);	// bottom dome translation
	G4ThreeVector OuterCryostatDomeTopZTranslation(OuterCryostatDomeXOffset,OuterCryostatDomeYOffset,OuterCryostatDomeTopZOffset);	// top dome translation

	// ******** TUBE ******** //
	G4Tubs *OuterCryostatTubs = new G4Tubs("OuterCryostatTubs",0.0*cm,dOuterCryostatOuterRadius,0.5*dOuterCryostatHeight,0*deg,360*deg);

	G4UnionSolid* OuterCryostatUnion1	= new G4UnionSolid("OuterCryostatUnion1",OuterCryostatTubs,OuterCryostatDomeTopSphere,0,OuterCryostatDomeTopZTranslation);
	G4UnionSolid* OuterCryostatUnion	= new G4UnionSolid("OuterCryostatUnion",OuterCryostatUnion1,OuterCryostatDomeBottomSphere,RotP180x,OuterCryostatDomeBottomZTranslation);

	// ******** TOP PIPE ******** //
	const G4double dTopPipeOuterRadius = GetGeometryParameter("TopPipeOuterRadius");
	const G4double dTopPipeThickness = GetGeometryParameter("TopPipeThickness");
	G4double dTopPipeThorusLength = 35.*cm;
	G4double dTopPipeTubeLength = 50*cm;

	// hollow
	G4Torus *TopPipeTorus = new G4Torus("TopPipeTorus",dTopPipeOuterRadius-dTopPipeThickness,dTopPipeOuterRadius,dTopPipeThorusLength,0.*deg,90.*deg);
	G4Tubs *TopPipeTubs = new G4Tubs("TopPipeTubs",dTopPipeOuterRadius-dTopPipeThickness,dTopPipeOuterRadius,dTopPipeTubeLength,0*deg,360*deg);

	G4ThreeVector TranslateTopTubs(-1.0*dTopPipeThorusLength-dTopPipeTubeLength+dTopPipeThorusLength,0.5*dTopPipeTubeLength+dTopPipeOuterRadius,0.0*cm);
	G4UnionSolid* TopPipe   = new G4UnionSolid("TopPipe",TopPipeTorus,TopPipeTubs,RotP90y,TranslateTopTubs);

	G4double TopPipeXOffset = -1.0*dTopPipeThorusLength;
	G4double TopPipeYOffset = 0*cm;
	G4double TopPipeZOffset = OuterCryostatDomeTopZOffset + dOuterCryostatDomeOuterRadius - 2.7*cm;

	G4UnionSolid* OuterCryostatAndPipe = new G4UnionSolid("OuterCryostatAndPipe",OuterCryostatUnion,TopPipe,
			RotP270x,G4ThreeVector(TopPipeXOffset,TopPipeYOffset,TopPipeZOffset));

	G4double OuterCryostatXOffset = 0.0*cm;
	G4double OuterCryostatYOffset = 0.0*cm;
	G4double OuterCryostatZOffset = 0.0*cm;
	//m_pOuterCryostatLogicalVolume = new G4LogicalVolume(OuterCryostatAndPipe, Titanium, "OuterCryostatLogicalVolume");
	m_pOuterCryostatLogicalVolume = new G4LogicalVolume(OuterCryostatAndPipe, Copper, "OuterCryostatLogicalVolume");

	//m_pOuterCryostatLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	m_pOuterCryostatPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(OuterCryostatXOffset, OuterCryostatYOffset, OuterCryostatZOffset),
		m_pOuterCryostatLogicalVolume, "OuterCryostat", m_pWaterLogicalVolume, false, 0);

	// ================================== Insulation Vacuum ====================================
	// ******** ROUNDED PARTS ********
	const G4double dVacuumOuterRadius = GetGeometryParameter("VacuumOuterRadius");
	const G4double dVacuumHeight = GetGeometryParameter("VacuumHeight");
	G4double dVacuumDomeOuterRadius = 2*dVacuumOuterRadius;
	G4double dVacuumDomeInnerRadius = 0*cm;

	G4Sphere* VacuumDomeTopSphere = new G4Sphere("VacuumDomeTopSphere",		// const G4String& pName,
							dVacuumDomeInnerRadius, dVacuumDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4Sphere* VacuumDomeBottomSphere = new G4Sphere("VacuumDomeBottomSphere",		// const G4String& pName,
							dVacuumDomeInnerRadius, dVacuumDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4double VacuumDomeXOffset = 0.0*cm;
	G4double VacuumDomeYOffset = 0.0*cm;
	G4double VacuumDomeBottomZOffset = -1*dVacuumHeight*0.5001 + sqrt(3)*dVacuumOuterRadius;
	G4double VacuumDomeTopZOffset = -1*(dVacuumOuterRadius*sqrt(3) - dVacuumHeight*0.5001);

	G4ThreeVector TranslationVacuumDomeBottom(VacuumDomeXOffset,VacuumDomeYOffset,VacuumDomeBottomZOffset);	// bottom dome translation
	G4ThreeVector TranslationVacuumDomeTop(VacuumDomeXOffset,VacuumDomeYOffset,VacuumDomeTopZOffset);	// top dome translation

	// ******** TUBE ******** //
	G4Tubs *VacuumTubs = new G4Tubs("VacuumTubs",0.0*cm,dVacuumOuterRadius,0.5*dVacuumHeight,0*deg,360*deg);

	G4UnionSolid* VacuumUnion = new G4UnionSolid("VacuumUnion",VacuumTubs,VacuumDomeTopSphere,0,TranslationVacuumDomeTop);
	//G4Transform3D TransformationVacuum(RotP180x,TranslationVacuumDomeBottom);
	G4UnionSolid* VacuumSolid   = new G4UnionSolid("VacuumSolid",VacuumUnion,VacuumDomeBottomSphere,RotP180x,TranslationVacuumDomeBottom);

	G4double VacuumXOffset = 0.0*cm;
	G4double VacuumYOffset = 0.0*cm;
	G4double VacuumZOffset = 0.0*cm;

	m_pCryostatVacuumLogicalVolume = new G4LogicalVolume(VacuumSolid, Vacuum, "CryostatVacuumLogicalVolume");

	m_pCryostatVacuumPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(VacuumXOffset, VacuumYOffset, VacuumZOffset),
			m_pCryostatVacuumLogicalVolume, "CryostatVacuum", m_pOuterCryostatLogicalVolume, false, 0);

	// ================================== Inner Cryostat ====================================
	// ******** ROUNDED PARTS ********

	G4double dInnerCryostatOuterRadius = GetGeometryParameter("InnerCryostatOuterRadius");
	G4double dInnerCryostatHeight = GetGeometryParameter("InnerCryostatHeight");
	G4double dInnerCryostatDomeOuterRadius = 2*dInnerCryostatOuterRadius;
	G4double dInnerCryostatDomeInnerRadius = 0*cm;
	G4Sphere* InnerCryostatDomeTopSphere = new G4Sphere("InnerCryostatDomeTopSphere",		// const G4String& pName,
							dInnerCryostatDomeInnerRadius, dInnerCryostatDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4Sphere* InnerCryostatDomeBottomSphere = new G4Sphere("InnerCryostatDomeBottomSphere",		// const G4String& pName,
							dInnerCryostatDomeInnerRadius, dInnerCryostatDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4double InnerCryostatDomeXOffset = 0.0*cm;
	G4double InnerCryostatDomeYOffset = 0.0*cm;
	G4double InnerCryostatDomeBottomZOffset = -1*dInnerCryostatHeight*0.5001 + sqrt(3)*dInnerCryostatOuterRadius;
	G4double InnerCryostatDomeTopZOffset = -1*(dInnerCryostatOuterRadius*sqrt(3) - dInnerCryostatHeight*0.5001);

	// ************ TUBE ************//
	G4ThreeVector TranslationInnerCryostatDomeBottom(InnerCryostatDomeXOffset,InnerCryostatDomeYOffset,InnerCryostatDomeBottomZOffset);	// bottom dome translation
	G4ThreeVector TranslationInnerCryostatDomeTop(InnerCryostatDomeXOffset,InnerCryostatDomeYOffset,InnerCryostatDomeTopZOffset);	// top dome translation

	G4Tubs *InnerCryostatTubs = new G4Tubs("InnerCryostatTubs",0.0*cm,dInnerCryostatOuterRadius,0.5*dInnerCryostatHeight,0*deg,360*deg);
	G4UnionSolid* InnerCryostatUnion1	= new G4UnionSolid("InnerCryostatUnion1",InnerCryostatTubs,InnerCryostatDomeTopSphere,0,TranslationInnerCryostatDomeTop);
	G4UnionSolid* InnerCryostatUnion  	= new G4UnionSolid("InnerCryostatUnion",InnerCryostatUnion1,InnerCryostatDomeBottomSphere,RotP180x,TranslationInnerCryostatDomeBottom);
	//m_pInnerCryostatLogicalVolume = new G4LogicalVolume(InnerCryostatUnion, Titanium, "InnerCryostatLogicalVolume");
	m_pInnerCryostatLogicalVolume = new G4LogicalVolume(InnerCryostatUnion, Copper, "InnerCryostatLogicalVolume");

	G4double InnerCryostatXOffset = 0.0*cm;
	G4double InnerCryostatYOffset = 0.0*cm;
	G4double InnerCryostatZOffset = 0.0*cm;


	m_pInnerCryostatPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(InnerCryostatXOffset, InnerCryostatYOffset, InnerCryostatZOffset),
		"InnerCryostat", m_pInnerCryostatLogicalVolume, m_pCryostatVacuumPhysicalVolume, false, 0);

	//================================== attributes =================================
	G4Colour hTitaniumColor(0.600, 0.600, 0.600, 0.1);
	G4VisAttributes *pTitaniumVisAtt = new G4VisAttributes(hTitaniumColor);
	pTitaniumVisAtt->SetVisibility(true);
	m_pOuterCryostatLogicalVolume->SetVisAttributes(pTitaniumVisAtt);
	//m_pOuterCryostatLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pCryostatVacuumLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pInnerCryostatLogicalVolume->SetVisAttributes(pTitaniumVisAtt);
	//m_pInnerCryostatLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);


	G4Colour hVacuumColor(1.0, 0.0, 0.0);
	G4VisAttributes *pVacuumVisAtt = new G4VisAttributes(hVacuumColor);
	pVacuumVisAtt->SetVisibility(true);
	m_pCryostatVacuumLogicalVolume->SetVisAttributes(pVacuumVisAtt);
	m_pInnerCryostatLogicalVolume->SetVisAttributes(pTitaniumVisAtt);

	//---------------------------------- optical surfaces -------------------------------
	//--------- between Water and OuterCryostat
	G4OpticalSurface* OpInnerWaterCryoSurface = new G4OpticalSurface("InnerWaterCryoSurface");
	OpInnerWaterCryoSurface->SetType(dielectric_metal);
	OpInnerWaterCryoSurface->SetFinish(ground);
	OpInnerWaterCryoSurface->SetModel(unified);

	G4LogicalBorderSurface* InnerWaterCryoSurface = 
		new G4LogicalBorderSurface("InnerWaterCryoSurface",m_pWaterPhysicalVolume,m_pOuterCryostatPhysicalVolume,OpInnerWaterCryoSurface);
	
	if(InnerWaterCryoSurface->GetVolume1() == m_pWaterPhysicalVolume) G4cout << "Vol 1 equal to Water" << G4endl;
	if(InnerWaterCryoSurface->GetVolume2() == m_pOuterCryostatPhysicalVolume ) G4cout << "Vol 2 equal to OuterCryostat" << G4endl;

}

void
DARWINDetectorConstruction::ConstructXenon()
{
	G4Material *GXe = G4Material::GetMaterial("GXe");
	G4Material *LXe = G4Material::GetMaterial("LXe");

	const G4double dOuterLXeOuterRadius = GetGeometryParameter("OuterLXeOuterRadius");
	const G4double dOuterLXeOuterHeight = GetGeometryParameter("OuterLXeOuterHeight");
	const G4double dOuterLXeHeight = GetGeometryParameter("OuterLXeHeight");
	const G4double dBellInnerRadius = GetGeometryParameter("BellInnerRadius");
	const G4double dBellTopThickness = GetGeometryParameter("BellTopThickness");
	const G4double dBellHeight = GetGeometryParameter("BellHeight");

	// ========================= Top GXe ==============================
	// ******** ROUNDED PARTS ********//
	G4double dGXeDomeOuterRadius = 2*dOuterLXeOuterRadius;
	G4double dGXeDomeInnerRadius = 0*cm;
	G4Sphere* GXeDomeTopSphere = new G4Sphere("GXeDomeTopSphere",		// const G4String& pName,
							dGXeDomeInnerRadius, dGXeDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4Sphere* GXeDomeBottomSphere = new G4Sphere("GXeDomeBottomSphere",		// const G4String& pName,
							dGXeDomeInnerRadius, dGXeDomeOuterRadius,		// G4double pRmin, G4double pRmax,
							0*deg, 360*deg,				// G4double pSPhi, G4double pDPhi,
							0*deg, 30*deg);				// G4double pSTheta, G4double pDTheta);

	G4double dGXeDomeXOffset = 0.0*cm;
	G4double dGXeDomeYOffset = 0.0*cm;
	G4double dGXeDomeBottomZOffset = -1*dOuterLXeOuterHeight*0.5001 + sqrt(3)*dOuterLXeOuterRadius;
	G4double dGXeDomeTopZOffset = -1*(dOuterLXeOuterRadius*sqrt(3) - dOuterLXeOuterHeight*0.5001);

	G4ThreeVector GXeDomeBottomTranslation(dGXeDomeXOffset,dGXeDomeYOffset,dGXeDomeBottomZOffset);	// bottom dome translation
	G4ThreeVector GXeDomeTopTranslation(dGXeDomeXOffset,dGXeDomeYOffset,dGXeDomeTopZOffset);	// top dome translation

	// ************ TUBE ************//
	G4Tubs *GXeTubs = new G4Tubs("GXeTubs",0.0*cm,dOuterLXeOuterRadius,0.5*dOuterLXeOuterHeight,0.0*deg,360.0*deg);

	G4UnionSolid* GXeUnion1 = new G4UnionSolid("GXeUnion1",GXeTubs,GXeDomeTopSphere,0,GXeDomeTopTranslation);
	G4UnionSolid* GXeUnion = new G4UnionSolid("GXeUnion",GXeUnion1,GXeDomeBottomSphere,RotP180x,GXeDomeBottomTranslation);
	m_pGXeLogicalVolume = new G4LogicalVolume(GXeUnion, GXe, "GXeLogicalVolume");

	G4double GXeXOffset = 0.0*cm;
	G4double GXeYOffset = 0.0*cm;
	G4double GXeZOffset = 0.0*cm;

	m_pGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(GXeXOffset, GXeYOffset, GXeZOffset),"GXePhysicalVolume", m_pGXeLogicalVolume, m_pInnerCryostatPhysicalVolume, false, 0);

	// =========================== Outer LXe cylinder ============================
	// ******** ROUNDED PARTS FROM GXe ********//
	// ************ TUBE ************//
	G4Tubs *OuterLXeTubs = new G4Tubs("OuterLXeTubs",0.0*cm,dOuterLXeOuterRadius,0.5*dOuterLXeOuterHeight,0.0*deg,360.0*deg);
	G4UnionSolid* OuterLXeUnion   = new G4UnionSolid("OuterLXeUnion",OuterLXeTubs,GXeDomeBottomSphere,RotP180x,GXeDomeBottomTranslation);

	m_pOuterLXeLogicalVolume = new G4LogicalVolume(OuterLXeUnion, LXe, "m_pOuterLXeLogicalVolume");

	G4double OuterLXeXOffset = 0.0*cm;
	G4double OuterLXeYOffset = 0.0*cm;
	G4double OuterLXeZOffset = 0.0*cm;

	m_pOuterLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(OuterLXeXOffset, OuterLXeYOffset, OuterLXeZOffset),
			"m_pOuterLXePhysicalVolume", m_pOuterLXeLogicalVolume, m_pGXePhysicalVolume, false, 0);

	// =============================== Top Gas inside the Diving Bell ====================================
	G4double dInnerGasHeight = dBellHeight;
	G4Tubs *InnerGXeTubs = new G4Tubs("InnerGXeTubs",0.0*cm,dBellInnerRadius,0.5*(dInnerGasHeight)-0.1*mm,0.0*deg,360.0*deg);

	m_pInnerGXeLogicalVolume = new G4LogicalVolume(InnerGXeTubs,GXe,"InnerGXeLogicalVolume");
	G4double InnerGXeXOffset = 0.0*cm;
	G4double InnerGXeYOffset = 0.0*cm;
	G4double InnerGXeZOffset = dOuterLXeOuterHeight*0.5 - (dInnerGasHeight)*0.5 - dBellTopThickness - dOuterLXeHeight;

	m_pInnerGXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(InnerGXeXOffset, InnerGXeYOffset, InnerGXeZOffset),
			"InnerGXePhysicalVolume", m_pInnerGXeLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	//============================== xenon sensitivity ==============================
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	DARWINLXeSensitiveDetector *pLXeSD = new DARWINLXeSensitiveDetector("DARWIN/LXeSD");
	pSDManager->AddNewDetector(pLXeSD);
	m_pOuterLXeLogicalVolume->SetSensitiveDetector(pLXeSD);
	m_pInnerGXeLogicalVolume->SetSensitiveDetector(pLXeSD);

	//================================== attributes =================================
	//G4Colour hLXeColor(0.094, 0.718, 0.812, 0.05);
	G4Colour hLXeColor(1.0, 1.0, 1.0);
	G4VisAttributes *pLXeVisAtt = new G4VisAttributes(hLXeColor);
	pLXeVisAtt->SetVisibility(true);
	m_pOuterLXeLogicalVolume->SetVisAttributes(pLXeVisAtt);
	//m_pOuterLXeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hGXeColor(0.239, 0.918, 0.878, 0.01);
	G4VisAttributes *pGXeVisAtt = new G4VisAttributes(hGXeColor);
	pGXeVisAtt->SetVisibility(true);
	//m_pGXeLogicalVolume->SetVisAttributes(pGXeVisAtt);
	//m_pInnerGXeLogicalVolume->SetVisAttributes(pGXeVisAtt);
	m_pGXeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pInnerGXeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

}

void
DARWINDetectorConstruction::ConstructTPC()
{
	G4Material *Teflon 					= G4Material::GetMaterial("Teflon");
	//G4Material *Titanium 				= G4Material::GetMaterial("Titanium");
	G4Material *Copper 					= G4Material::GetMaterial("Copper");

	const G4double dBellOuterRadius 	= GetGeometryParameter("BellOuterRadius");
	const G4double dBellInnerRadius 	= GetGeometryParameter("BellInnerRadius");
	const G4double dBellHeight 			= GetGeometryParameter("BellHeight");
	const G4double dBellTopThickness 	= GetGeometryParameter("BellTopThickness");
	const G4double dTPCOuterRadius 		= GetGeometryParameter("TPCOuterRadius");
	const G4double dTPCInnerRadius 		= GetGeometryParameter("TPCInnerRadius");
	const G4double dTPCHeight 			= GetGeometryParameter("TPCHeight");
	const G4double dOuterLXeOuterHeight = GetGeometryParameter("OuterLXeOuterHeight");
	const G4double dOuterLXeHeight 		= GetGeometryParameter("OuterLXeHeight");

	// =============================== Diving Bell ====================================
	G4Tubs *BellTubs = new G4Tubs("BellTubs",dBellInnerRadius,dBellOuterRadius,0.5*dBellHeight,0.0*deg,360.0*deg);
	G4Tubs *BellDisk = new G4Tubs("BellDisk",0.0*cm,dBellOuterRadius,0.5*dBellTopThickness,0.0*deg,360.0*deg);

	G4double BellDiskZOffset = 0.5*dBellTopThickness+0.5*dBellHeight;
	G4ThreeVector BellTopTranslation(0.0,0.0,BellDiskZOffset);	// top disk translation

	G4UnionSolid* BellSolid = new G4UnionSolid("BellSolid",BellTubs,BellDisk,0,BellTopTranslation);

	//m_pBellLogicalVolume = new G4LogicalVolume(BellSolid,Titanium,"BellLogicalVolume");
	m_pBellLogicalVolume = new G4LogicalVolume(BellSolid, Copper,"BellLogicalVolume");
	G4double BellXOffset = 0.0*cm;
	G4double BellYOffset = 0.0*cm;
	G4double BellZOffset = dOuterLXeOuterHeight*0.5 - (dBellHeight+2*dBellTopThickness)*0.5 - dOuterLXeHeight;

	G4Colour grey (0.5, 0.5, 0.5);
	G4VisAttributes* m_pBellVisibility = new G4VisAttributes(grey);
	m_pBellVisibility->SetVisibility(true);
	m_pBellVisibility->SetForceSolid(true);
	m_pBellLogicalVolume->SetVisAttributes(m_pBellVisibility);
	//m_pBellLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	m_pBellPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(BellXOffset, BellYOffset, BellZOffset),
			"Bell", m_pBellLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	// =================================== TPC ==========================================
	G4Tubs *TPCTubs = new G4Tubs("TPCTubs",dTPCInnerRadius,dTPCOuterRadius,dTPCHeight*0.5,0.0*deg,360.0*deg);
	m_pTPCLogicalVolume = new G4LogicalVolume(TPCTubs,Teflon,"TPCLogicalVolume");
	G4double TPCXOffset = 0.0*cm;
	G4double TPCYOffset = 0.0*cm;
	G4double TPCZOffset = -0.5*(dBellHeight+dBellTopThickness);

	G4cout<<dTPCHeight<<"\t"<<dOuterLXeOuterHeight<<"\t"<<dBellHeight<<"\t"<<dOuterLXeHeight<<G4endl;

	m_pTPCPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(TPCXOffset, TPCYOffset, TPCZOffset),
			"TPC", m_pTPCLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	G4Colour hPTFEColor(1., 1., 1., 0.05);
	G4VisAttributes *pPTFEVisAtt = new G4VisAttributes(hPTFEColor);
	pPTFEVisAtt->SetVisibility(true);
	m_pTPCLogicalVolume->SetVisAttributes(pPTFEVisAtt);
	//m_pTPCLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	//================================== optical surfaces =================================	
	G4double dSigmaAlpha = 0.1;
	G4OpticalSurface *pTeflonOpticalSurface = new G4OpticalSurface("TeflonOpticalSurface",
		unified, groundbackpainted, dielectric_dielectric, dSigmaAlpha);
	pTeflonOpticalSurface->SetMaterialPropertiesTable(Teflon->GetMaterialPropertiesTable());	
	
	new G4LogicalBorderSurface("SidePTFELogicalBorderSurface",
		m_pOuterLXePhysicalVolume, m_pTPCPhysicalVolume, pTeflonOpticalSurface);

}

void
DARWINDetectorConstruction::ConstructFieldCage()
{
	G4Material *Titanium = G4Material::GetMaterial("Titanium");

	G4Material *GridMeshAluminium = G4Material::GetMaterial("GridMeshAluminium");

	const G4double dTPCHeight = GetGeometryParameter("TPCHeight");
	const G4double dTPCInnerRadius = GetGeometryParameter("TPCInnerRadius");
	const G4double dGridsOuterRadius = GetGeometryParameter("GridsOuterRadius");
	const G4double dGridsInnerRadius = GetGeometryParameter("GridsInnerRadius");
	const G4double dTopGridsHeight = GetGeometryParameter("TopGridsHeight");
	const G4double dBottomGridsHeight = GetGeometryParameter("BottomGridsHeight");
	const G4double dGridMeshHalfZ = 0.5*GetGeometryParameter("GridMeshThickness");
	const G4double dPhotoSensorsHeight = GetGeometryParameter("PhotoSensorsHeight");
	const G4double dScreeningMeshToAnode = GetGeometryParameter("ScreeningMeshToAnode");
	const G4double dAnodeToBelowLiquidMesh = GetGeometryParameter("AnodeToBelowLiquidMesh");
	const G4double dCathodeToVeryBottomMesh = GetGeometryParameter("CathodeToVeryBottomMesh");
	const G4double dVeryBottomMeshToPhotoSensors = GetGeometryParameter("VeryBottomMeshToPhotoSensors");
	const G4double dBellHeight = GetGeometryParameter("BellHeight");
	const G4double dBellTopThickness = GetGeometryParameter("BellTopThickness");
	const G4double dOuterLXeOuterHeight = GetGeometryParameter("OuterLXeOuterHeight");
	const G4double dOuterLXeHeight = GetGeometryParameter("OuterLXeHeight");

	G4double TPCZOffset = -.5*(dBellHeight+dBellTopThickness);

	// ===================================== Grids =======================================
	G4Tubs *TopGridRingTubs = new G4Tubs("TopGridRingTubs",dGridsInnerRadius,dGridsOuterRadius,dTopGridsHeight*0.5,0.0*deg,360.0*deg);
	G4Tubs *BottomGridRingTubs = new G4Tubs("BottomGridRingTubs",dGridsInnerRadius,dGridsOuterRadius,dBottomGridsHeight*0.5,0.0*deg,360.0*deg);

	G4Tubs *GridMeshTubs = new G4Tubs("GridMeshTubs", 0., dTPCInnerRadius-0.1*mm, dGridMeshHalfZ,0.0*deg,360.0*deg);

	m_pTopGridsRingLogicalVolume = new G4LogicalVolume(TopGridRingTubs,Titanium,"TopGridsRingLogicalVolume");
	m_pBottomGridsRingLogicalVolume = new G4LogicalVolume(BottomGridRingTubs,Titanium,"BottomGridsRingLogicalVolume");

	m_pGridMeshLogicalVolume = new G4LogicalVolume(GridMeshTubs, GridMeshAluminium, "GridMeshLogicalVolume");

	G4double dGridsXOffset = 0.0*cm;
	G4double dGridsYOffset = 0.0*cm;

	//================================== Top grid ring + mesh =================================
	G4double dTopScreeningGridZOffset = -0.5*dBellHeight+dScreeningMeshToAnode+0.5*dAnodeToBelowLiquidMesh;
	m_pTopScreeningGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dTopScreeningGridZOffset),
		"TopScreeningGridRing", m_pTopGridsRingLogicalVolume, m_pInnerGXePhysicalVolume, false, 0);

	m_pTopGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dTopScreeningGridZOffset),
		"TopGridMesh", m_pGridMeshLogicalVolume, m_pInnerGXePhysicalVolume, false, 0);

	//================================== Anode grid ring + mesh =================================
	G4double dAnodeGridZOffset = dTopScreeningGridZOffset-dScreeningMeshToAnode;
	m_pAnodeGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dAnodeGridZOffset),
		"AnodeGridRing", m_pTopGridsRingLogicalVolume, m_pInnerGXePhysicalVolume, false, 0);

	m_pAnodeGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dAnodeGridZOffset),
		"AnodeMesh", m_pGridMeshLogicalVolume, m_pInnerGXePhysicalVolume, false, 0);

	//================================== Below liquid grid ring + mesh =================================
	G4double dBelowLiquidGridZOffset = 0.5*dTPCHeight-0.5*dAnodeToBelowLiquidMesh;
	m_pBelowLiquidGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dBelowLiquidGridZOffset),
		"BelowLiquidGridRing", m_pTopGridsRingLogicalVolume, m_pTPCPhysicalVolume, false, 0);

	G4double dBelowLiquidGridMeshZOffset = 0.5*dTPCHeight+TPCZOffset-0.5*dAnodeToBelowLiquidMesh;
	m_pBelowLiquidGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dBelowLiquidGridMeshZOffset),
		"BelowLiquidGridMesh", m_pGridMeshLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	//================================== Cathode grid ring + mesh =================================
	G4double dCathodeGridZOffset = -0.5*dTPCHeight +
								dPhotoSensorsHeight +
								dCathodeToVeryBottomMesh +
								dVeryBottomMeshToPhotoSensors;
	m_pCathodeGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dCathodeGridZOffset),
		"CathodeGridRing", m_pBottomGridsRingLogicalVolume, m_pTPCPhysicalVolume, false, 0);

	G4double dCathodeGridMeshZOffset = -0.5*dOuterLXeOuterHeight+
										dOuterLXeHeight+
										dPhotoSensorsHeight+
										dCathodeToVeryBottomMesh+
										dVeryBottomMeshToPhotoSensors;
	m_pCathodeGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dCathodeGridMeshZOffset),
		"CathodeMesh", m_pGridMeshLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	//================================== Very Bottom grid ring + mesh =================================
	G4double dVeryBottomGridZOffset = -0.5*dTPCHeight+
									dPhotoSensorsHeight+
									dVeryBottomMeshToPhotoSensors;
	m_pVeryBottomGridRingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dVeryBottomGridZOffset),
		"VeryBottomGridRing", m_pBottomGridsRingLogicalVolume, m_pTPCPhysicalVolume, false, 0);

	G4double dVeryBottomGridMeshZOffset = -0.5*dOuterLXeOuterHeight+dOuterLXeHeight+dPhotoSensorsHeight+dVeryBottomMeshToPhotoSensors;
	m_pVeryBottomGridMeshPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(dGridsXOffset, dGridsYOffset, dVeryBottomGridMeshZOffset),
		"VeryBottomGridMesh", m_pGridMeshLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	G4Colour hGridMeshColor(0.42, 0.72, 0.0, 1.);
	G4VisAttributes *pGridMeshVisAtt = new G4VisAttributes(hGridMeshColor);
	pGridMeshVisAtt->SetVisibility(true);
	//m_pTopGridsRingLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pBottomGridsRingLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pGridMeshLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pTopGridsRingLogicalVolume	->SetVisAttributes(pGridMeshVisAtt);
	m_pBottomGridsRingLogicalVolume	->SetVisAttributes(pGridMeshVisAtt);
	m_pGridMeshLogicalVolume		->SetVisAttributes(pGridMeshVisAtt);

}

void
DARWINDetectorConstruction::ConstructQUPIDArrays()
{
	G4Material *Quartz = G4Material::GetMaterial("Quartz");
	G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
	G4Material *PhotoCathodeAluminium = G4Material::GetMaterial("PhotoCathodeAluminium");
	G4Material *CoatingAluminium = G4Material::GetMaterial("CoatingAluminium");

	//==================================== QUPIDs =====================================
 	const G4double dQUPIDWindowOuterRadius = GetGeometryParameter("QUPIDWindowOuterRadius");
 	const G4double dQUPIDPhotocathodeOuterRadius = GetGeometryParameter("QUPIDPhotocathodeOuterRadius");
 	const G4double dQUPIDPhotocathodeInnerRadius = GetGeometryParameter("QUPIDPhotocathodeInnerRadius");
 	const G4double dQUPIDWindowOuterHeight = GetGeometryParameter("QUPIDWindowOuterHeight");
 	const G4double dQUPIDPhotocathodeOuterHeight = GetGeometryParameter("QUPIDPhotocathodeOuterHeight");
 	const G4double dQUPIDPhotocathodeInnerHeight = GetGeometryParameter("QUPIDPhotocathodeInnerHeight");

 	const G4double dQUPIDWindowZCut = dQUPIDWindowOuterRadius - dQUPIDWindowOuterHeight;
 	const G4double dQUPIDPhotocathodeOuterZCut = dQUPIDPhotocathodeOuterRadius - dQUPIDPhotocathodeOuterHeight;
 	const G4double dQUPIDPhotocathodeInnerZCut = dQUPIDPhotocathodeInnerRadius - dQUPIDPhotocathodeInnerHeight;

	const G4double dQUPIDBodyOuterRadius = GetGeometryParameter("QUPIDBodyOuterRadius");
	const G4double dQUPIDBodyInnerRadius = GetGeometryParameter("QUPIDBodyInnerRadius");
	const G4double dQUPIDBodyHeight = GetGeometryParameter("QUPIDBodyHeight");

	const G4double dQUPIDBodyAluminiumCoatingHeight = GetGeometryParameter("QUPIDBodyAluminiumCoatingHeight");
	const G4double dQUPIDBodyAluminiumCoatingOffsetZ = 0.5*(dQUPIDBodyHeight-dQUPIDBodyAluminiumCoatingHeight);
	const G4double dQUPIDAluminiumCoatingThickness = GetGeometryParameter("QUPIDAluminiumCoatingThickness");

	const G4double dQUPIDAPDRadius = GetGeometryParameter("QUPIDAPDRadius");
	const G4double dQUPIDAPDHeight = GetGeometryParameter("QUPIDAPDHeight");
	const G4double dQUPIDAPDOffsetZ = -0.5*(dQUPIDBodyHeight-dQUPIDAPDHeight);

	const G4double dQUPIDBaseRadius = GetGeometryParameter("QUPIDBaseRadius");
	const G4double dQUPIDBaseThickness = GetGeometryParameter("QUPIDBaseThickness");
	const G4double dQUPIDBaseAluminiumCoatingOffsetZ = 0.5*(dQUPIDBaseThickness-dQUPIDAluminiumCoatingThickness);

	const G4double dQUPIDsMinimumDisplacement = GetGeometryParameter("QUPIDsMinimumAllowedDistance");

	//const G4double dTPCInnerRadius = GetGeometryParameter("TPCInnerRadius");
	const G4double dBellHeight = GetGeometryParameter("BellHeight");
	const G4double dPSArrayOuterRadius = GetGeometryParameter("PSArrayOuterRadius");

	const G4double dOuterLXeHeight = GetGeometryParameter("OuterLXeHeight");
	const G4double dOuterLXeOuterHeight = GetGeometryParameter("OuterLXeOuterHeight");
	const G4double dPhotoSensorsHeight = GetGeometryParameter("PhotoSensorsHeight");

	G4cout<<"dQUPIDWindowZCut: "<<dQUPIDWindowZCut<<G4endl;
	//--------------------------------- QUPID window ----------------------------------
	G4Ellipsoid *pQUPIDWindowEllipsoid = new G4Ellipsoid("QUPIDWindowEllipsoid",
		dQUPIDWindowOuterRadius, dQUPIDWindowOuterRadius, dQUPIDWindowOuterRadius, dQUPIDWindowZCut, dQUPIDWindowOuterRadius);

	m_pQUPIDWindowLogicalVolume = new G4LogicalVolume(pQUPIDWindowEllipsoid, Quartz, "QUPIDWindowVolume", 0, 0, 0);

	//--------------------------------- QUPID photocathode ----------------------------------
	G4Ellipsoid *pQUPIDPhotocathodeEllipsoid = new G4Ellipsoid("QUPIDPhotocathodeEllipsoid", dQUPIDPhotocathodeOuterRadius,
		dQUPIDPhotocathodeOuterRadius, dQUPIDPhotocathodeOuterRadius, dQUPIDPhotocathodeOuterZCut, dQUPIDPhotocathodeOuterRadius);

	m_pQUPIDPhotocathodeLogicalVolume = new G4LogicalVolume(pQUPIDPhotocathodeEllipsoid, PhotoCathodeAluminium, "QUPIDPhotocathodeVolume", 0, 0, 0);

	m_pQUPIDPhotocathodePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pQUPIDPhotocathodeLogicalVolume, "QUPIDPhotocathode", m_pQUPIDWindowLogicalVolume, false, 0);

	//--------------------------------- QUPID photocathode interior----------------------------------
	G4Ellipsoid *pQUPIDPhotocathodeInteriorEllipsoid = new G4Ellipsoid("QUPIDPhotocathodeInteriorEllipsoid", dQUPIDPhotocathodeInnerRadius,
		dQUPIDPhotocathodeInnerRadius, dQUPIDPhotocathodeInnerRadius, dQUPIDPhotocathodeInnerZCut, dQUPIDPhotocathodeInnerRadius);

	m_pQUPIDPhotocathodeInteriorLogicalVolume = new G4LogicalVolume(pQUPIDPhotocathodeInteriorEllipsoid,
		Vacuum, "QUPIDPhotocathodeInteriorVolume", 0, 0, 0);

	m_pQUPIDPhotocathodeInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pQUPIDPhotocathodeInteriorLogicalVolume, "QUPIDPhotocathodeInterior", m_pQUPIDPhotocathodeLogicalVolume, false, 0);

	//--------------------------------- QUPID body ----------------------------------
	const G4double dQUPIDBodyHalfZ = 0.5*dQUPIDBodyHeight;

	G4Tubs *pQUPIDBodyTubs = new G4Tubs("QUPIDBodyTubs", 0, dQUPIDBodyOuterRadius, dQUPIDBodyHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDBodyLogicalVolume = new G4LogicalVolume(pQUPIDBodyTubs, Quartz, "QUPIDBodyVolume", 0, 0, 0);

	//--------------------------------- QUPID body Aluminium coating ----------------------------------
	const G4double dQUPIDBodyAluminiumCoatingHalfZ = 0.5*dQUPIDBodyAluminiumCoatingHeight;

	G4Tubs *pQUPIDBodyAluminiumCoatingTubs = new G4Tubs("QUPIDBodyAluminiumCoatingTubs", dQUPIDBodyInnerRadius,
		dQUPIDBodyInnerRadius+dQUPIDAluminiumCoatingThickness, dQUPIDBodyAluminiumCoatingHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDBodyAluminiumCoatingLogicalVolume = new G4LogicalVolume(pQUPIDBodyAluminiumCoatingTubs,
		CoatingAluminium, "QUPIDBodyAluminiumCoatingVolume", 0, 0, 0);

	m_pQUPIDBodyAluminiumCoatingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dQUPIDBodyAluminiumCoatingOffsetZ),
		m_pQUPIDBodyAluminiumCoatingLogicalVolume, "QUPIDBodyAluminiumCoating", m_pQUPIDBodyLogicalVolume, false, 0);

	//--------------------------------- QUPID body interior ----------------------------------
	G4Tubs *pQUPIDBodyInteriorTubs = new G4Tubs("QUPIDBodyInteriorTubs", 0, dQUPIDBodyInnerRadius, dQUPIDBodyHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDBodyInteriorLogicalVolume = new G4LogicalVolume(pQUPIDBodyInteriorTubs, Vacuum, "QUPIDBodyInteriorVolume", 0, 0, 0);

	m_pQUPIDBodyInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pQUPIDBodyInteriorLogicalVolume, "QUPIDBodyInterior", m_pQUPIDBodyLogicalVolume, false, 0);

	//--------------------------------- QUPID APD ----------------------------------
	const G4double dQUPIDAPDHalfZ = 0.5*dQUPIDAPDHeight;

	G4Tubs *pQUPIDAPDTubs = new G4Tubs("QUPIDAPDTubs", 0, dQUPIDAPDRadius, dQUPIDAPDHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDAPDLogicalVolume = new G4LogicalVolume(pQUPIDAPDTubs, Quartz, "QUPIDAPDVolume", 0, 0, 0);

	m_pQUPIDAPDPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dQUPIDAPDOffsetZ),
		m_pQUPIDAPDLogicalVolume, "QUPIDAPD", m_pQUPIDBodyInteriorLogicalVolume, false, 0);

	//--------------------------------- QUPID base ----------------------------------
	const G4double dQUPIDBaseHalfZ = 0.5*dQUPIDBaseThickness;

	G4Tubs *pQUPIDBaseTubs = new G4Tubs("QUPIDBaseTubs", 0, dQUPIDBaseRadius, dQUPIDBaseHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDBaseLogicalVolume = new G4LogicalVolume(pQUPIDBaseTubs, Quartz, "QUPIDBaseVolume", 0, 0, 0);

	//--------------------------------- QUPID base Aluminium coating ----------------------------------
	const G4double dQUPIDBaseAluminiumCoatingHalfZ = 0.5*dQUPIDAluminiumCoatingThickness;

	G4Tubs *pQUPIDBaseAluminiumCoatingTubs = new G4Tubs("QUPIDBaseAluminiumCoatingTubs",
		0, dQUPIDBodyInnerRadius, dQUPIDBaseAluminiumCoatingHalfZ, 0.*deg, 360.*deg);

	m_pQUPIDBaseAluminiumCoatingLogicalVolume = new G4LogicalVolume(pQUPIDBaseAluminiumCoatingTubs,
		CoatingAluminium, "QUPIDBaseAluminiumCoatingVolume", 0, 0, 0);

	m_pQUPIDBaseAluminiumCoatingPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., dQUPIDBaseAluminiumCoatingOffsetZ),
		m_pQUPIDBaseAluminiumCoatingLogicalVolume, "QUPIDBaseAluminiumCoating", m_pQUPIDBaseLogicalVolume, false, 0);

	//================================== Define Offsets ==================================
	stringstream hVolumeName;

	G4double TopPhotoSensorsZOffset = 0.5*dBellHeight-dPhotoSensorsHeight; // lower ooffset
	G4double BottomPhotoSensorsZOffset = -0.5*dOuterLXeOuterHeight+dOuterLXeHeight+dPhotoSensorsHeight;
	G4double dQUPIDTopWindowZOffset = TopPhotoSensorsZOffset+dQUPIDWindowOuterHeight;//-dQUPIDWindowZCut;
	G4double dQUPIDTopBodyZOffset = dQUPIDTopWindowZOffset+0.5*dQUPIDBodyHeight-dQUPIDWindowZCut;
	G4double dQUPIDTopBaseZOffset = dQUPIDTopBodyZOffset+0.5*(dQUPIDBodyHeight+dQUPIDBaseThickness);
	G4cout<<dQUPIDWindowOuterHeight<<"\t"<<dQUPIDBodyHeight<<"\t"<<dQUPIDBaseThickness<<G4endl;
	G4double dQUPIDBottomWindowZOffset = BottomPhotoSensorsZOffset-dQUPIDWindowOuterHeight;//-dQUPIDWindowZCut;
	G4double dQUPIDBottomBodyZOffset = dQUPIDBottomWindowZOffset-0.5*dQUPIDBodyHeight+dQUPIDWindowZCut;
	G4double dQUPIDBottomBaseZOffset = dQUPIDBottomBodyZOffset-0.5*(dQUPIDBodyHeight+dQUPIDBaseThickness);
	G4cout<<dQUPIDWindowOuterHeight<<"\t"<<dQUPIDBodyHeight<<"\t"<<dQUPIDBaseThickness<<G4endl;
	G4cout<<dQUPIDTopWindowZOffset<<"\t"<<dQUPIDTopBodyZOffset<<"\t"<<dQUPIDTopBaseZOffset<<"\t"<<TopPhotoSensorsZOffset<<"\t"<<dPhotoSensorsHeight<<G4endl;
	G4cout<<dQUPIDBottomWindowZOffset<<"\t"<<dQUPIDBottomBodyZOffset<<"\t"<<dQUPIDBottomBaseZOffset<<"\t"<<BottomPhotoSensorsZOffset<<"\t"<<dPhotoSensorsHeight<<G4endl;

	//================================== top array ==================================
	G4double Rmax = dPSArrayOuterRadius;
	G4double Rsmall = dQUPIDBaseRadius;
	G4double rad_dist = dQUPIDsMinimumDisplacement;
	G4double Rmin = Rsmall + rad_dist;
	G4double Rbig = Rmax - Rsmall - rad_dist;
	G4int j=0;
	G4int i=0;
	G4int counter=0;
	//G4double Rnorm = 2*Rmax;
	G4double pi = 2*acos(0);
	G4double R;
	G4double theta; // angle between two tangent circunferences
	G4double l; // arc length
	G4double x, y;
	G4double remainder;
	if(Rmax<Rsmall)
	{
		G4cout<<"Rmax has to be bigger than Rsmall"<<G4endl; exit(1);
	}
	while((Rbig - 2*j*Rmin)>0)
	{
		R = Rbig - 2*j*Rmin;
		if(Rmin < R) 
		{
			theta = 2*asin(Rmin/R);
			l = R*theta;
			for(i=0;i<(int)(2*pi*R/l);i++)
			{
				remainder = (G4double) (((2*pi*R - ((G4int)(2*pi*R/l))*l)/R)/((G4int)(2*pi*R/l)));
				x = R*cos((theta+remainder)*i);
				y = R*sin((theta+remainder)*i);
				//G4cout<<counter<<"\t"<<x<<"\t"<<y<<G4endl;

				hVolumeName.str(""); hVolumeName << "QUPIDWindowNo" << counter+1;

				m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
				//m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDTopWindowZOffset), m_pQUPIDWindowLogicalVolume,
					hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
				hVolumeName.str(""); hVolumeName << "QUPIDBodyNo" << counter+1;
		
				m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
				//m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDTopBodyZOffset), m_pQUPIDBodyLogicalVolume,
					hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
				hVolumeName.str(""); hVolumeName << "QUPIDBaseNo" << counter+1;
		
				m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
				//m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDTopBaseZOffset), m_pQUPIDBaseLogicalVolume,
					hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
				counter++;
			}
		}
		else
		{
			hVolumeName.str(""); hVolumeName << "QUPIDWindowNo" << counter+1;

			m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
			//m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDTopWindowZOffset), m_pQUPIDWindowLogicalVolume,
				hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
			hVolumeName.str(""); hVolumeName << "QUPIDBodyNo" << counter+1;
		
			m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
			//m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDTopBodyZOffset), m_pQUPIDBodyLogicalVolume,
				hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
			hVolumeName.str(""); hVolumeName << "QUPIDBaseNo" << counter+1;
		
			m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(RotP180x,
			//m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDTopBaseZOffset), m_pQUPIDBaseLogicalVolume,
				hVolumeName.str(), m_pInnerGXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
				counter++;
		}
		j++;
	}
	G4cout<<"total number of QUPIDs per array: "<<counter<<G4endl;
	j=0;

	while((Rbig - 2*j*Rmin)>0)
	{
		R = Rbig - 2*j*Rmin;
		if(Rmin < R) 
		{
			theta = 2*asin(Rmin/R);
			l = R*theta;
			for(i=0;i<(int)(2*pi*R/l);i++)
			{
				remainder = (G4double) (((2*pi*R - ((G4int)(2*pi*R/l))*l)/R)/((G4int)(2*pi*R/l)));
				x = R*cos((theta+remainder)*i);
				y = R*sin((theta+remainder)*i);
				//G4cout<<counter<<"\t"<<x<<"\t"<<y<<G4endl;

				hVolumeName.str(""); hVolumeName << "QUPIDWindowNo" << counter+1;

				m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
				//m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDBottomWindowZOffset), m_pQUPIDWindowLogicalVolume,
					hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
				hVolumeName.str(""); hVolumeName << "QUPIDBodyNo" << counter+1;
		
				m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
				//m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDBottomBodyZOffset), m_pQUPIDBodyLogicalVolume,
					hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
				hVolumeName.str(""); hVolumeName << "QUPIDBaseNo" << counter+1;
		
				m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
				//m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
					G4ThreeVector(x,y,dQUPIDBottomBaseZOffset), m_pQUPIDBaseLogicalVolume,
					hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
					//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
				counter++;
			}
		}
		else
		{
			hVolumeName.str(""); hVolumeName << "QUPIDWindowNo" << counter+1;

			m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
			//m_hQUPIDWindowPhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDBottomWindowZOffset), m_pQUPIDWindowLogicalVolume,
				hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
			hVolumeName.str(""); hVolumeName << "QUPIDBodyNo" << counter+1;
		
			m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
			//m_hQUPIDBodyPhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDBottomBodyZOffset), m_pQUPIDBodyLogicalVolume,
				hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
		
			hVolumeName.str(""); hVolumeName << "QUPIDBaseNo" << counter+1;
		
			m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
			//m_hQUPIDBasePhysicalVolumes.push_back(new G4PVPlacement(0,
				G4ThreeVector(0.,0.,dQUPIDBottomBaseZOffset), m_pQUPIDBaseLogicalVolume,
				hVolumeName.str(), m_pOuterLXeLogicalVolume, false, counter));
				//hVolumeName.str(), m_pLabLogicalVolume, false, counter));
				counter++;
		}
		j++;
	}
	G4cout<<"total number of QUPIDs: "<<counter<<G4endl;

	//------------------------------- PMT sensitivity -------------------------------
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	pPmtSD = new DARWINPmtSensitiveDetector("DARWIN/PmtSD");
	pSDManager->AddNewDetector(pPmtSD);
	m_pQUPIDPhotocathodeLogicalVolume->SetSensitiveDetector(pPmtSD);

	//---------------------------------- attributes ---------------------------------

	m_pQUPIDPhotocathodeInteriorLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pQUPIDBodyInteriorLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hQUPIDWindowColor(1., 0.757, 0.024);
	G4VisAttributes *pQUPIDWindowVisAtt = new G4VisAttributes(hQUPIDWindowColor);
	pQUPIDWindowVisAtt->SetVisibility(true);
	//m_pQUPIDWindowLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pQUPIDBodyLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pQUPIDBaseLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pQUPIDWindowLogicalVolume->SetVisAttributes(pQUPIDWindowVisAtt);
	m_pQUPIDBodyLogicalVolume->SetVisAttributes(pQUPIDWindowVisAtt);
	m_pQUPIDBaseLogicalVolume->SetVisAttributes(pQUPIDWindowVisAtt);

	G4Colour hQUPIDPhotocathodeColor(1., 0.082, 0.011);
	G4VisAttributes *pQUPIDPhotocathodeVisAtt = new G4VisAttributes(hQUPIDPhotocathodeColor);
	pQUPIDPhotocathodeVisAtt->SetVisibility(true);
	m_pQUPIDPhotocathodeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pQUPIDPhotocathodeLogicalVolume->SetVisAttributes(pQUPIDPhotocathodeVisAtt);

	G4Colour hQUPIDAluminiumCoatingColor(1., 0.486, 0.027);
	G4VisAttributes *pQUPIDAluminiumCoatingVisAtt = new G4VisAttributes(hQUPIDAluminiumCoatingColor);
	pQUPIDAluminiumCoatingVisAtt->SetVisibility(true);
	//m_pQUPIDBodyAluminiumCoatingLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pQUPIDBaseAluminiumCoatingLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pQUPIDBodyAluminiumCoatingLogicalVolume->SetVisAttributes(pQUPIDAluminiumCoatingVisAtt);
	m_pQUPIDBaseAluminiumCoatingLogicalVolume->SetVisAttributes(pQUPIDAluminiumCoatingVisAtt);

	G4Colour hQUPIDAPDColor(0.788, 0.188, 0.780);
	G4VisAttributes *pQUPIDAPDVisAtt = new G4VisAttributes(hQUPIDAPDColor);
	pQUPIDAPDVisAtt->SetVisibility(true);
	//m_pQUPIDAPDLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pQUPIDAPDLogicalVolume->SetVisAttributes(pQUPIDAPDVisAtt);
}

void
DARWINDetectorConstruction::ConstructVetoPMTArrays()
{
	G4Material *Quartz = G4Material::GetMaterial("Quartz");
	G4Material *SS304LSteel = G4Material::GetMaterial("SS304LSteel");
	G4Material *Vacuum = G4Material::GetMaterial("Vacuum");
	G4Material *PhotoCathodeAluminium = G4Material::GetMaterial("PhotoCathodeAluminium");

	//==================================== 10" PMTs ===================================== (Hamamatsu R7081MOD-ASSY)
 	const G4double dPMTWindowOuterRadius = GetGeometryParameter("PMTWindowOuterRadius");
 	const G4double dPMTWindowOuterHalfZ = GetGeometryParameter("PMTWindowOuterHalfZ");
 	const G4double dPMTWindowTopZ = GetGeometryParameter("PMTWindowTopZ");
 	const G4double dPMTPhotocathodeOuterRadius = GetGeometryParameter("PMTPhotocathodeOuterRadius");
 	const G4double dPMTPhotocathodeOuterHalfZ = GetGeometryParameter("PMTPhotocathodeOuterHalfZ");
 	const G4double dPMTPhotocathodeTopZ = GetGeometryParameter("PMTPhotocathodeTopZ");
 	const G4double dPMTPhotocathodeInnerRadius = GetGeometryParameter("PMTPhotocathodeInnerRadius");
 	const G4double dPMTPhotocathodeInnerHalfZ = GetGeometryParameter("PMTPhotocathodeInnerHalfZ");
 	const G4double dPMTBodyOuterRadius = GetGeometryParameter("PMTBodyOuterRadius");
 	const G4double dPMTBodyInnerRadius = GetGeometryParameter("PMTBodyInnerRadius");
 	const G4double dPMTBodyHeight = GetGeometryParameter("PMTBodyHeight");
 	const G4double dPMTBaseOuterRadius = GetGeometryParameter("PMTBaseOuterRadius");
 	const G4double dPMTBaseInnerRadius = GetGeometryParameter("PMTBaseInnerRadius");
 	const G4double dPMTBaseHeight = GetGeometryParameter("PMTBaseHeight");
 	const G4double dPMTBaseInteriorHeight = GetGeometryParameter("PMTBaseInteriorHeight");

	//--------------------------------- PMT window ----------------------------------
	G4Ellipsoid *pPMTWindowEllipsoid = new G4Ellipsoid("PMTWindowEllipsoid", dPMTWindowOuterRadius, dPMTWindowOuterRadius, dPMTWindowOuterHalfZ, -dPMTWindowOuterHalfZ, dPMTWindowTopZ);

	m_pPMTWindowLogicalVolume = new G4LogicalVolume(pPMTWindowEllipsoid, Quartz, "PMTWindowVolume", 0, 0, 0);

	//--------------------------------- PMT photocathode ----------------------------------
	G4Ellipsoid *pPMTPhotocathodeEllipsoid = new G4Ellipsoid("PMTPhotocathodeEllipsoid", dPMTPhotocathodeOuterRadius, dPMTPhotocathodeOuterRadius, dPMTPhotocathodeOuterHalfZ, -dPMTPhotocathodeOuterHalfZ, dPMTPhotocathodeTopZ);

	m_pPMTPhotocathodeLogicalVolume = new G4LogicalVolume(pPMTPhotocathodeEllipsoid, PhotoCathodeAluminium, "PMTPhotocathodeVolume", 0, 0, 0);

	m_pPMTPhotocathodePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPMTPhotocathodeLogicalVolume, "PMTPhotocathode", m_pPMTWindowLogicalVolume, false, 0);

	//--------------------------------- PMT photocathode interior 1 ----------------------------------
	G4Ellipsoid *pPMTPhotocathodeInterior1Ellipsoid = new G4Ellipsoid("PMTPhotocathodeInterior1Ellipsoid", dPMTPhotocathodeInnerRadius, dPMTPhotocathodeInnerRadius, dPMTPhotocathodeInnerHalfZ, -dPMTPhotocathodeInnerHalfZ, dPMTPhotocathodeTopZ);

	m_pPMTPhotocathodeInterior1LogicalVolume = new G4LogicalVolume(pPMTPhotocathodeInterior1Ellipsoid, Vacuum, "PMTPhotocathodeInterior1LogicalVolume", 0, 0, 0);

	m_pPMTPhotocathodeInterior1PhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPMTPhotocathodeInterior1LogicalVolume, "PMTPhotocathodeInterior1", m_pPMTPhotocathodeLogicalVolume, false, 0);

	//--------------------------------- PMT photocathode interior 2 ----------------------------------
	G4Ellipsoid *pPMTPhotocathodeInterior2Ellipsoid = new G4Ellipsoid("PMTPhotocathodeInterior2Ellipsoid", dPMTPhotocathodeInnerRadius, dPMTPhotocathodeInnerRadius, dPMTPhotocathodeInnerHalfZ, dPMTPhotocathodeTopZ, dPMTWindowTopZ);

	m_pPMTPhotocathodeInterior2LogicalVolume = new G4LogicalVolume(pPMTPhotocathodeInterior2Ellipsoid, Vacuum, "PMTPhotocathodeInterior2LogicalVolume", 0, 0, 0);

	m_pPMTPhotocathodeInterior2PhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPMTPhotocathodeInterior2LogicalVolume, "PMTPhotocathodeInterior2", m_pPMTWindowLogicalVolume, false, 0);

	//--------------------------------- PMT body ----------------------------------
	const G4double dPMTBodyHalfZ = 0.5*dPMTBodyHeight;

	G4Tubs *pPMTBodyTubs = new G4Tubs("PMTBodyTubs", 0, dPMTBodyOuterRadius, dPMTBodyHalfZ, 0.*deg, 360.*deg);

	m_pPMTBodyLogicalVolume = new G4LogicalVolume(pPMTBodyTubs, SS304LSteel, "PMTBodyVolume", 0, 0, 0);

	//--------------------------------- PMT body interior ----------------------------------
	G4Tubs *pPMTBodyInteriorTubs = new G4Tubs("PMTBodyInteriorTubs", 0, dPMTBodyInnerRadius, dPMTBodyHalfZ, 0.*deg, 360.*deg);

	m_pPMTBodyInteriorLogicalVolume = new G4LogicalVolume(pPMTBodyInteriorTubs, Vacuum, "PMTBodyInteriorVolume", 0, 0, 0);

	m_pPMTBodyInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPMTBodyInteriorLogicalVolume, "PMTBodyInterior", m_pPMTBodyLogicalVolume, false, 0);

	//--------------------------------- PMT base (not really the base but the bottom part of the pmt)----------------------------------
	const G4double dPMTBaseHalfZ = 0.5*dPMTBaseHeight;

	G4Tubs *pPMTBaseTubs = new G4Tubs("PMTBaseTubs", 0, dPMTBaseOuterRadius, dPMTBaseHalfZ, 0.*deg, 360.*deg);

	m_pPMTBaseLogicalVolume = new G4LogicalVolume(pPMTBaseTubs, SS304LSteel, "PMTBaseVolume", 0, 0, 0);

	//--------------------------------- PMT base interior ----------------------------------
	const G4double dPMTBaseInteriorHalfZ = 0.5*dPMTBaseInteriorHeight;

	G4Tubs *pPMTBaseInteriorTubs = new G4Tubs("PMTBaseInteriorTubs", 0, dPMTBaseInnerRadius, dPMTBaseInteriorHalfZ, 0.*deg, 360.*deg);

	m_pPMTBaseInteriorLogicalVolume = new G4LogicalVolume(pPMTBaseInteriorTubs, Vacuum, "PMTBaseInteriorVolume", 0, 0, 0);

	m_pPMTBaseInteriorPhysicalVolume = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.),
		m_pPMTBaseInteriorLogicalVolume, "PMTBaseInterior", m_pPMTBaseLogicalVolume, false, 0);

	//================================== Liquid scintillator arrays ==================================
	G4int iNbTopPMTs = (G4int) GetGeometryParameter("NbTopPMTs");
	G4int iNbBottomPMTs = (G4int) GetGeometryParameter("NbBottomPMTs");
	G4int iNbLSPMTs = (G4int) GetGeometryParameter("NbLSPMTs");

	stringstream hVolumeName;

	for(G4int iPMTNb=iNbTopPMTs+iNbBottomPMTs; iPMTNb<iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs; iPMTNb++)
	{
	}

	//================================== Water arrays ==================================
	G4int iNbWaterPMTs = (G4int) GetGeometryParameter("NbWaterPMTs");

	for(G4int iPMTNb=iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs; iPMTNb<iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs+iNbWaterPMTs; iPMTNb++)
	{
		hVolumeName.str(""); hVolumeName << "WaterPMTWindowNo" << iPMTNb;

		m_hPMTWindowPhysicalVolumes.push_back(new G4PVPlacement(GetPMTRotation(iPMTNb),
			GetPMTPosition(iPMTNb, PMT_WINDOW), m_pPMTWindowLogicalVolume,
			hVolumeName.str(), m_pWaterLogicalVolume, false, iPMTNb));

		hVolumeName.str(""); hVolumeName << "WaterPMTBodyNo" << iPMTNb;

		m_hPMTBodyPhysicalVolumes.push_back(new G4PVPlacement(GetPMTRotation(iPMTNb),
			GetPMTPosition(iPMTNb, PMT_BODY), m_pPMTBodyLogicalVolume,
			hVolumeName.str(), m_pWaterLogicalVolume, false, iPMTNb));

		hVolumeName.str(""); hVolumeName << "WaterPMTBaseNo" << iPMTNb;

		m_hPMTBasePhysicalVolumes.push_back(new G4PVPlacement(GetPMTRotation(iPMTNb),
			GetPMTPosition(iPMTNb, PMT_BASE), m_pPMTBaseLogicalVolume,
			hVolumeName.str(), m_pWaterLogicalVolume, false, iPMTNb));
	}

	//------------------------------- PMT sensitivity -------------------------------
	m_pPMTPhotocathodeLogicalVolume->SetSensitiveDetector(pPmtSD);

	//---------------------------------- attributes ---------------------------------

	m_pPMTPhotocathodeInterior1LogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pPMTPhotocathodeInterior2LogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pPMTBodyInteriorLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	m_pPMTBaseInteriorLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hPMTWindowColor(1., 0.757, 0.024);
	G4VisAttributes *pPMTWindowVisAtt = new G4VisAttributes(hPMTWindowColor);
	pPMTWindowVisAtt->SetVisibility(false);
	m_pPMTWindowLogicalVolume->SetVisAttributes(pPMTWindowVisAtt);
	//m_pPMTWindowLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hPMTPhotocathodeColor(1., 0.082, 0.011);
	G4VisAttributes *pPMTPhotocathodeVisAtt = new G4VisAttributes(hPMTPhotocathodeColor);
	pPMTPhotocathodeVisAtt->SetVisibility(false);
	m_pPMTPhotocathodeLogicalVolume->SetVisAttributes(pPMTPhotocathodeVisAtt);
	//m_pPMTPhotocathodeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	G4Colour hPMTCasingColor(1., 0.486, 0.027);
	G4VisAttributes *pPMTCasingVisAtt = new G4VisAttributes(hPMTCasingColor);
	pPMTCasingVisAtt->SetVisibility(false);
	m_pPMTBodyLogicalVolume->SetVisAttributes(pPMTCasingVisAtt);
	m_pPMTBaseLogicalVolume->SetVisAttributes(pPMTCasingVisAtt);
	//m_pPMTBodyLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
	//m_pPMTBaseLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);
}

void
DARWINDetectorConstruction::ConstructSensitiveLXe()
{
	G4Material *LXe = G4Material::GetMaterial("LXe");
	G4Material *Teflon = G4Material::GetMaterial("Teflon");

	const G4double dSensitiveLXeHeight = GetGeometryParameter("SensitiveLXeHeight");
	const G4double dSensitiveLXeOuterRadius = GetGeometryParameter("SensitiveLXeOuterRadius");
	const G4double dCathodeToVeryBottomMesh = GetGeometryParameter("CathodeToVeryBottomMesh");
	const G4double dVeryBottomMeshToPhotoSensors = GetGeometryParameter("VeryBottomMeshToPhotoSensors");
	const G4double dPhotoSensorsHeight = GetGeometryParameter("PhotoSensorsHeight");
	const G4double dOuterLXeOuterHeight = GetGeometryParameter("OuterLXeOuterHeight");
	const G4double dOuterLXeHeight = GetGeometryParameter("OuterLXeHeight");

	G4double BottomPhotoSensorsZOffset = -0.5*dOuterLXeOuterHeight+dOuterLXeHeight+0.5*dPhotoSensorsHeight;
	// ============================ Sensitive LXe volume =============================
	G4Tubs *SensitiveLXeTubs = new G4Tubs("SensitiveLXeTubs",0.0*cm,dSensitiveLXeOuterRadius-0.1*cm,dSensitiveLXeHeight*0.5-0.1*cm,0.0*deg,360.0*deg);
	m_pSensitiveLXeLogicalVolume = new G4LogicalVolume(SensitiveLXeTubs,LXe,"m_pSensitiveLXeLogicalVolume");
	G4double SensitiveLXeXOffset = 0.0*cm;
	G4double SensitiveLXeYOffset = 0.0*cm;
	G4double SensitiveLXeZOffset = BottomPhotoSensorsZOffset + 0.5*dPhotoSensorsHeight + dCathodeToVeryBottomMesh + dVeryBottomMeshToPhotoSensors + 0.5*dSensitiveLXeHeight;

	G4Colour hLXeColor(0.094, 0.718, 0.812, 0.05);
	G4VisAttributes *pLXeVisAtt = new G4VisAttributes(hLXeColor);
	pLXeVisAtt->SetVisibility(true);
//	m_pLXeLogicalVolume->SetVisAttributes(pLXeVisAtt);
	m_pSensitiveLXeLogicalVolume->SetVisAttributes(G4VisAttributes::Invisible);

	m_pSensitiveLXePhysicalVolume = new G4PVPlacement(0, G4ThreeVector(SensitiveLXeXOffset, SensitiveLXeYOffset, SensitiveLXeZOffset),
		"m_pSensitiveLXePhysicalVolume", m_pSensitiveLXeLogicalVolume, m_pOuterLXePhysicalVolume, false, 0);

	//------------------------------ xenon sensitivity ------------------------------
	G4SDManager *pSDManager = G4SDManager::GetSDMpointer();

	DARWINLXeSensitiveDetector *pLXeSD = new DARWINLXeSensitiveDetector("DARWIN/LXeSD");
	pSDManager->AddNewDetector(pLXeSD);
	m_pSensitiveLXeLogicalVolume->SetSensitiveDetector(pLXeSD);
	//m_pGXeLogicalVolume->SetSensitiveDetector(pLXeSD);

	//================================== optical surfaces =================================	
	G4double dSigmaAlpha = 0.1;
	G4OpticalSurface *pTeflonOpticalSurface = new G4OpticalSurface("TeflonOpticalSurface",
		unified, groundbackpainted, dielectric_dielectric, dSigmaAlpha);
	pTeflonOpticalSurface->SetMaterialPropertiesTable(Teflon->GetMaterialPropertiesTable());	
	
	new G4LogicalBorderSurface("SidePTFELogicalBorderSurface",
		m_pSensitiveLXePhysicalVolume, m_pTPCPhysicalVolume, pTeflonOpticalSurface);

}

/*
void
DARWINDetectorConstruction::CheckOverlapping()
{
	G4cout << "Checking for overlapping... " << G4endl;
	int iOverlapped = 0;

	G4int res = 100000;
	G4bool verbose = true;

	iOverlapped += m_pLabPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pWaterTankPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pWaterPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pOuterCryostatPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pCryostatVacuumPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pInnerCryostatPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pOuterLXePhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pGXePhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pTPCPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pTopGridMeshPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pTopScreeningGridRingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pAnodeGridRingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pAnodeGridMeshPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pBelowLiquidGridRingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pBelowLiquidGridMeshPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pCathodeGridRingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pCathodeGridMeshPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pVeryBottomGridRingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pVeryBottomGridMeshPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pQUPIDPhotocathodePhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pQUPIDPhotocathodeInteriorPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pQUPIDBodyAluminiumCoatingPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pQUPIDBodyInteriorPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pQUPIDAPDPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pQUPIDBaseAluminiumCoatingPhysicalVolume->CheckOverlaps(res, 0, verbose);

	iOverlapped += m_pPMTPhotocathodePhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pPMTPhotocathodeInterior1PhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pPMTPhotocathodeInterior2PhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pPMTBodyInteriorPhysicalVolume->CheckOverlaps(res, 0, verbose);
	iOverlapped += m_pPMTBaseInteriorPhysicalVolume->CheckOverlaps(res, 0, verbose);

	const G4int iNbQUPIDs = 500;
	for(G4int iPMTNb=0; iPMTNb<iNbQUPIDs; iPMTNb++) iOverlapped += m_hQUPIDWindowPhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);
	for(G4int iPMTNb=0; iPMTNb<iNbQUPIDs; iPMTNb++) iOverlapped += m_hQUPIDBodyPhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);
	for(G4int iPMTNb=0; iPMTNb<iNbQUPIDs; iPMTNb++) iOverlapped += m_hQUPIDBasePhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);

	const G4int iNbPMTs = (G4int) GetGeometryParameter("NbPMTs");
	for(G4int iPMTNb=0; iPMTNb<iNbPMTs-iNbQUPIDs; iPMTNb++) iOverlapped += m_hPMTWindowPhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);
	for(G4int iPMTNb=0; iPMTNb<iNbPMTs-iNbQUPIDs; iPMTNb++) iOverlapped += m_hPMTBodyPhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);
	for(G4int iPMTNb=0; iPMTNb<iNbPMTs-iNbQUPIDs; iPMTNb++) iOverlapped += m_hPMTBasePhysicalVolumes[iPMTNb]->CheckOverlaps(res, 0, verbose);

	G4cout << "Done." << G4endl << G4endl;
	assert(iOverlapped==0);
}
*/

G4ThreeVector
DARWINDetectorConstruction::GetPMTPosition(G4int iPMTNb, PMTPart ePMTPart)
{
	const G4int iNbTopPMTs = (G4int) GetGeometryParameter("NbTopPMTs");
	const G4int iNbBottomPMTs = (G4int) GetGeometryParameter("NbBottomPMTs");
	const G4int iNbLSPMTs = (G4int) GetGeometryParameter("NbLSPMTs");
	const G4int iNbWaterPMTs = (G4int) GetGeometryParameter("NbWaterPMTs");

	G4ThreeVector hPos;

	if(iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs)
		hPos = GetPMTPositionLSArray(iPMTNb, ePMTPart);
	else if(iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs+iNbWaterPMTs)
		hPos = GetPMTPositionWaterArray(iPMTNb, ePMTPart);

	return hPos;
}


G4RotationMatrix *
DARWINDetectorConstruction::GetPMTRotation(G4int iPMTNb)
{
	const G4int iNbTopPMTs = (G4int) GetGeometryParameter("NbTopPMTs");
	const G4int iNbBottomPMTs = (G4int) GetGeometryParameter("NbBottomPMTs");

	const G4int iNbLSPMTs = (G4int) GetGeometryParameter("NbLSPMTs");
	const G4int iNbLSTopPMTs = (G4int) GetGeometryParameter("NbLSTopPMTs");
	const G4int iNbLSBottomPMTs = (G4int) GetGeometryParameter("NbLSBottomPMTs");
	const G4int iNbLSSidePMTColumns = (G4int) GetGeometryParameter("NbLSSidePMTColumns");

	const G4int iNbWaterPMTs = (G4int) GetGeometryParameter("NbWaterPMTs");
	const G4int iNbWaterTopPMTs = (G4int) GetGeometryParameter("NbWaterTopPMTs");
	const G4int iNbWaterBottomPMTs = (G4int) GetGeometryParameter("NbWaterBottomPMTs");
	const G4int iNbWaterSidePMTColumns = (G4int) GetGeometryParameter("NbWaterSidePMTColumns");

	G4RotationMatrix *pRotationMatrix = new G4RotationMatrix();

	if (iPMTNb < iNbTopPMTs)
		pRotationMatrix->rotateX(180.*deg);
	else if (iPMTNb < iNbTopPMTs+iNbBottomPMTs)
		pRotationMatrix->rotateX(0.*deg);
	else if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs)
	{
		if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSTopPMTs)
			pRotationMatrix->rotateX(0.*deg);
		else if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSTopPMTs+iNbLSBottomPMTs)
			pRotationMatrix->rotateX(180.*deg);
		else
		{
			pRotationMatrix->rotateY(-90.*deg);
			pRotationMatrix->rotateX(((iPMTNb-iNbTopPMTs-iNbBottomPMTs-iNbLSTopPMTs-iNbLSBottomPMTs)%iNbLSSidePMTColumns)*(360./iNbLSSidePMTColumns)*deg);
		}
	}
	else if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs+iNbWaterPMTs)
	{
		if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs+iNbWaterTopPMTs)
			pRotationMatrix->rotateX(0.*deg);
		else if (iPMTNb < iNbTopPMTs+iNbBottomPMTs+iNbLSPMTs+iNbWaterTopPMTs+iNbWaterBottomPMTs)
			pRotationMatrix->rotateX(180.*deg);
		else
		{
			pRotationMatrix->rotateY(-90.*deg);
			pRotationMatrix->rotateX(((iPMTNb-iNbTopPMTs-iNbBottomPMTs-iNbLSPMTs-iNbWaterTopPMTs-iNbWaterBottomPMTs)%iNbWaterSidePMTColumns)*(360./iNbWaterSidePMTColumns)*deg);
		}
	}

	return pRotationMatrix;
}

G4ThreeVector
DARWINDetectorConstruction::GetPMTPositionLSArray(G4int iPMTNb, PMTPart ePMTPart)
{
	const G4int iNbTopPMTs = (G4int) GetGeometryParameter("NbTopPMTs");
	const G4int iNbBottomPMTs = (G4int) GetGeometryParameter("NbBottomPMTs");
	const G4int iNbLSTopPMTs = (G4int) GetGeometryParameter("NbLSTopPMTs");
	const G4int iNbLSBottomPMTs = (G4int) GetGeometryParameter("NbLSBottomPMTs");
	const G4int iNbLSSidePMTColumns = (G4int) GetGeometryParameter("NbLSSidePMTColumns");
	const G4int iNbLSSidePMTRows = (G4int) GetGeometryParameter("NbLSSidePMTRows");

	const G4double dLSTopPMTWindowZ = GetGeometryParameter("LSTopPMTWindowZ");
	const G4double dLSBottomPMTWindowZ = GetGeometryParameter("LSBottomPMTWindowZ");
	const G4double dLSSidePMTWindowR = GetGeometryParameter("LSSidePMTWindowR");
	const G4double dPMTBodyOffset = GetGeometryParameter("PMTWindowTopZ") + 0.5*GetGeometryParameter("PMTBodyHeight");
	const G4double dPMTBaseOffset = GetGeometryParameter("PMTWindowTopZ") + GetGeometryParameter("PMTBodyHeight") + 0.5*GetGeometryParameter("PMTBaseHeight");

	const G4double dLSTopPMTDistance = GetGeometryParameter("LSTopPMTDistance");
	const G4double dLSBottomPMTDistance = GetGeometryParameter("LSBottomPMTDistance");
	const G4double dLSSidePMTRowDistance = GetGeometryParameter("LSSidePMTRowDistance");
//	const G4double dLSSidePMTColumnDistance = 2.*M_PI * dLSSidePMTWindowR / iNbLSSidePMTColumns;

	const G4int iNbLSTopPMTRows = 5;
	const G4int hLSTopPMTsPerRow[iNbLSTopPMTRows] = {1, 3, 5, 3, 1};
	const G4int iNbLSBottomPMTRows = 5;
	const G4int hLSBottomPMTsPerRow[iNbLSBottomPMTRows] = {1, 3, 5, 3, 1};

	G4ThreeVector hPos;

	iPMTNb -= iNbTopPMTs + iNbBottomPMTs;

	if (iPMTNb<iNbLSTopPMTs)
	{
		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setZ(dLSTopPMTWindowZ);
				break;

			case PMT_BODY:
				hPos.setZ(dLSTopPMTWindowZ + dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setZ(dLSTopPMTWindowZ + dPMTBaseOffset);
				break;
		}

		vector<G4double> hLSTopPMTsRowOffset;
		for(G4int iLSTopPMTRowNb=0; iLSTopPMTRowNb<iNbLSTopPMTRows; iLSTopPMTRowNb++)
			hLSTopPMTsRowOffset.push_back((G4double)(-0.5*(hLSTopPMTsPerRow[iLSTopPMTRowNb]-1)*dLSTopPMTDistance));

		G4int iPMTNb1 = 0;
		G4int iTotal = hLSTopPMTsPerRow[0];
		while(iPMTNb+1>iTotal)
			iTotal += hLSTopPMTsPerRow[++iPMTNb1];

		G4int iPMTNb2 = iPMTNb + hLSTopPMTsPerRow[iPMTNb1] - iTotal;

		hPos.setX(hLSTopPMTsRowOffset[iPMTNb1] + iPMTNb2*dLSTopPMTDistance);
		hPos.setY((0.5*(iNbLSTopPMTRows-1)-iPMTNb1)*dLSTopPMTDistance);
	}

	else if (iPMTNb<iNbLSTopPMTs+iNbLSBottomPMTs)
	{
		iPMTNb -= iNbLSTopPMTs;

		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setZ(dLSBottomPMTWindowZ);
				break;

			case PMT_BODY:
				hPos.setZ(dLSBottomPMTWindowZ - dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setZ(dLSBottomPMTWindowZ - dPMTBaseOffset);
				break;
		}

		vector<G4double> hLSBottomPMTsRowOffset;
		for(G4int iLSBottomPMTRowNb=0; iLSBottomPMTRowNb<iNbLSBottomPMTRows; iLSBottomPMTRowNb++)
			hLSBottomPMTsRowOffset.push_back((G4double)(-0.5*(hLSBottomPMTsPerRow[iLSBottomPMTRowNb]-1)*dLSBottomPMTDistance));

		G4int iPMTNb1 = 0;
		G4int iTotal = hLSBottomPMTsPerRow[0];
		while(iPMTNb+1>iTotal)
			iTotal += hLSBottomPMTsPerRow[++iPMTNb1];

		G4int iPMTNb2 = iPMTNb + hLSBottomPMTsPerRow[iPMTNb1] - iTotal;

		hPos.setX(hLSBottomPMTsRowOffset[iPMTNb1] + iPMTNb2*dLSBottomPMTDistance);
		hPos.setY((0.5*(iNbLSBottomPMTRows-1)-iPMTNb1)*dLSBottomPMTDistance);
	}

	else
	{
		iPMTNb -= iNbLSTopPMTs+iNbLSBottomPMTs;

		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setX(dLSSidePMTWindowR);
				break;

			case PMT_BODY:
				hPos.setX(dLSSidePMTWindowR + dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setX(dLSSidePMTWindowR + dPMTBaseOffset);
				break;
		}

		hPos.setY(0);
		hPos.setZ((0.5*(iNbLSSidePMTRows-1)-(iPMTNb/iNbLSSidePMTColumns))*dLSSidePMTRowDistance);
		hPos.rotateZ((iPMTNb%iNbLSSidePMTColumns)*(360./iNbLSSidePMTColumns)*deg);
	}

	return hPos;
}

G4ThreeVector
DARWINDetectorConstruction::GetPMTPositionWaterArray(G4int iPMTNb, PMTPart ePMTPart)
{
	const G4int iNbTopPMTs = (G4int) GetGeometryParameter("NbTopPMTs");
	const G4int iNbBottomPMTs = (G4int) GetGeometryParameter("NbBottomPMTs");
	const G4int iNbLSPMTs = (G4int) GetGeometryParameter("NbLSPMTs");
	const G4int iNbWaterTopPMTs = (G4int) GetGeometryParameter("NbWaterTopPMTs");
	const G4int iNbWaterBottomPMTs = (G4int) GetGeometryParameter("NbWaterBottomPMTs");
	const G4int iNbWaterSidePMTColumns = (G4int) GetGeometryParameter("NbWaterSidePMTColumns");
	const G4int iNbWaterSidePMTRows = (G4int) GetGeometryParameter("NbWaterSidePMTRows");

	const G4double dWaterTopPMTWindowZ = GetGeometryParameter("WaterTopPMTWindowZ");
	const G4double dWaterBottomPMTWindowZ = GetGeometryParameter("WaterBottomPMTWindowZ");
	const G4double dWaterSidePMTWindowR = GetGeometryParameter("WaterSidePMTWindowR");
	const G4double dPMTBodyOffset = GetGeometryParameter("PMTWindowTopZ") + 0.5*GetGeometryParameter("PMTBodyHeight");
	const G4double dPMTBaseOffset = GetGeometryParameter("PMTWindowTopZ") + GetGeometryParameter("PMTBodyHeight") + 0.5*GetGeometryParameter("PMTBaseHeight");

	const G4double dWaterTopPMTDistance = GetGeometryParameter("WaterTopPMTDistance");
	const G4double dWaterBottomPMTDistance = GetGeometryParameter("WaterBottomPMTDistance");
	const G4double dWaterSidePMTRowDistance = GetGeometryParameter("WaterSidePMTRowDistance");
//	const G4double dWaterSidePMTColumnDistance = 2.*M_PI * dWaterSidePMTWindowR / iNbWaterSidePMTColumns;

	const G4int iNbWaterTopPMTRows = 3;
	const G4int hWaterTopPMTsPerRow[iNbWaterTopPMTRows] = {3, 3, 3};
	const G4int iNbWaterBottomPMTRows = 5;
	const G4int hWaterBottomPMTsPerRow[iNbWaterBottomPMTRows] = {5, 5, 5, 5, 5};

	G4ThreeVector hPos;

	iPMTNb -= iNbTopPMTs + iNbBottomPMTs + iNbLSPMTs;

	if (iPMTNb<iNbWaterTopPMTs)
	{
		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setZ(dWaterTopPMTWindowZ);
				break;

			case PMT_BODY:
				hPos.setZ(dWaterTopPMTWindowZ + dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setZ(dWaterTopPMTWindowZ + dPMTBaseOffset);
				break;
		}

		vector<G4double> hWaterTopPMTsRowOffset;
		for(G4int iWaterTopPMTRowNb=0; iWaterTopPMTRowNb<iNbWaterTopPMTRows; iWaterTopPMTRowNb++)
			hWaterTopPMTsRowOffset.push_back((G4double)(-0.5*(hWaterTopPMTsPerRow[iWaterTopPMTRowNb]-1)*dWaterTopPMTDistance));

		G4int iPMTNb1 = 0;
		G4int iTotal = hWaterTopPMTsPerRow[0];
		while(iPMTNb+1>iTotal)
			iTotal += hWaterTopPMTsPerRow[++iPMTNb1];

		G4int iPMTNb2 = iPMTNb + hWaterTopPMTsPerRow[iPMTNb1] - iTotal;

		hPos.setX(hWaterTopPMTsRowOffset[iPMTNb1] + iPMTNb2*dWaterTopPMTDistance);
		hPos.setY((0.5*(iNbWaterTopPMTRows-1)-iPMTNb1)*dWaterTopPMTDistance);
	}

	else if (iPMTNb<iNbWaterTopPMTs+iNbWaterBottomPMTs) // WaterBottomPMTs
	{
		iPMTNb -= iNbWaterTopPMTs;

		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setZ(dWaterBottomPMTWindowZ);
				break;

			case PMT_BODY:
				hPos.setZ(dWaterBottomPMTWindowZ - dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setZ(dWaterBottomPMTWindowZ - dPMTBaseOffset);
				break;
		}

		vector<G4double> hWaterBottomPMTsRowOffset;
		for(G4int iWaterBottomPMTRowNb=0; iWaterBottomPMTRowNb<iNbWaterBottomPMTRows; iWaterBottomPMTRowNb++)
		  hWaterBottomPMTsRowOffset.push_back((G4double)(-0.5*(hWaterBottomPMTsPerRow[iWaterBottomPMTRowNb]-1)*dWaterBottomPMTDistance));			
		  //  G4cout << (G4double)(-0.5*(hWaterBottomPMTsPerRow[iWaterBottomPMTRowNb]-1)*dWaterBottomPMTDistance) << G4endl;}

		G4int iPMTNb1 = 0;
		G4int iTotal = hWaterBottomPMTsPerRow[0];
		while(iPMTNb+1>iTotal)
			iTotal += hWaterBottomPMTsPerRow[++iPMTNb1];

		G4int iPMTNb2 = iPMTNb + hWaterBottomPMTsPerRow[iPMTNb1] - iTotal;

		hPos.setX(hWaterBottomPMTsRowOffset[iPMTNb1] + iPMTNb2*dWaterBottomPMTDistance);
		hPos.setY((0.5*(iNbWaterBottomPMTRows-1)-iPMTNb1)*dWaterBottomPMTDistance);
		//G4cout << hWaterBottomPMTsRowOffset[iPMTNb1] + iPMTNb2*dWaterBottomPMTDistance << " " << (0.5*(iNbWaterBottomPMTRows-1)-iPMTNb1)*dWaterBottomPMTDistance << G4endl;
	}

	else // WaterLateralPMTs
	{
		iPMTNb -= iNbWaterTopPMTs+iNbWaterBottomPMTs;

		switch(ePMTPart)
		{
			case PMT_WINDOW:
				hPos.setX(dWaterSidePMTWindowR);
				break;

			case PMT_BODY:
				hPos.setX(dWaterSidePMTWindowR + dPMTBodyOffset);
				break;

			case PMT_BASE:
				hPos.setX(dWaterSidePMTWindowR + dPMTBaseOffset);
				break;
		}

		hPos.setY(0);
		hPos.setZ((0.5*(iNbWaterSidePMTRows-1)-(iPMTNb/iNbWaterSidePMTColumns))*dWaterSidePMTRowDistance);
		hPos.rotateZ((iPMTNb%iNbWaterSidePMTColumns)*(360./iNbWaterSidePMTColumns)*deg);
	}

	return hPos;
}

void DARWINDetectorConstruction::SetTeflonReflectivity(G4double dReflectivity)
{
	G4Material *pTeflonMaterial = G4Material::GetMaterial(G4String("Teflon"));

	if(pTeflonMaterial)
	{
		G4cout << "\n----> Setting Teflon reflectivity to " << dReflectivity << G4endl;

		G4MaterialPropertiesTable *pTeflonPropertiesTable = pTeflonMaterial->GetMaterialPropertiesTable();
		
		const G4int iNbEntries = 3;

		G4double teflon_PP[iNbEntries] = { 6.91 * eV, 6.98 * eV, 7.05 * eV };
		G4double teflon_REFL[iNbEntries] = {dReflectivity, dReflectivity, dReflectivity};
		pTeflonPropertiesTable->RemoveProperty("REFLECTIVITY");
		pTeflonPropertiesTable->AddProperty("REFLECTIVITY", teflon_PP, teflon_REFL, iNbEntries);
	}
	else
	{
		G4cout << "!!!!> Teflon material not found!" << G4endl;
		exit(-1);
	}
}


void DARWINDetectorConstruction::SetLXeScintillation(G4bool bScintillation)
{
	G4cout << "----> Setting LXe(GXe) scintillation to " << bScintillation << G4endl;
			
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
	if(pLXeMaterial)
	{	
	
		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
		if(bScintillation)
			pLXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1000./(1.0*keV));
	}
	else
	{
		G4cout << "ls!> LXe materials not found!" << G4endl;
		exit(-1);
	}
	
	G4Material *pGXeMaterial = G4Material::GetMaterial(G4String("GXe"));
	if(pGXeMaterial)
	{	
	
		G4MaterialPropertiesTable *pGXePropertiesTable = pGXeMaterial->GetMaterialPropertiesTable();
		if(bScintillation)
			pGXePropertiesTable->AddConstProperty("SCINTILLATIONYIELD", 1000./(1.0*keV));
	}
	else
	{
		G4cout << "ls!> GXe materials not found!" << G4endl;
		exit(-1);
	}


}



void DARWINDetectorConstruction::SetLXeAbsorbtionLength(G4double dAbsorbtionLength)
{
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));

	if(pLXeMaterial)
	{
		G4cout << "----> Setting LXe absorbtion length to " << dAbsorbtionLength/cm << " cm" << G4endl;

		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();
			
		const G4int iNbEntries = 3;

		G4double LXe_PP[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
		G4double LXe_ABSL[iNbEntries] = {dAbsorbtionLength, dAbsorbtionLength, dAbsorbtionLength};
		pLXePropertiesTable->RemoveProperty("ABSLENGTH");
		pLXePropertiesTable->AddProperty("ABSLENGTH", LXe_PP, LXe_ABSL, iNbEntries);
	}
	else
	{
		G4cout << "ls!> LXe materials not found!" << G4endl;
		exit(-1);
	}
}

void DARWINDetectorConstruction::SetLXeRayScatterLength(G4double dRayScatterLength)
{
	G4Material *pLXeMaterial = G4Material::GetMaterial(G4String("LXe"));
  
	if(pLXeMaterial)
	{

		G4cout << "----> Setting LXe scattering length to " << dRayScatterLength/cm << " cm" << G4endl;

		G4MaterialPropertiesTable *pLXePropertiesTable = pLXeMaterial->GetMaterialPropertiesTable();

		const G4int iNbEntries = 3;

		G4double LXe_PP[iNbEntries] = {6.91*eV, 6.98*eV, 7.05*eV};
		G4double LXe_SCAT[iNbEntries] = {dRayScatterLength, dRayScatterLength, dRayScatterLength};
		pLXePropertiesTable->RemoveProperty("RAYLEIGH");
		pLXePropertiesTable->AddProperty("RAYLEIGH", LXe_PP, LXe_SCAT, iNbEntries);
	}
	else
	{
		G4cout << "ls!> LXe materials not found!" << G4endl;
		exit(-1);
	}
}


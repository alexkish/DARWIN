#ifndef __DARWINDETECTORCONSTRUCTION_H__
#define __DARWINDETECTORCONSTRUCTION_H__

#include <globals.hh>

#include <vector>
#include <map>

using std::vector;
using std::map;

class G4Colour;
class G4LogicalVolume;
class G4VPhysicalVolume;
class DARWINDetectorMessenger;

#include "DARWINPmtSensitiveDetector.hh"

#include <G4VUserDetectorConstruction.hh>

class DARWINDetectorConstruction: public G4VUserDetectorConstruction
{
public:
	DARWINDetectorConstruction();
	~DARWINDetectorConstruction();

	G4VPhysicalVolume* Construct();

	void SetTeflonReflectivity(G4double dReflectivity);
	void SetLXeScintillation(G4bool dScintillation);
	void SetLXeAbsorbtionLength(G4double dAbsorbtionLength);
	void SetLXeRayScatterLength(G4double dRayScatterLength);

	static G4double GetGeometryParameter(const char *szParameter);


private:
	void DefineMaterials();
	void DefineGeometryParameters();

	void ConstructLaboratory();
	void ConstructVeto();
	void ConstructLSMVeto();
	void ConstructCryostat();
	void ConstructXenon();
	void ConstructTPC();
	void ConstructFieldCage();
	void ConstructQUPIDArrays();
	void ConstructVetoPMTArrays();
	void ConstructSensitiveLXe();

	void CheckOverlapping();
	void PrintGeometryInformation();

	typedef enum {PMT_WINDOW, PMT_BODY, PMT_BASE} PMTPart;

	G4ThreeVector GetPMTPosition(G4int iPMTNb, PMTPart ePMTPart);
	G4ThreeVector GetQUPIDPositionTopArray(G4int iPMTNb, PMTPart ePMTPart);
	G4ThreeVector GetQUPIDPositionBottomArray(G4int iPMTNb, PMTPart ePMTPart);
	G4ThreeVector GetPMTPositionLSArray(G4int iPMTNb, PMTPart ePMTPart);
	G4ThreeVector GetPMTPositionWaterArray(G4int iPMTNb, PMTPart ePMTPart);
	G4RotationMatrix *GetPMTRotation(G4int iPMTNb);

private:
	// rotation matrices
	G4RotationMatrix *m_pRotationXPlus225;
	G4RotationMatrix *m_pRotationXMinus225;
	G4RotationMatrix *m_pRotationXPlus45;
	G4RotationMatrix *m_pRotationXMinus45;
	G4RotationMatrix *m_pRotationXPlus90;
	G4RotationMatrix *m_pRotationXMinus90;
	G4RotationMatrix *m_pRotationX180;

	G4RotationMatrix *RotM90y;
	G4RotationMatrix *RotP90y;
	G4RotationMatrix *RotM90x;
	G4RotationMatrix *RotP270x;
	G4RotationMatrix *RotP90x;
	G4RotationMatrix *RotP90z;
	G4RotationMatrix *RotM90z;
	G4RotationMatrix *RotP0x;
	G4RotationMatrix *RotP180x;

	// logical volumes
	G4LogicalVolume *m_pMotherLogicalVolume;

	G4LogicalVolume *m_pLabLogicalVolume;

	G4LogicalVolume *m_pWaterTankLogicalVolume;
	G4LogicalVolume *m_pWaterLogicalVolume;
	G4LogicalVolume *m_pLSVesselLogicalVolume;
	G4LogicalVolume *m_pLSLogicalVolume;

	G4LogicalVolume *m_pOuterCryostatLogicalVolume;
	G4LogicalVolume *m_pCryostatVacuumLogicalVolume;
	G4LogicalVolume *m_pInnerCryostatLogicalVolume;
	G4LogicalVolume *m_pTopPipeLogicalVolume;

	G4LogicalVolume *m_pBellLogicalVolume;

	G4LogicalVolume *m_pOuterLXeLogicalVolume;
	G4LogicalVolume *m_pGXeLogicalVolume;
	G4LogicalVolume *m_pInnerGXeLogicalVolume;
	G4LogicalVolume *m_pSensitiveLXeLogicalVolume;

	G4LogicalVolume *m_pLXeLogicalVolume;

	G4LogicalVolume *m_pTPCLogicalVolume;

	G4LogicalVolume *m_pTopGridsRingLogicalVolume;
	G4LogicalVolume *m_pBottomGridsRingLogicalVolume;
	G4LogicalVolume *m_pGridMeshLogicalVolume;

	G4LogicalVolume *m_pQUPIDWindowLogicalVolume;
	G4LogicalVolume *m_pQUPIDPhotocathodeLogicalVolume;
	G4LogicalVolume *m_pQUPIDPhotocathodeInteriorLogicalVolume;
	G4LogicalVolume *m_pQUPIDBodyLogicalVolume;
	G4LogicalVolume *m_pQUPIDBodyAluminiumCoatingLogicalVolume;
	G4LogicalVolume *m_pQUPIDBodyInteriorLogicalVolume;
	G4LogicalVolume *m_pQUPIDAPDLogicalVolume;
	G4LogicalVolume *m_pQUPIDBaseLogicalVolume;
	G4LogicalVolume *m_pQUPIDBaseAluminiumCoatingLogicalVolume;
	G4LogicalVolume *m_pPhotoSensorsLogicalVolume;

	G4LogicalVolume *m_pPMTWindowLogicalVolume;
	G4LogicalVolume *m_pPMTPhotocathodeLogicalVolume;
	G4LogicalVolume *m_pPMTPhotocathodeInterior1LogicalVolume;
	G4LogicalVolume *m_pPMTPhotocathodeInterior2LogicalVolume;
	G4LogicalVolume *m_pPMTBodyLogicalVolume;
	G4LogicalVolume *m_pPMTBodyInteriorLogicalVolume;
	G4LogicalVolume *m_pPMTBaseLogicalVolume;
	G4LogicalVolume *m_pPMTBaseInteriorLogicalVolume;

	// physical volumes
	G4VPhysicalVolume *m_pLabPhysicalVolume;

	G4VPhysicalVolume *m_pWaterTankPhysicalVolume;
	G4VPhysicalVolume *m_pWaterPhysicalVolume;
	G4VPhysicalVolume *m_pLSVesselPhysicalVolume;
	G4VPhysicalVolume *m_pLSPhysicalVolume;

	G4VPhysicalVolume *m_pOuterCryostatPhysicalVolume;
	G4VPhysicalVolume *m_pCryostatVacuumPhysicalVolume;
	G4VPhysicalVolume *m_pInnerCryostatPhysicalVolume;
	G4VPhysicalVolume *m_pTopPipePhysicalVolume;

	G4VPhysicalVolume *m_pOuterLXePhysicalVolume;
	G4VPhysicalVolume *m_pGXePhysicalVolume;
	G4VPhysicalVolume *m_pInnerGXePhysicalVolume;

	G4VPhysicalVolume *m_pBellPhysicalVolume;

	G4VPhysicalVolume *m_pTPCPhysicalVolume;

	G4VPhysicalVolume *m_pTopGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pTopScreeningGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pAnodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pBelowLiquidGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pBelowLiquidGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pCathodeGridMeshPhysicalVolume;
	G4VPhysicalVolume *m_pVeryBottomGridRingPhysicalVolume;
	G4VPhysicalVolume *m_pVeryBottomGridMeshPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hQUPIDWindowPhysicalVolumes;
	G4VPhysicalVolume *m_pQUPIDPhotocathodePhysicalVolume;
	G4VPhysicalVolume *m_pQUPIDPhotocathodeInteriorPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hQUPIDBodyPhysicalVolumes;
	G4VPhysicalVolume *m_pQUPIDBodyAluminiumCoatingPhysicalVolume;
	G4VPhysicalVolume *m_pQUPIDBodyInteriorPhysicalVolume;
	G4VPhysicalVolume *m_pQUPIDAPDPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hQUPIDBasePhysicalVolumes;
	G4VPhysicalVolume *m_pQUPIDBaseAluminiumCoatingPhysicalVolume;
	G4VPhysicalVolume *m_pPhotoSensorsPhysicalVolume;

	vector<G4VPhysicalVolume *> m_hPMTWindowPhysicalVolumes;
	G4VPhysicalVolume *m_pPMTPhotocathodePhysicalVolume;
	G4VPhysicalVolume *m_pPMTPhotocathodeInterior1PhysicalVolume;
	G4VPhysicalVolume *m_pPMTPhotocathodeInterior2PhysicalVolume;
	vector<G4VPhysicalVolume *> m_hPMTBodyPhysicalVolumes;
	G4VPhysicalVolume *m_pPMTBodyInteriorPhysicalVolume;
	vector<G4VPhysicalVolume *> m_hPMTBasePhysicalVolumes;
	G4VPhysicalVolume *m_pPMTBaseInteriorPhysicalVolume;

	G4VPhysicalVolume *m_pSensitiveLXePhysicalVolume;

	static map<G4String, G4double> m_hGeometryParameters;
	
	DARWINDetectorMessenger *m_pDetectorMessenger;

	DARWINPmtSensitiveDetector *pPmtSD;
};

#endif // __XENON10PDETECTORCONSTRUCTION_H__


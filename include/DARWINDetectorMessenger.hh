#ifndef DARWINDetectorMessenger_h
#define DARWINDetectorMessenger_h 1

#include <G4UImessenger.hh>
#include <globals.hh>

class DARWINDetectorConstruction;

class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;
class G4UIcmdWithoutParameter;

class DARWINDetectorMessenger: public G4UImessenger
{
public:
	DARWINDetectorMessenger(DARWINDetectorConstruction *pXeDetector);
	~DARWINDetectorMessenger();

	void SetNewValue(G4UIcommand *pUIcommand, G4String hString);

private:
	DARWINDetectorConstruction* m_pXeDetector;

	G4UIdirectory *m_pDetectorDir;

	G4UIcmdWithADouble *m_pTeflonReflectivityCmd;
	G4UIcmdWithABool *m_pLXeScintillationCmd;
	G4UIcmdWithADoubleAndUnit *m_pLXeAbsorbtionLengthCmd;
	G4UIcmdWithADoubleAndUnit *m_pLXeRayScatterLengthCmd;

};

#endif


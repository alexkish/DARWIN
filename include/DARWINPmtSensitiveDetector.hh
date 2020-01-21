#ifndef __XENON10PPMTSENSITIVEDETECTOR_H__
#define __XENON10PPMTSENSITIVEDETECTOR_H__

#include <G4VSensitiveDetector.hh>

#include "DARWINPmtHit.hh"

class G4Step;
class G4HCofThisEvent;

class DARWINPmtSensitiveDetector: public G4VSensitiveDetector
{
public:
	DARWINPmtSensitiveDetector(G4String hName);
	~DARWINPmtSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	DARWINPmtHitsCollection* m_pPmtHitsCollection;
};

#endif // __XENON10PPMTSENSITIVEDETECTOR_H__


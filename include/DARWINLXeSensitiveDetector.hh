#ifndef __XENON10PLXESENSITIVEDETECTOR_H__
#define __XENON10PLXESENSITIVEDETECTOR_H__

#include <map>
#include <G4VSensitiveDetector.hh>

#include "DARWINLXeHit.hh"

using std::map;

class G4Step;
class G4HCofThisEvent;

class DARWINLXeSensitiveDetector: public G4VSensitiveDetector
{
public:
	DARWINLXeSensitiveDetector(G4String hName);
	~DARWINLXeSensitiveDetector();

	void Initialize(G4HCofThisEvent *pHitsCollectionOfThisEvent);
	G4bool ProcessHits(G4Step *pStep, G4TouchableHistory *pHistory);
	void EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent);

private:
	DARWINLXeHitsCollection* m_pLXeHitsCollection;

	map<int,G4String> m_hParticleTypes;
};

#endif // __XENON10PLXESENSITIVEDETECTOR_H__


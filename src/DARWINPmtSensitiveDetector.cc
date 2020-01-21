#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4VProcess.hh>
#include <G4ThreeVector.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>

#include <map>

using namespace std;

#include "DARWINPmtSensitiveDetector.hh"

DARWINPmtSensitiveDetector::DARWINPmtSensitiveDetector(G4String hName): G4VSensitiveDetector(hName)
{
	collectionName.insert("PmtHitsCollection");
}

DARWINPmtSensitiveDetector::~DARWINPmtSensitiveDetector()
{
}

void DARWINPmtSensitiveDetector::Initialize(G4HCofThisEvent* pHitsCollectionOfThisEvent)
{
	m_pPmtHitsCollection = new DARWINPmtHitsCollection(SensitiveDetectorName, collectionName[0]);

	static G4int iHitsCollectionID = -1;

	if(iHitsCollectionID < 0)
		iHitsCollectionID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	
	pHitsCollectionOfThisEvent->AddHitsCollection(iHitsCollectionID, m_pPmtHitsCollection); 
}

G4bool DARWINPmtSensitiveDetector::ProcessHits(G4Step* pStep, G4TouchableHistory *pHistory)
{
	G4Track *pTrack = pStep->GetTrack();

	if(pTrack->GetDefinition()->GetParticleName() == "opticalphoton")
	{
		DARWINPmtHit* pHit = new DARWINPmtHit();

		pHit->SetPosition(pStep->GetPreStepPoint()->GetPosition());
		pHit->SetTime(pTrack->GetGlobalTime());
		pHit->SetPmtNb(pTrack->GetTouchable()->GetVolume(1)->GetCopyNo());

		m_pPmtHitsCollection->insert(pHit);

		//		pHit->Print();
//		pHit->Draw();

		return true;
	}
	else
	{
		return false;
	}
}

void DARWINPmtSensitiveDetector::EndOfEvent(G4HCofThisEvent *pHitsCollectionOfThisEvent)
{
//  if (verboseLevel>0) { 
//     G4int NbHits = trackerCollection->entries();
//     G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
//            << " hits in the tracker chambers: " << G4endl;
//     for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
//    } 
}


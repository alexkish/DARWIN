#include <G4ios.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4Track.hh>
#include <G4Event.hh>
#include <G4VProcess.hh>
#include <G4StackManager.hh>

#include "DARWINAnalysisManager.hh"

#include "DARWINStackingAction.hh"

DARWINStackingAction::DARWINStackingAction(DARWINAnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

DARWINStackingAction::~DARWINStackingAction()
{
}

G4ClassificationOfNewTrack
DARWINStackingAction::ClassifyNewTrack(const G4Track *pTrack)
{
	G4ClassificationOfNewTrack hTrackClassification = fUrgent;

	if(pTrack->GetDefinition()->GetParticleType() == "nucleus" && !pTrack->GetDefinition()->GetPDGStable())
	{
		if(pTrack->GetParentID() > 0 && pTrack->GetCreatorProcess()->GetProcessName() == "RadioactiveDecay")
			hTrackClassification = fPostpone;
	}

	return hTrackClassification;
}

void
DARWINStackingAction::NewStage()
{
}

void
DARWINStackingAction::PrepareNewEvent()
{ 
}









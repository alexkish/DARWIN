#include <G4Event.hh>

#include "DARWINEventAction.hh"

DARWINEventAction::DARWINEventAction(DARWINAnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

DARWINEventAction::~DARWINEventAction()
{
}

void
DARWINEventAction::BeginOfEventAction(const G4Event *pEvent)
{
	if(pEvent->GetEventID() % 10000 == 0)
	{
		G4cout << G4endl;
		G4cout << "------ Begin event " << pEvent->GetEventID()
			<< "------" << G4endl;
	}
	
	if(m_pAnalysisManager)
		m_pAnalysisManager->BeginOfEvent(pEvent);
}

void DARWINEventAction::EndOfEventAction(const G4Event *pEvent)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->EndOfEvent(pEvent);
}



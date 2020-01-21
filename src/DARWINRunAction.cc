#include <Randomize.hh>

#include <sys/time.h>

#include "DARWINAnalysisManager.hh"

#include "DARWINRunAction.hh"

DARWINRunAction::DARWINRunAction(DARWINAnalysisManager *pAnalysisManager)
{
	m_pAnalysisManager = pAnalysisManager;
}

DARWINRunAction::~DARWINRunAction()
{

}

void
DARWINRunAction::BeginOfRunAction(const G4Run *pRun)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->BeginOfRun(pRun);

	struct timeval hTimeValue;
	gettimeofday(&hTimeValue, NULL);

	CLHEP::HepRandom::setTheEngine(new CLHEP::DRand48Engine);
	CLHEP::HepRandom::setTheSeed(hTimeValue.tv_usec);
}

void
DARWINRunAction::EndOfRunAction(const G4Run *pRun)
{
	if(m_pAnalysisManager)
		m_pAnalysisManager->EndOfRun(pRun);
}


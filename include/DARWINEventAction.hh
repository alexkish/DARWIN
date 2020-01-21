#ifndef __XENON10PEVENTACTION_H__
#define __XENON10PEVENTACTION_H__

#include <G4UserEventAction.hh>

#include "DARWINAnalysisManager.hh"

class G4Event;

class DARWINEventAction : public G4UserEventAction
{
public:
	DARWINEventAction(DARWINAnalysisManager *pAnalysisManager = 0);
	~DARWINEventAction();

public:
	void BeginOfEventAction(const G4Event *pEvent);
	void EndOfEventAction(const G4Event *pEvent);

private:
	DARWINAnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PEVENTACTION_H__


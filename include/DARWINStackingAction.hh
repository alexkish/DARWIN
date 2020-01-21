#ifndef __XENON10PSTACKINGACTION_H__
#define __XENON10PSTACKINGACTION_H__

#include <globals.hh>
#include <G4UserStackingAction.hh>

class DARWINAnalysisManager;

class DARWINStackingAction: public G4UserStackingAction
{
public:
	DARWINStackingAction(DARWINAnalysisManager *pAnalysisManager=0);
	~DARWINStackingAction();
  
	virtual G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* aTrack);
	virtual void NewStage();
	virtual void PrepareNewEvent();

private:
	DARWINAnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PSTACKINGACTION_H__


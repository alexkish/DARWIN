#ifndef __XENON10PRUNACTION_H__
#define __XENON10PRUNACTION_H__

#include <G4UserRunAction.hh>

class G4Run;

class DARWINAnalysisManager;

class DARWINRunAction: public G4UserRunAction
{
public:
	DARWINRunAction(DARWINAnalysisManager *pAnalysisManager=0);
	~DARWINRunAction();

public:
	void BeginOfRunAction(const G4Run *pRun);
	void EndOfRunAction(const G4Run *pRun);

private:
	DARWINAnalysisManager *m_pAnalysisManager;
};

#endif // __XENON10PRUNACTION_H__


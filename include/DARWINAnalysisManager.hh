#ifndef __DARWINANALYSISMANAGER_H__
#define __DARWINANALYSISMANAGER_H__

#include <globals.hh>

#include <TParameter.h>

class G4Run;
class G4Event;
class G4Step;

class TFile;
class TTree;

class DARWINEventData;
class DARWINPrimaryGeneratorAction;

class DARWINAnalysisManager
{
public:
	DARWINAnalysisManager(DARWINPrimaryGeneratorAction *pPrimaryGeneratorAction);
	virtual ~DARWINAnalysisManager();

public:
	virtual void BeginOfRun(const G4Run *pRun); 
	virtual void EndOfRun(const G4Run *pRun); 
	virtual void BeginOfEvent(const G4Event *pEvent); 
	virtual void EndOfEvent(const G4Event *pEvent); 
	virtual void Step(const G4Step *pStep);	

	void SetDataFilename(const G4String &hFilename) { m_hDataFilename = hFilename; }
	void SetNbEventsToSimulate(G4int iNbEventsToSimulate) { m_iNbEventsToSimulate = iNbEventsToSimulate; }

private:
	G4bool FilterEvent(DARWINEventData *pEventData);

private:
	G4int m_iLXeHitsCollectionID;
	G4int m_iPmtHitsCollectionID;

	G4String m_hDataFilename;
	G4int m_iNbEventsToSimulate;

	TFile *m_pTreeFile;
	TTree *m_pTree;
	TParameter<int> *m_pNbEventsToSimulateParameter;

	DARWINPrimaryGeneratorAction *m_pPrimaryGeneratorAction;

	DARWINEventData *m_pEventData;
};

#endif // __DARWINPANALYSISMANAGER_H__


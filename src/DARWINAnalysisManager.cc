#include <G4SDManager.hh>
#include <G4Run.hh>
#include <G4Event.hh>
#include <G4HCofThisEvent.hh>

#include <numeric>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TParameter.h>

#include "DARWINDetectorConstruction.hh"
#include "DARWINLXeHit.hh"
#include "DARWINPmtHit.hh"
#include "DARWINPrimaryGeneratorAction.hh"
#include "DARWINEventData.hh"

#include "DARWINAnalysisManager.hh"

DARWINAnalysisManager::DARWINAnalysisManager(DARWINPrimaryGeneratorAction *pPrimaryGeneratorAction)
{
	m_iLXeHitsCollectionID = -1;
	m_iPmtHitsCollectionID = -1;

	m_hDataFilename = "events.root";

	m_pPrimaryGeneratorAction = pPrimaryGeneratorAction;

	m_pEventData = new DARWINEventData();
}

DARWINAnalysisManager::~DARWINAnalysisManager()
{
}

void
DARWINAnalysisManager::BeginOfRun(const G4Run *pRun)
{
	m_pTreeFile = new TFile(m_hDataFilename.c_str(), "RECREATE", "File containing event data for DARWIN");
	m_pTree = new TTree("t1", "Tree containing event data for DARWIN");

	gROOT->ProcessLine("#include <vector>");

	m_pTree->Branch("eventid", &m_pEventData->m_iEventId, "eventid/I");
	m_pTree->Branch("ntpmthits", &m_pEventData->m_iNbTopPmtHits, "ntpmthits/I");
	m_pTree->Branch("nbpmthits", &m_pEventData->m_iNbBottomPmtHits, "nbpmthits/I");
	m_pTree->Branch("pmthits", "vector<int>", &m_pEventData->m_pPmtHits);
	m_pTree->Branch("etot", &m_pEventData->m_fTotalEnergyDeposited, "etot/F");
	m_pTree->Branch("nsteps", &m_pEventData->m_iNbSteps, "nsteps/I");
	
	m_pTree->Branch("trackid", "vector<int>", &m_pEventData->m_pTrackId);
	m_pTree->Branch("type", "vector<string>", &m_pEventData->m_pParticleType);
	m_pTree->Branch("parentid", "vector<int>", &m_pEventData->m_pParentId);
	m_pTree->Branch("parenttype", "vector<string>", &m_pEventData->m_pParentType);
	m_pTree->Branch("creaproc", "vector<string>", &m_pEventData->m_pCreatorProcess);
	m_pTree->Branch("edproc", "vector<string>", &m_pEventData->m_pDepositingProcess);
	m_pTree->Branch("xp", "vector<float>", &m_pEventData->m_pX);
	m_pTree->Branch("yp", "vector<float>", &m_pEventData->m_pY);
	m_pTree->Branch("zp", "vector<float>", &m_pEventData->m_pZ);
	m_pTree->Branch("ed", "vector<float>", &m_pEventData->m_pEnergyDeposited);
	m_pTree->Branch("time", "vector<float>", &m_pEventData->m_pTime);

	m_pTree->Branch("type_pri", "vector<string>", &m_pEventData->m_pPrimaryParticleType);
	m_pTree->Branch("xp_pri", &m_pEventData->m_fPrimaryX, 	"xp_pri/F");
	m_pTree->Branch("yp_pri", &m_pEventData->m_fPrimaryY, 	"yp_pri/F");
	m_pTree->Branch("zp_pri", &m_pEventData->m_fPrimaryZ, 	"zp_pri/F");
	m_pTree->Branch("e_pri",  &m_pEventData->m_fPrimaryE,	"e_pri/F");

	// write everything to one file, do not switch
	m_pTree->SetMaxTreeSize(10737418240LL); // 10G bytes, don't split file automatically

	m_pNbEventsToSimulateParameter = new TParameter<int>("nbevents", m_iNbEventsToSimulate);
	m_pNbEventsToSimulateParameter->Write();

}

void
DARWINAnalysisManager::EndOfRun(const G4Run *pRun)
{
	m_pTreeFile->Write();
	m_pTreeFile->Close();
}

void
DARWINAnalysisManager::BeginOfEvent(const G4Event *pEvent)
{
	if(m_iLXeHitsCollectionID == -1)
	{
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		m_iLXeHitsCollectionID = pSDManager->GetCollectionID("LXeHitsCollection");
	} 

	if(m_iPmtHitsCollectionID == -1)
	{
		G4SDManager *pSDManager = G4SDManager::GetSDMpointer();
		m_iPmtHitsCollectionID = pSDManager->GetCollectionID("PmtHitsCollection");
	}
}

void
DARWINAnalysisManager::EndOfEvent(const G4Event *pEvent)
{
	G4HCofThisEvent* pHCofThisEvent = pEvent->GetHCofThisEvent();
	DARWINLXeHitsCollection* pLXeHitsCollection = 0;
	DARWINPmtHitsCollection* pPmtHitsCollection = 0;

	G4int iNbLXeHits = 0, iNbPmtHits = 0;
	
	if(pHCofThisEvent)
	{
		if(m_iLXeHitsCollectionID != -1)
		{
			pLXeHitsCollection = (DARWINLXeHitsCollection *)(pHCofThisEvent->GetHC(m_iLXeHitsCollectionID));
			iNbLXeHits = (pLXeHitsCollection)?(pLXeHitsCollection->entries()):(0);
		}

		if(m_iPmtHitsCollectionID != -1)
		{
			pPmtHitsCollection = (DARWINPmtHitsCollection *)(pHCofThisEvent->GetHC(m_iPmtHitsCollectionID));
			iNbPmtHits = (pPmtHitsCollection)?(pPmtHitsCollection->entries()):(0);
		}
	}

	if(iNbLXeHits || iNbPmtHits)
	{
		m_pEventData->m_iEventId = pEvent->GetEventID();

		m_pEventData->m_pPrimaryParticleType->push_back(m_pPrimaryGeneratorAction->GetParticleTypeOfPrimary());

		m_pEventData->m_fPrimaryE = m_pPrimaryGeneratorAction->GetEnergyOfPrimary();

		m_pEventData->m_fPrimaryX = m_pPrimaryGeneratorAction->GetPositionOfPrimary().x();
		m_pEventData->m_fPrimaryY = m_pPrimaryGeneratorAction->GetPositionOfPrimary().y();
		m_pEventData->m_fPrimaryZ = m_pPrimaryGeneratorAction->GetPositionOfPrimary().z();

		G4int iNbSteps = 0;
		G4float fTotalEnergyDeposited = 0.;

		// LXe hits
		for(G4int i=0; i<iNbLXeHits; i++)
		{
			DARWINLXeHit *pHit = (*pLXeHitsCollection)[i];

			if(pHit->GetParticleType() != "opticalphoton")
			{
				m_pEventData->m_pTrackId->push_back(pHit->GetTrackId());
				m_pEventData->m_pParentId->push_back(pHit->GetParentId());

				m_pEventData->m_pParticleType->push_back(pHit->GetParticleType());
				m_pEventData->m_pParentType->push_back(pHit->GetParentType());
				m_pEventData->m_pCreatorProcess->push_back(pHit->GetCreatorProcess());
				m_pEventData->m_pDepositingProcess->push_back(pHit->GetDepositingProcess());

				m_pEventData->m_pX->push_back(pHit->GetPosition().x()/mm);
				m_pEventData->m_pY->push_back(pHit->GetPosition().y()/mm);
				m_pEventData->m_pZ->push_back(pHit->GetPosition().z()/mm);

				fTotalEnergyDeposited += pHit->GetEnergyDeposited()/keV;
				m_pEventData->m_pEnergyDeposited->push_back(pHit->GetEnergyDeposited()/keV);

				m_pEventData->m_pKineticEnergy->push_back(pHit->GetKineticEnergy()/keV);
				m_pEventData->m_pTime->push_back(pHit->GetTime()/second);

				iNbSteps++;
			}
		};

		m_pEventData->m_iNbSteps = iNbSteps;
		m_pEventData->m_fTotalEnergyDeposited = fTotalEnergyDeposited;

		//G4int iNbTopPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbTopPmts");
		//G4int iNbBottomPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbBottomPmts");
		//G4int iNbTopVetoPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbTopVetoPmts");
		//G4int iNbBottomVetoPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbBottomVetoPmts");
		//m_pEventData->m_pPmtHits->resize(iNbTopPmts+iNbBottomPmts+iNbTopVetoPmts+iNbBottomVetoPmts, 0);

		G4int iNbTopPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbTopPMTs");
		G4int iNbBottomPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbBottomPMTs");
		G4int iNbLSPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbLSPMTs");
		G4int iNbWaterPmts = (G4int) DARWINDetectorConstruction::GetGeometryParameter("NbWaterPMTs");
		m_pEventData->m_pPmtHits->resize(iNbTopPmts+iNbBottomPmts+iNbLSPmts+iNbWaterPmts, 0);

		// Pmt hits
		// Pmt hits
		for(G4int i=0; i<iNbPmtHits; i++)
			(*(m_pEventData->m_pPmtHits))[(*pPmtHitsCollection)[i]->GetPmtNb()]++;

		m_pEventData->m_iNbTopPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin(), m_pEventData->m_pPmtHits->begin()+iNbTopPmts, 0);
		m_pEventData->m_iNbBottomPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts, m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts, 0);
		m_pEventData->m_iNbLSPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts, m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts+iNbLSPmts, 0);
		m_pEventData->m_iNbWaterPmtHits =
			accumulate(m_pEventData->m_pPmtHits->begin()+iNbTopPmts+iNbBottomPmts+iNbLSPmts, m_pEventData->m_pPmtHits->end(), 0);

//      if((fTotalEnergyDeposited > 0. || iNbPmtHits > 0) && !FilterEvent(m_pEventData))
		//if(fTotalEnergyDeposited > 0. || iNbPmtHits > 0)
		if(fTotalEnergyDeposited > 0.)
			m_pTree->Fill();

		m_pEventData->Clear();
	}
}

void
DARWINAnalysisManager::Step(const G4Step *pStep)
{
}

/*
G4bool
DARWINAnalysisManager::FilterEvent(DARWINEventData *pEventData)
{
	G4double dEnergyDepositedSensitiveRegion = 0.;

	vector<float> *pX = pEventData->m_pX;
	vector<float> *pY = pEventData->m_pY;
	vector<float> *pZ = pEventData->m_pZ;
	vector<float> *pEnergyDeposited = pEventData->m_pEnergyDeposited;

	const G4double dDriftLength = DARWINDetectorConstruction::GetGeometryParameter("DriftLength");
	const G4double dRadius = DARWINDetectorConstruction::GetGeometryParameter("TeflonCylinderInnerRadius");

	for(G4int i=0; i<pEnergyDeposited->size(); i++)
	{
		if((*pZ)[i] < 0.-21.5 && (*pZ)[i] > -dDriftLength-21.5 && std::sqrt((*pX)[i]*(*pX)[i] + (*pY)[i]*(*pY)[i]) < dRadius)
			dEnergyDepositedSensitiveRegion += (*pEnergyDeposited)[i];
	}

//    if(dEnergyDepositedSensitiveRegion > 0. && dEnergyDepositedSensitiveRegion < 100.)
	if(dEnergyDepositedSensitiveRegion > 0.)
		return false;
	else
		return true;
}
*/
	

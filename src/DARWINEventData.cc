#include "DARWINEventData.hh"

DARWINEventData::DARWINEventData()
{
	m_iEventId = 0;
	m_iNbTopPmtHits = 0;
	m_iNbBottomPmtHits = 0;
	m_iNbLSPmtHits = 0;
	m_iNbWaterPmtHits = 0;
	m_pPmtHits = new vector<int>;

	m_fTotalEnergyDeposited = 0.;
	m_iNbSteps = 0;

	m_pTrackId = new vector<int>;
	m_pParentId = new vector<int>;
	m_pParticleType = new vector<string>;
	m_pParticlePdg = new vector<int>;
	m_pParentType = new vector<string>;
	m_pCreatorProcess = new vector<string>;
	m_pDepositingProcess = new vector<string>;
	m_pX = new vector<float>;
	m_pY = new vector<float>;
	m_pZ = new vector<float>;
	m_pEnergyDeposited = new vector<float>;
	m_pKineticEnergy = new vector<float>;
	m_pTime = new vector<float>;
	m_pNr = new vector<int>;

	m_pPrimaryParticleType = new vector<string>;
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
	m_fPrimaryCx = 0.;
	m_fPrimaryCy = 0.;
	m_fPrimaryCz = 0.;	
	m_fPrimaryE = 0.;	

	m_iNSave = 0;
	m_pSave_flag = new vector<int>;
	m_pSave_type = new vector<int>;
	m_pSave_x = new vector<float>;
	m_pSave_y = new vector<float>;
	m_pSave_z = new vector<float>;
	m_pSave_cx = new vector<float>;
	m_pSave_cy = new vector<float>;
	m_pSave_cz = new vector<float>;
	m_pSave_e = new vector<float>;

	m_iTotOptPhot = 0;
}

DARWINEventData::~DARWINEventData()
{
	delete m_pPmtHits;
	delete m_pTrackId;
	delete m_pParentId;
	delete m_pParticleType;
	delete m_pParticlePdg;
	delete m_pParentType;
	delete m_pCreatorProcess;
	delete m_pDepositingProcess;
	delete m_pX;
	delete m_pY;
	delete m_pZ;
	delete m_pEnergyDeposited;
	delete m_pKineticEnergy;
	delete m_pTime;
	delete m_pNr;

	delete m_pPrimaryParticleType;

	delete m_pSave_flag ;
	delete m_pSave_type ;
	delete m_pSave_x ;
	delete m_pSave_y ;
	delete m_pSave_z ;
	delete m_pSave_cx ;
	delete m_pSave_cy ;
	delete m_pSave_cz ;
	delete m_pSave_e ;
}

void
DARWINEventData::Clear()
{
	m_iEventId = 0;
	m_iNbTopPmtHits = 0;
	m_iNbBottomPmtHits = 0;
	m_iNbLSPmtHits = 0;
	m_iNbWaterPmtHits = 0;

	m_pPmtHits->clear();

	m_fTotalEnergyDeposited = 0.0;
	m_iNbSteps = 0;

	m_pTrackId->clear();
	m_pParentId->clear();
	m_pParticleType->clear();
	m_pParticlePdg->clear();
	m_pParentType->clear();
	m_pCreatorProcess->clear();
	m_pDepositingProcess->clear();
	m_pX->clear();
	m_pY->clear();
	m_pZ->clear();
	m_pEnergyDeposited->clear();
	m_pKineticEnergy->clear();
	m_pTime->clear();
	m_pNr->clear();

	m_pPrimaryParticleType->clear();
	m_fPrimaryX = 0.;
	m_fPrimaryY = 0.;
	m_fPrimaryZ = 0.;	
	m_fPrimaryCx = 0.;
	m_fPrimaryCy = 0.;
	m_fPrimaryCz = 0.;	
	m_fPrimaryE = 0.;	

	m_iNSave = 0;
	m_pSave_flag->clear();
	m_pSave_type->clear();
	m_pSave_x->clear(); 
	m_pSave_y->clear(); 
	m_pSave_z->clear(); 
	m_pSave_cx->clear(); 
	m_pSave_cy->clear(); 
	m_pSave_cz->clear(); 
	m_pSave_e->clear();

	m_iTotOptPhot = 0;
}


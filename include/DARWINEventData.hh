#ifndef __XENON10PEVENTDATA_H__
#define __XENON10PEVENTDATA_H__

#include <string>
#include <vector>

using std::string;
using std::vector;

class DARWINEventData
{
public:
	 DARWINEventData();
	~DARWINEventData();

public:
	void Clear();

public:
	int m_iEventId;								// the event ID
	int m_iNbTopPmtHits;						// number of top pmt hits
	int m_iNbBottomPmtHits;						// number of bottom pmt hits
	int m_iNbLSPmtHits;						// number of LS pmt hits
	int m_iNbWaterPmtHits;						// number of water pmt hits
	vector<int> *m_pPmtHits;					// number of photon hits per pmt
	float m_fTotalEnergyDeposited;				// total energy deposited in the ScintSD
	int m_iNbSteps;								// number of energy depositing steps
	vector<int> *m_pTrackId;					// id of the particle
	vector<int> *m_pParentId;					// id of the parent particle
	vector<string> *m_pParticleType;			// type of particle
	vector<int> *m_pParticlePdg;			// PDG code of particle
	vector<string> *m_pParentType;				// type of particle
	vector<string> *m_pCreatorProcess;			// interaction
	vector<string> *m_pDepositingProcess;		// energy depositing process
	vector<float> *m_pX;						// position of the step
	vector<float> *m_pY;
	vector<float> *m_pZ;
	vector<float> *m_pEnergyDeposited; 			// energy deposited in the step
	vector<float> *m_pKineticEnergy;			// particle kinetic energy after the step			
	vector<float> *m_pTime;						// time of the step
	vector<int> *m_pNr;						// NuclearRecoil (1) or EMrecoil (0)
	vector<string> *m_pPrimaryParticleType;		// type of particle
	float m_fPrimaryX;							// position of the primary particle
	float m_fPrimaryY;
	float m_fPrimaryZ;	
	float m_fPrimaryCx;							// direction of the primary particle
	float m_fPrimaryCy;
	float m_fPrimaryCz;	
        float m_fPrimaryE;							// Initial energy of the primary particle

  //MS In the following bank we save particle information in various positions, for ex. 
  //- entering the OuterCryostat from outside ... i.e. those crossing the whole shield
  //- optical photons when touching the tank surfaces
  //- active LXe volume (between anode and cathode)
  int	m_iNSave;              //
  vector<int> *m_pSave_flag;
  vector<int> *m_pSave_type;
  vector<float> *m_pSave_x;
  vector<float> *m_pSave_y;
  vector<float> *m_pSave_z;
  vector<float> *m_pSave_cx;
  vector<float> *m_pSave_cy;
  vector<float> *m_pSave_cz;
  vector<float> *m_pSave_e;

  int m_iTotOptPhot;
  float m_fTotPathWater;
};

#endif // __XENON10PEVENTDATA_H__


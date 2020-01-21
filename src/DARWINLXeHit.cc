#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "DARWINLXeHit.hh"

G4Allocator<DARWINLXeHit> DARWINLXeHitAllocator;

DARWINLXeHit::DARWINLXeHit() {}

DARWINLXeHit::~DARWINLXeHit()
{
	if(m_pParticleType) delete m_pParticleType;
	if(m_pParentType) delete m_pParentType;
	if(m_pCreatorProcess) delete m_pCreatorProcess;
	if(m_pDepositingProcess) delete m_pDepositingProcess;
}

DARWINLXeHit::DARWINLXeHit(const DARWINLXeHit &hDARWINLXeHit):G4VHit()
{
	m_iTrackId = hDARWINLXeHit.m_iTrackId;
	m_iParentId = hDARWINLXeHit.m_iParentId;
	m_pParticleType = hDARWINLXeHit.m_pParticleType;
	m_pParticlePdg = hDARWINLXeHit.m_pParticlePdg;
	m_pParentType = hDARWINLXeHit.m_pParentType ;
	m_pCreatorProcess = hDARWINLXeHit.m_pCreatorProcess ;
	m_pDepositingProcess = hDARWINLXeHit.m_pDepositingProcess ;
	m_hPosition = hDARWINLXeHit.m_hPosition;
	m_dEnergyDeposited = hDARWINLXeHit.m_dEnergyDeposited;
	m_dKineticEnergy = hDARWINLXeHit.m_dKineticEnergy ;
	m_dTime = hDARWINLXeHit.m_dTime;
}

const DARWINLXeHit &
DARWINLXeHit::operator=(const DARWINLXeHit &hDARWINLXeHit)
{
	m_iTrackId = hDARWINLXeHit.m_iTrackId;
	m_iParentId = hDARWINLXeHit.m_iParentId;
	m_pParticleType = hDARWINLXeHit.m_pParticleType;
	m_pParticlePdg = hDARWINLXeHit.m_pParticlePdg;
	m_pParentType = hDARWINLXeHit.m_pParentType ;
	m_pCreatorProcess = hDARWINLXeHit.m_pCreatorProcess ;
	m_pDepositingProcess = hDARWINLXeHit.m_pDepositingProcess ;
	m_hPosition = hDARWINLXeHit.m_hPosition;
	m_dEnergyDeposited = hDARWINLXeHit.m_dEnergyDeposited;
	m_dKineticEnergy = hDARWINLXeHit.m_dKineticEnergy ;
	m_dTime = hDARWINLXeHit.m_dTime;

	return *this;
}

G4int
DARWINLXeHit::operator==(const DARWINLXeHit &hDARWINLXeHit) const
{
	return ((this == &hDARWINLXeHit) ? (1) : (0));
}

void DARWINLXeHit::Draw()
{
	G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

	if(pVVisManager)
	{
		G4Circle hCircle(m_hPosition);
		G4Colour hColour(1.000, 0.973, 0.184);
		G4VisAttributes hVisAttributes(hColour);

		hCircle.SetScreenSize(0.1);
		hCircle.SetFillStyle(G4Circle::filled);
		hCircle.SetVisAttributes(hVisAttributes);
		pVVisManager->Draw(hCircle);
	}
}

void DARWINLXeHit::Print()
{
	G4cout << "-------------------- LXe hit --------------------" 
		<< "Id: " << m_iTrackId
		<< " Particle: " << *m_pParticleType
		<< " ParticlePdgCode: " << m_pParticlePdg
		<< " ParentId: " << m_iParentId
		<< " ParentType: " << *m_pParentType << G4endl
	       << "CreatorProcess: " << *m_pCreatorProcess << G4endl
		<< " DepositingProcess: " << *m_pDepositingProcess << G4endl
		<< "Position: " << m_hPosition.x()/mm
		<< " " << m_hPosition.y()/mm
		<< " " << m_hPosition.z()/mm
		<< " mm" << G4endl
		<< "EnergyDeposited: " << m_dEnergyDeposited/keV << " keV"
		<< " KineticEnergyLeft: " << m_dKineticEnergy/keV << " keV"
		<< " Time: " << m_dTime/s << " s" << G4endl;
}


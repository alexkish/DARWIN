#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

#include "DARWINPmtHit.hh"

G4Allocator<DARWINPmtHit> DARWINPmtHitAllocator;

DARWINPmtHit::DARWINPmtHit() {}

DARWINPmtHit::~DARWINPmtHit() {}

DARWINPmtHit::DARWINPmtHit(const DARWINPmtHit &hDARWINPmtHit):G4VHit()
{
	m_hPosition = hDARWINPmtHit.m_hPosition;
	m_dTime = hDARWINPmtHit.m_dTime;
	m_iPmtNb = hDARWINPmtHit.m_iPmtNb;
}

const DARWINPmtHit &
DARWINPmtHit::operator=(const DARWINPmtHit &hDARWINPmtHit)
{
	m_hPosition = hDARWINPmtHit.m_hPosition;
	m_dTime = hDARWINPmtHit.m_dTime;
	m_iPmtNb = hDARWINPmtHit.m_iPmtNb;
	
	return *this;
}

G4int
DARWINPmtHit::operator==(const DARWINPmtHit &hDARWINPmtHit) const
{
	return ((this == &hDARWINPmtHit) ? (1) : (0));
}

void DARWINPmtHit::Draw()
{
//    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
//    
//    if(pVVisManager)
//    {
//        G4Circle hCircle(m_hPosition);
//        G4Colour hColour(1.000, 0.973, 0.184);
//        G4VisAttributes hVisAttributes(hColour);
//        
//        hCircle.SetScreenSize(0.1);
//        hCircle.SetFillStyle(G4Circle::filled);
//        hCircle.SetVisAttributes(hVisAttributes);
//        pVVisManager->Draw(hCircle);
//    }
}

void DARWINPmtHit::Print()
{
	G4cout << "Pmt hit ---> " 
		<< "Pmt#" << m_iPmtNb
		<< " Position: " << m_hPosition.x()/mm
		<< " " << m_hPosition.y()/mm
		<< " " << m_hPosition.z()/mm
		<< " mm"
		<< " Time: " << m_dTime/s << " s" << G4endl;
}


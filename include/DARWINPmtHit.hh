#ifndef __XENON10PPMTHIT_H__
#define __XENON10PPMTHIT_H__

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>

class DARWINPmtHit: public G4VHit
{
public:
	DARWINPmtHit();
	~DARWINPmtHit();
	DARWINPmtHit(const DARWINPmtHit &);
	const DARWINPmtHit & operator=(const DARWINPmtHit &);
	G4int operator==(const DARWINPmtHit &) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	void Draw();
	void Print();

public:
	void SetPosition(G4ThreeVector hPosition) { m_hPosition = hPosition; }
	void SetTime(G4double dTime) { m_dTime = dTime; }
	void SetPmtNb(G4int iPmtNb) { m_iPmtNb = iPmtNb; }

	G4ThreeVector GetPosition() { return m_hPosition; }
	G4double GetTime() { return m_dTime; }
	G4int GetPmtNb() { return m_iPmtNb; }

private:
	G4ThreeVector m_hPosition;
	G4double m_dTime;
	G4int m_iPmtNb;
};

typedef G4THitsCollection<DARWINPmtHit> DARWINPmtHitsCollection;

extern G4Allocator<DARWINPmtHit> DARWINPmtHitAllocator;

inline void*
DARWINPmtHit::operator new(size_t)
{
	return((void *) DARWINPmtHitAllocator.MallocSingle());
}

inline void
DARWINPmtHit::operator delete(void *pDARWINPmtHit)
{
	DARWINPmtHitAllocator.FreeSingle((DARWINPmtHit*) pDARWINPmtHit);
}

#endif // __XENON10PPMTHIT_H__


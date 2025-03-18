#ifndef MJM_CONSTANTS_H__
#define MJM_CONSTANTS_H__
#include <math.h>
#include <cmath>


// wrong, check nist 
#define MJM_HART_EV 27.212
#define MJM_EV_CMINV 8065.544
#define MJM_ANG_EV 12398.4 
namespace mjm_constants
{


typedef double D;


/*
 should use nist etal
// http://www.ucl.ac.uk/~uccaati/Energy.html
// nist has diff value
const D hartree_ev=27.212;
const D ev_cminv= 8065.544;
//http://physics.nist.gov/Pubs/AtSpec/node01.html
const D ang_ev=12398.41;
2.417 989 40(21) × 1014 Hz
8 065.544 45(69) cm-1
1 239.841 91(11) nm 
*/
typedef const D (*CFUNC )(const D ) ;

const D hart2ev(const D x ) { return x*MJM_HART_EV; } 
const D hart2icm(const D x ) { return x*MJM_HART_EV*MJM_EV_CMINV; } 
const D hart2ang( const D x) { return MJM_ANG_EV/MJM_HART_EV/x; }   


class converter
{
D m_pow,m_factor;
public:
converter(const D f, const D e) : m_pow(e),m_factor(f) {}

D convert( const D& x )
{
return ::pow(x,m_pow)*m_factor;

}

};


//template<int e, D  f> D convert(const D x) { return ::pow(x,e)*f; } 



}; // mjm_constants




#endif


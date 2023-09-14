#ifndef MJM_TEXT_DATA_H__
#define MJM_TEXT_DATA_H__


#include <stdlib.h>
//#include <tcl.h>
//#include "mjm_tcl_base.h"
#include "mjm_solvers.h"
#include "mjm_jdftx_base.h"
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
// this should be moved and unified.
//#include <cblas.h>
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>

//#include "mjm_sample_spectrum.h"
//#include "mjm_io.h"
// the crap that needs this shold be moved.
#if NEED_COMP_STUFF
#include "prob_current.h"
#include "mjm_interchanges.h"
#include "mjm_dft_vibs.h"
#endif


// this really is specialized to tcl but not really

class mjm_text_data : public mjm_jdftx_typedefs
{
public:
typedef mjm_text_data Myt;


typedef std::vector<D> P;

// just for ease for now  
// argh......
//typedef std::vector<P> V;
//typedef prob_current ProbTy;
//typedef prob_current PC;
//typedef PC::Set_Selector WSelTy;
//typedef mjm_jdftx_ionpos IonsTy;
//typedef mjm_jdftx_forces ForcesTy;
//typedef mjm_2D_data SloMatTy;
//typedef mjm_2D_data_complex SloCoMatTy;
//typedef mjm_solvers_slow_mode SloMoTy;

/*
WSelTy m_sel;
PC m_p;
V m_data;
IdxTy m_n1, m_n2, m_n3;

IonsTy m_ions;
ForcesTy m_forces;

SloMatTy m_port;

// add a map of matricies for taking pert output and manipulating here.
// right now want to multiple a mass matrix term by term with force constants
// and generate frequencies. 

typedef std::map<StrTy,SloMatTy> MatMap;
MatMap m_mat_operands;


typedef pert_class Pert;
//typedef std::vector<Pert> Perts;
typedef perts_class Perts;
// where the f is this messe dup?
typedef simple_complex CoTy;
Perts m_perts;
// Hartrees/au^2/am -> mks finally s-1 to ev 
// amu-> grms/mole / 6e23 -> grams/atom * 1e-3 -> kg/atom
const D UC; // =6.023e23/1e-3/.529e-10/.529e-10*4.36e-18;
//const D hbar_ev=4.136e-15/6.28*1.03; // ??? kluge this is wrong
//const D hbar_ev=6.582119e-16*1.03; // ??? kluge this is wrong
const D HBAR_EV; //=6.582119e-16*1.0; // ??? kluge this is wrong
const D EV_INVERSE_CM; // 8065.5;
//ProbTy * m_p;
public:
mjm_jdftx_data():
//UC(6.023e23/1e-3/.529e-10/.529e-10*4.36e-18),
UC(6.022e23/1e-3/.529e-10/.529e-10*4.36e-18),
HBAR_EV(6.582119e-16*1.0), // ??? kluge this is wrong
 EV_INVERSE_CM(8065.5)
 
 {}
*/

~mjm_text_data() { }

#if 0 

IdxTy ion_count() { return 0;}
void load_ions(IsTy *is ) {}// { m_ions.load(is); 
//m_ions.compile(); 

//} 
void load_matrix(const StrTy & nm, IsTy *is ) { 
//SloMatTy m;
//m.load_lines(is);
//m_mat_operands[nm]=m;

}


void find_eigenvalues(const StrTy & nm , const StrTy & masses) 
{ 
/*`
SloMatTy  m=m_mat_operands[nm];
m.dump_matrix(std::cout,"MFCONSTS");
//m_port.dump_matrix(std::cout,"DYNMAT");
info()<<MM_MARK<<" now finding eigenvalues "<<CRLF;
if ( masses.length()!=0)
{
info()<<MM_MARK<<" now finding masses  "<<masses<<CRLF;
SloMatTy   mass=m_mat_operands[masses];
mass.dump_matrix(std::cout,"MASSES");
m.element_multiply(mass);
}
m.dump_matrix(std::cout,"MINPUTM");

SloMoTy smt;
const IdxTy modes=m.size();
// what the 
//smt.set_size(modes,modes);
info()<<MM_MARK<<" calling solve real  "<<modes<<CRLF;
info().flush();
mjm_solvers::solve_real(smt,m);
dump_omega_cm(smt, m.size() );
*/

}

template <class Os> void store_ions(Os &os, const StrTy & hdr ) 
{} // { m_ions.store(os,hdr); } 

void perturb_ion(const IdxTy atom, const IdxTy coord, const D delta)
{}// { m_ions.perturb(atom,coord,delta); } 
void perturb_ion_natural(const IdxTy atom, const IdxTy coord, const D delta)
{}
/*{
// this is stuypiud, should get both and return coords... 
//const IdxTy nearest=m_ions.nearest(atom,0);
//const IdxTy second_nearest=m_ions.nearest(atom,1);

m_ions.perturb_natural(atom,coord,delta);  
}
*/

void load_forces(IsTy *is ) {} // { m_forces.load(is); } 
template <class Os> void store_forces(Os &os, const StrTy & hdr ) 
{} // { m_forces.store(os,hdr); } 
void push_perturbed(IsTy * ions, IsTy * forces) 
{} // { m_perts.push_back(Pert(ions,forces)); }

/*
void make_phase_shifts(SloCoMatTy & retards, const IonsTy & ions,
const D omega, const D k, const IdxTy algo )
{} // {
*/

/*
const IdxTy modes=retards.size();
IdxTy ai=0;
IdxTy ci=0;
for ( IdxTy atomi=0; atomi<modes; ++atomi)
{
IdxTy aj=0;
IdxTy cj=0;
for ( IdxTy atomj=0; atomj<modes; ++atomj)
{
const D dx=m_ions.get(ai,0)-m_ions.get(aj,0);
const D dy=m_ions.get(ai,1)-m_ions.get(aj,1);
const D dz=m_ions.get(ai,2)-m_ions.get(aj,2);
const D dr=::sqrt(dx*dx+dy*dy+dz*dz);
const D c=3e8/.529e-10;
D dphi=dr/c*omega*k; // *1e5;
//if (ai<aj) dphi=-dphi;
info()<<MM_MARK<< "dhpi="<<dphi<<" omega="<<omega<< " c= " <<c<<CRLF;
retards(atomi,atomj)=CoTy(::cos(dphi),::sin(dphi));
++cj; if ( cj==3) { cj=0; ++aj; } 
}
++ci; if ( ci==3) { ci=0; ++ai; } 
} // atomi

}
*/

#endif


#if 0 
// this code was hacked terribly due to use of c instead oz 
void include_retardation()
{
//if ( m_port.size()==0) m_perts.find_perturbations(m_port,m_ions,m_forces);

//m_port.dump_sparse(std::cout);
info()<<MM_MARK<<" now finding non-retarded eigenvalues "<<CRLF;
SloCoMatTy static_omega,latest_omega;
// these are the eigenvalues returned directly
find_perturbed(&static_omega);
const IdxTy modes=static_omega.size();
info()<<MM_MARK<<" now finding non-retarded eigenvalues modes= "<<modes<<CRLF;
const IdxTy nfac=3;
D ff[nfac];
ff[0]=1e4;
ff[1]=2e4;
ff[2]=5e4;
for ( IdxTy ffi=0; ffi<nfac; ++ffi)
{
for ( IdxTy ei=0; ei<modes; ++ei)
{
info()<<MM_MARK<<" now finding non-retarded eigenvalues mode ei= "<<ei<<CRLF;
CoTy mode=static_omega(ei,0);
//D omega=::sqrt(::fabs(mode.r()));
D omega2=mode.r();
const D ev=::sqrt(::fabs(omega2*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
const D cminverse=ev*EV_INVERSE_CM;
const D omega=cminverse*3e10;
SloCoMatTy retards;
retards.set_size(modes,modes);

make_phase_shifts(retards,m_ions,omega,ff[ffi],0);
SloCoMatTy fc=m_port;
fc.element_multiply(retards);
SloMoTy smt;
mjm_solvers::solve_complex(smt,fc);
smt.export_eigenvalues(latest_omega);
//solve_complex(latest_omega,fc);
info()<<MM_MARK<<"static cf="<<ff[ffi]<<" o="<< cminverse<<" :  ";
IdxTy einew=ei;
for (einew=0; einew<modes; ++einew)
{
const D omeganew=latest_omega(einew,0).r();
//const D omeganewi=latest_omega(einew,0).i();
const D evc=::sqrt(::fabs(omeganew*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
const D cminversec=evc*EV_INVERSE_CM;
//const D omegac=cminversec*3e10;
info()<<cminversec<<" ";
//info()<<MM_MARK<<" now  static "<< cminverse<<" new "<<cminversec<<" "<<omeganew<<" "<<omeganewi<<CRLF;
}
info()<<CRLF;
//static_omega.dump_matrix(info());
} // ei
} //ff
info()<<MM_MARK<<" now DUN finding non-retarded eigenvalues "<<CRLF;
}


SloMoTy find_perturbed(SloCoMatTy * eireal=0, const StrTy param="")


{
// 

info()<<MM_MARK<<" ions compile should be done alread  COMPILE reference matrix  "<<CRLF;
//m_ions.compile();

m_perts.find_perturbations(m_port,m_ions,m_forces,param);
m_port.dump_matrix(std::cout,"FCONST");
// now there is no mass yet... 
SloMatTy masses;
m_perts.make_masses(masses,m_ions);

info()<<MM_MARK<<" now trying to clean up matrix  "<<CRLF;
//m_perts.clean_matrix(m_port,m_ions);
//m_perts.kluge_matrix(m_port,m_ions);

info()<<MM_MARK<<" now adjusting masses  "<<CRLF;
m_port.element_multiply(masses);

m_port.dump_matrix(std::cout,"DYNMAT");
info()<<MM_MARK<<" now finding eigenvalues "<<CRLF;
SloMoTy smt;
const IdxTy modes=m_port.size();
// what the 
//smt.set_size(modes,modes);
info()<<MM_MARK<<" calling solve real  "<<modes<<CRLF;
info().flush();
mjm_solvers::solve_real(smt,m_port);
if(eireal!=0)
{
info()<<MM_MARK<<" exporting  eigenvalues "<<CRLF;
info().flush();
 smt.export_eigenvalues(*eireal);
}

info()<<MM_MARK<<" dun finding eigenvalues "<<CRLF;
info().flush();
const IdxTy xN=m_port.size();
double xWR[xN],xWI[xN];
smt.export_eigenvalues(xWR,xWI);
for (IdxTy  i=0; i<xN; ++i)
{
const D ev=::sqrt(::fabs(xWR[i]*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
const D cminverse=ev*EV_INVERSE_CM;
 info()<<MM_MARK<<" "<<i<<" "<<ev<<" ( "<<cminverse<<" ) "<<xWR[i]<<" "<<xWI[i]<<CRLF;
}
//smt.dump_modes(info(),m_ions);
return smt;

}

void dump_omega_cm(  SloMoTy & smt, const IdxTy & n ) const
{
const IdxTy xN=n; //m_port.size();
double xWR[xN],xWI[xN];
smt.export_eigenvalues(xWR,xWI); // should be a const function?
for (IdxTy  i=0; i<xN; ++i)
{
const D ev=::sqrt(::fabs(xWR[i]*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
const D cminverse=ev*EV_INVERSE_CM;
 info()<<MM_MARK<<" "<<i<<" "<<ev<<" ( "<<cminverse<<" ) "<<xWR[i]<<" "<<xWI[i]<<CRLF;
}
}



void analyze_modes(const SloMatTy & m_port,const IonsTy & m_ions,const SloMatTy  &vrd)
{
const IdxTy sz=vrd.size();
typedef std::vector<D> ModalTy;
ModalTy m;
// for each mode, ZZ
for (IdxTy mode=0; mode<sz; ++mode)
{
m.clear();
for (IdxTy q=0; q<sz; ++q) m.push_back(vrd(mode,q));

// this should be sz/3
const IdxTy atoms=m_ions.size();
for (IdxTy atom=0; atom<atoms; ++atom)
{



} // atom



} // mode

}
// when loading, only take from s and count of n
// if wfn is odd, ajust energy range to match
void select_range(const IdxTy s, const IdxTy n)
{
m_sel.select_range(s,n);

}
// tis contains the grid dims that we never see
void load_selected_wfn(IsTy * is)
{
std::cout<<MM_MARK<<" loading selected wfn"<<CRLF;
m_p.add(is,m_sel);
}
// when using separate files
void set_grid_dims( const IdxTy n1, const IdxTy n2, const IdxTy n3)
{ m_n1=n1; m_n2=n2; m_n3=n3;
info()<<MM_MARK<<" setting gird dim "<<n1<<" "<<n2<<" "<<n3<<CRLF;
 }
bool load_selected_wfn_raw(const StrTy nm, const IdxTy i)
{
info()<<MM_MARK<<" loading selected wfn"<<nm<<CRLF;
if ( m_n1*m_n2*m_n3==0) 
	info()<<MM_MARK<<" should set dims first now zero"<<CRLF;
return m_p.add_raw_wfn(nm,i,m_sel,m_n1,m_n2,m_n3);
}
// these layouts must match each other or inconsistent results occur
void load_occupancies(IsTy * is)
{
std::cout<<MM_MARK<<" loading occipancies  "<<" "<<CRLF;
m_p.load_selected_weights(is,m_sel);

}

void load_selected_energies(IsTy * is)
{
m_p.add_energies(is,0,m_sel);
}
void find_dipoles()
{
m_p.setup_noithing(m_sel);
std::cout<<MM_MARK<<" doing dipole  "<<" "<<CRLF;
m_p.dipole_etc();
}



operator const V () const { return m_data; } 
void add(const D x, const D y)
{
P p;
p.push_back(x);
p.push_back(y);
m_data.push_back(p);
}
// arrgghhhh   
D & x( const IdxTy i) { return m_data[i][0]; } 
D & y( const IdxTy i) { return m_data[i][1]; } 
IdxTy size() const { return m_data.size(); } 

#endif
public:
void load(IsTy * is){} 
void store(OsTy & os) {}
#if 0 
{
IdxTy i=0;
D x,y;
while ( !is->eof()&&is->good())
{
(*is)>>x>>y;
bool add_this_one=true;
if ( !(!is->eof()&&is->good())) break;
if ( add_this_one)
{
//std::cout<<MM_MARK<<" adding energfy "<<i<<" energ "<<e<<" "<<i<<CRLF;
add(x,y);
++i;
} // add
}

}
void store(OsTy & os)
{
const IdxTy sz = m_data.size();
std::cout<<MM_MARK<<" storing size is "<<sz<<CRLF;
IdxTy i=0;
while ( i<sz)
{
os<<m_data[i][0]<<" "<<m_data[i][1]<<CRLF;
++i;
}


}
/*
// output name, kernel, start, delta, end, kernel param
template < class Kernel > void regular(Myt & dest, const Kernel & k, const D start, const D delta, const D end)
{
std::cout<<MM_MARK<<" daa size is "<<m_data.size()<<" start="<<start<<
" end "<<end<< " delta="<<delta<<CRLF;
for ( D t=start; t<end; t+=delta)
{
D ygrid= k.eval(*this,t);
dest.add(t,ygrid);

} // t loop 



}
*/

#endif // zero 



};
/*
template <class SigClass> class mjm_tcl_base


*/
#endif


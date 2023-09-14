#ifndef MJM_SOLVERS_BASE_H__
#define MJM_SOLVERS_BASE_H__
 

#include <stdlib.h>
//#include <tcl.h>
//#include "mjm_tcl_base.h"
//#include "mjm_jdftx_base.h"
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
// need a non-template signatujre wtf
// #include <cblas.h>
// not sure this works with linpack?
#include <complex>
#include <map>
#include <vector>

//#include "mjm_sample_spectrum.h"
#include "mjm_io.h"
// the crap that needs this shold be moved.
//#include "prob_current.h"
#include "mjm_interchanges.h"

#if 0 
// these really need to go into a c++ file and be hidden from all esle
extern "C"
{
// double eignen value
void dgeev_(char * JOBVL, char * JOBVR, int * N, double * A,
int * LDA, double * WR, double * WL, double * VL,int * LDVL,
double * VR, int * LDVR, double * WORK, int * LWORK, int * INFO); 
// complex eignen values tjhis does not work on most data wtfv
// complex is fing float not doujlb
// DO NOT USE C, USE Z complex uses FLOAT not DOUBLE WTF
void cgeev_(char * JOBVL, char * JOBVR, int * N, double * A,
int * LDA, double * W, double * VL,int * LDVL,
double * VR, int * LDVR, double * WORK, int * LWORK, double * RWORK, int * INFO); 
// us zx not fing c
void zgeev_(char * JOBVL, char * JOBVR, int * N, double * A,
int * LDA, double * W, double * VL,int * LDVL,
double * VR, int * LDVR, double * WORK, int * LWORK, double * RWORK, int * INFO); 


// DGEBRD real double SVD
void dgesvd_(char * JOBU, char * JOBVT,int * M, int * N, double * A,
 int *LDA, double * S, double * U, int * LDU, double * VT, int * LDVT, double * WORK, int * LWORK,
int * INFO);



// more robust??? No, do not USE c use Z 
void cgees_(char * JOBVS, char * SORT, int * SELECT,int * N,  double * A,
int * LDA, int  * SDIM, double * W, double * VS,int * LDVS,
 double * WORK, int * LWORK, double * RWORK, double * BWORK, int * INFO); 

void zgees_(char * JOBVS, char * SORT, int * SELECT,int * N,  double * A,
int * LDA, int  * SDIM, double * W, double * VS,int * LDVS,
 double * WORK, int * LWORK, double * RWORK, double * BWORK, int * INFO); 

void dgees_(char * JOBVS, char * SORT, int * SELECT,int * N,  double * A,
int * LDA, int  * SDIM, double * WR, double * WI, double * VS,int * LDVS,
 double * WORK, int * LWORK, double * BWORK, int * INFO); 


};

#endif


class mjm_solvers_typedefs
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef double D;
typedef simple_complex CoTy;

static OsTy & info() { return std::cout; }


};

// these use the stl junk ratjer than arrays with
// are own memory mamangers
// we should have our own pool for transient objects
class mjm_solvers_slow_mode : public mjm_solvers_typedefs
{
typedef mjm_solvers_slow_mode Myt;
typedef std::vector<CoTy> Eigenvalues;
typedef std::vector<CoTy> ModeShape;
typedef std::vector<ModeShape> ModeShapes;

Eigenvalues m_values;
ModeShapes m_shapes;

public:
// generally eigenvalues will exist but maybe not modes
IdxTy size() const { return m_values.size(); }
CoTy eigenvalue(const IdxTy i) { return m_values[i]; }
ModeShape &  eigenvector(const IdxTy i) { return m_shapes[i]; }
const ModeShape &  eigenvector(const IdxTy i) const { return m_shapes[i]; }

void add_eigenvalue(const CoTy & e) { m_values.push_back(e); }
void add_eigenvalue(const D & r, const D &i ) { m_values.push_back(CoTy(r,i)); } 
void add_eigenvalues(const D * const  r, const D * const i, const IdxTy n ) 
{ 
for ( IdxTy j=0; j<n; ++j) { add_eigenvalue(r[j],i[j]); }
} 
// complex in r,i format as douvbles not floats...
void add_eigenvalues(const D * const  w, const IdxTy n ) 
{
IdxTy p=0; 
for ( IdxTy j=0; j<n; ++j) { add_eigenvalue(w[p],w[p+1]);p+=2; }
} 
// add using the convention these are real unless
// omega is complex then theyh have real followed by imaginart
void add_eigenvector_columns_rc(const D * const p, const IdxTy n)
{

const D *  ptr=p;

for ( IdxTy i=0; i<n; ++i) 
{ 
bool is_real=(eigenvalue(i).i()==0);
// pick up real parts 
add_mode(ptr,n,n,1); 
if (!is_real)
{
// duplicate the real part
add_mode(ptr,n,n,1); 
// now get im and cc
mode(m_shapes[i],ptr,1,n,1); 
mode(m_shapes[i+1],ptr,-1,n,1); 
++i;
}
ptr+=1; 

}

}
void add_eigenvector_columns(const D * const p, const IdxTy n)
{
const D *  ptr=p;
for ( IdxTy i=0; i<n; ++i) { add_mode(ptr,n,n); ptr+=2;  }

}
// point to an array of complex mode things
void add_mode(const D * const p, const IdxTy n,const IdxTy pitch, const IdxTy pf=2)
{
D const *  ptr=p;
ModeShape ms;
for ( IdxTy i=0; i<n; ++i)
{
if (pf!=1)ms.push_back(CoTy(*ptr,*(ptr+1)));
else ms.push_back(CoTy(*ptr,0));
// pitch is in numbers, complex is 2x
ptr+=pitch*pf;

} 
m_shapes.push_back(ms);
}

void mode(ModeShape &ms, const D * const p, const D f, const IdxTy pitch, const IdxTy pf=2)
{
D const *  ptr=p;
const IdxTy n=ms.size();
for ( IdxTy i=0; i<n; ++i)
{
ms[i]=CoTy(ms[i].r(),*ptr*f);
ptr+=pitch*pf;
}

}


template <class Os, class Ty> void dump_modes(Os & os, const Ty & ions) const
{

const IdxTy N=size();
for (IdxTy i=0; i<N; ++i)
{
const ModeShape & m=m_shapes[i];
IdxTy atom=0;
for (IdxTy j=0; j<N; j+=3)
{
const D x=m[j].r();
const D y=m[j+1].r();
const D z=m[j+1].r();
const D l=sqrt(x*x+y*y+z*z);
const D rx=ions.get(atom,0);
const D ry=ions.get(atom,1);
const D rz=ions.get(atom,2);
const D r=sqrt(rx*rx+ry*ry+rz*rz);
const D dot=(x*rx+y*ry+z*rz)/l/r; 
os<<" MODE SHAPE "<<i<<" "<<atom<<" eigenv="<<m_values[i].mag()<<" l="<<l<<" r.l="<<dot<<CRLF;
++atom;
}

}

}

template <class Td> void export_eigenvalues(Td & dest)
{
const IdxTy N=size();
dest.set_size(N,1);
//for (int i=0; i<N; ++i) dest(i,0)=CoTy(W[i*2],W[i*2+1]);
for (IdxTy i=0; i<N; ++i) dest(i,0)=m_values[i];
}


void export_eigenvalues(D* re, D* im)
{
const IdxTy N=size();
//for (int i=0; i<N; ++i) dest(i,0)=CoTy(W[i*2],W[i*2+1]);
for (IdxTy i=0; i<N; ++i) {re[i]=m_values[i].r(); im[i]=m_values[i].i();} 
}




};

// this really is specialized to tcl but not really

class mjm_solvers : public mjm_solvers_typedefs
{
typedef mjm_solvers Myt;

typedef double Dd;
public:
//Perts m_perts;
// Hartrees/au^2/am -> mks finally s-1 to ev 
// amu-> grms/mole / 6e23 -> grams/atom * 1e-3 -> kg/atom
const D UC; // =6.023e23/1e-3/.529e-10/.529e-10*4.36e-18;
//const D hbar_ev=4.136e-15/6.28*1.03; // ??? kluge this is wrong
//const D hbar_ev=6.582119e-16*1.03; // ??? kluge this is wrong
const D HBAR_EV; //=6.582119e-16*1.0; // ??? kluge this is wrong
const D EV_INVERSE_CM; // 8065.5;
//ProbTy * m_p;
public:
mjm_solvers():
UC(6.023e23/1e-3/.529e-10/.529e-10*4.36e-18),
HBAR_EV(6.582119e-16*1.0), // ??? kluge this is wrong
 EV_INVERSE_CM(8065.5)
 
 {}
~mjm_solvers() { }







// input mat should be const?
//template <class Tmodes, class Tmat> static void 
//solve_complex(Tmodes & dest,Tmat & fc);
#if 1 
static int complex_eigens(int N,Dd * W, Dd * VR, Dd * VL, Dd * A);
template <class Tmodes, class Tmat> static void solve_complex(Tmodes & dest,Tmat & fc)
{
//char  JOBVL='N';
//char  JOBVR='V';
//char  JOBVR='N';
int N=fc.size();
D * A =(D*)fc.raw();
//info()<<MM_MARK; 
//m_port.dump_matrix(info());
double W[2*N];
//double * VL=0;
// complex eigenvalues return a real,imag pair of vectors LOL
double * VL= new double[N*N*2];
double * VR= new double[N*N*2];

int INFO=complex_eigens(N,W,VR,VL,A);

#if 0

int LDA=N;
int  LDVL=1;
int  LDVR=N;
//int LWORK=N*N*4;
int LWORK=N*4;
double WORK[LWORK*4];
double RWORK[N*4];
int INFO=0;
// these are actually omega SQUARED
info()<<MM_MARK<<" now calling cgeev_ "<<N<<CRLF;
// t5his no work
zgeev_(&JOBVL,&JOBVR,&N,A,&LDA,W,VL,&LDVL,VR,&LDVR,WORK,&LWORK,RWORK,&INFO);
const bool stupid=false&&true;
if ( stupid)
{
char  JOBVS='V';
char  SORT='N';
int * SELECT=0;
int  SDIM=0;
int LDVS=N;
double  VS[N*N*2];
double BWORK[1];
// c is wrong, this is z fff
zgees_(&JOBVS,&SORT,SELECT,&N,A,&LDA,&SDIM,W,VS,&LDVS, WORK,&LWORK, RWORK,BWORK,&INFO); 
}

#endif

info()<<MM_MARK<<" info= "<<INFO<<CRLF;
// what the ffff
//dest.load((simple_complex*)(W),IdxTy(N),IdxTy(1));
dest.add_eigenvalues(W,N);
dest.add_eigenvector_columns(VR,N);
/*
if ( 0)
{
dest.set_size(N,1);
//typedef simple_complex CoTy;
for (int i=0; i<N; ++i) dest(i,0)=CoTy(W[i*2],W[i*2+1]);
//info()<<MM_MARK<<" loaded ok  info= "<<INFO<<CRLF;
for (int i=0; i<N; ++i)
{
//const D ev=::sqrt(::fabs(dest(IdxTy(0),IdxTy(i)).r()*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
//const D cminverse=ev*EV_INVERSE_CM;
// info()<<MM_MARK<<" "<<i<<" "<<ev<<" ( "<<cminverse<<" ) "<<WR[i]<<" "<<WI[i]<<CRLF;
}

}
*/

// don't forget to put the eigens into dest LOL. see above
delete [] A;
delete [] VR;
delete [] VL;

}

#endif


//template <class Tmodes, class Tmat> 
//static void solve_real(Tmodes & dest,Tmat & fc);

#if 1 
template <class Tmodes, class Tmat> 
static void solve_real(Tmodes & dest,Tmat & fc,const StrTy & lbl="EIGENVECTORS")
 {


int N=fc.size();
info()<<MM_MARK<<" dumpding fc  step 1in solve_real "<<N<<CRLF;
info().flush();
double * A =fc.raw();
info()<<MM_MARK<<" dumpding fc in solve_real "<<N<<CRLF;
//fc.dump_matrix(info());
double WR[N],WI[N];
// complex eigenvalues return a real,imag pair of vectors LOL
double * VR= new double[N*N];
int INFO=real_eigens(N,WR,WI,VR,A);
// copied code
info()<<MM_MARK<<" dumping eigenvectors "<<INFO<<CRLF;
dest.add_eigenvalues(WR,WI,N);
dest.add_eigenvector_columns_rc(VR,N);
info()<<MM_MARK<<" done adding ok info= "<<INFO<<CRLF;
Tmat vrd;
vrd.load(VR,N,N);
vrd.dump_matrix(info(),lbl);
delete [] A;
delete [] VR;
} 
static int real_eigens(int N,Dd * WR, Dd * WI, Dd * VR, Dd * A);
#if 0
static int real_eigens(int N,Dd * WR, Dd * WI, Dd * VR, Dd * A)
{char  JOBVL='N';
char  JOBVR='V';
//char  JOBVR='N';
int LDA=N;
double * VL=0;
int  LDVL=1;
int  LDVR=N;
int LWORK=N*N;
double WORK[LWORK];
int INFO=0;
const bool stupid=false&&true;
//const bool stupid=true;
if ( stupid)
{
char  JOBVS='V';
char  SORT='N';
int * SELECT=0;
int  SDIM=0;
int LDVS=N;
//double  VS[N*N*2];
double BWORK[1];
// c is wrong, this is z fff
// these are actually omega SQUARED
dgees_(&JOBVS,&SORT,SELECT,&N,A,&LDA,&SDIM,WR,WI,VR,&LDVS, WORK,&LWORK, BWORK,&INFO); }
else { dgeev_(&JOBVL,&JOBVR,&N, A, &LDA,WR,WI,VL,&LDVL,VR,&LDVR,WORK, &LWORK , &INFO); }
return INFO;
}

/*
info()<<MM_MARK<<" info= "<<INFO<<CRLF;
dest.add_eigenvalues(WR,WI,N);
// tjhis faiols with complex ev
dest.add_eigenvector_columns_rc(VR,N);
info()<<MM_MARK<<" done adding ok info= "<<INFO<<CRLF;

*/

/*
if ( eireal!=0)
{
eireal->set_size(N,1);
typedef simple_complex CoTy;
for (int i=0; i<N; ++i) (*eireal)(i,0)=CoTy(WR[i],WI[i]);



}
for (int i=0; i<N; ++i)
{
const D ev=::sqrt(::fabs(WR[i]*UC))*HBAR_EV;
//const D cminverse=ev*8065.5;
const D cminverse=ev*EV_INVERSE_CM;
 info()<<MM_MARK<<" "<<i<<" "<<ev<<" ( "<<cminverse<<" ) "<<WR[i]<<" "<<WI[i]<<CRLF;
}
*/

/*
delete [] A;

//SloMatTy vrd;


info()<<MM_MARK<<" dumping eigenvectors "<<INFO<<CRLF;
Tmat vrd;
vrd.load(VR,N,N);
vrd.dump_matrix(info(),lbl);


// try to examine the mode shapes in vrd in light of force constants
// and ion positions
//analyze_modes(m_port,m_ions,vrd);


delete [] VR;


}

*/


#endif


#endif


} ;


#endif


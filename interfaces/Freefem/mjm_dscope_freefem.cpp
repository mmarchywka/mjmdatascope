// -*- Mode : c++ -*-
//
// SUMMARY  : This code is an interface to mjm_datascope derived from the Freefem++ codebase as descirbed below.See the include file mjm_dscope_ff_iface.h for details. 
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

/*

 This file only has use as  part of Freefem++
See the Freefem++ site for details. 
 */


// TODO most of the mjm_ds stuff is not contained 
// but more needs to be done

#include "mjm_dscope_ff_iface.h"

// moved to interface file could of course be repseated here
#if 0
#include "ff++.hpp"
#include "AFunction_ext.hpp"
#endif

// right now not used in interfce 
#include "lgmesh.hpp"
using Fem2D::Mesh;
using Fem2D::MeshPoint;
extern bool NoWait;
typedef Mesh const * pmesh;

// I hate doing this but most of the names have a dscope 
// in them similar to old c code 
using namespace mjm_dscope_globals;

// This was a first attempt and should still work but really a mess 
class MjmDscope :  public E_F0 { public:
   typedef pmesh  Result;
   Expression getmesh,getmeshL;
   Expression names; // filename;
std::vector<Expression> m_exp;
// TODO  need to hide this in dscope
Ragged m_dscope_header;
Ragged m_evct_header;
//bool have_evct;

   MjmDscope(const basicAC_F0 & args)
    {
      args.SetNameParam();
      getmesh=to<pmesh>(args[0]);
      names=to<string*>(args[1]);

// TODO  remove old code in commentcs once tested  
dscope_mesh_header(m_dscope_header,"");
/*
Dscope::traits_type::Ss ss;
ss<<"# id test pattern chunks 0 "<<CRLF;
ss<<"# ffmesh 2D "<<CRLF;
//ss<<"# trace "<<CRLF;
//ss<<""<<CRLF;
m_dscope_header.load(ss,false);
*/

dscope_start();      
for(int i=2; i<args.size(); ++i) m_exp.push_back( to< double >(args[i]));
   }
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmesh>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MjmDscope(args);}
    AnyType operator()(Stack s) const ;

};

AnyType MjmDscope::operator()(Stack stack) const
{
  using  Fem2D::MeshPointStack;
  const  Fem2D::Mesh * Thh = GetAny<pmesh>((*getmesh)(stack));
   string * params =  GetAny<string*>((*names)(stack));
//   if (!xx && !yy ) 
if (true)
{
MM_ERR(" je suis dans "<<MMPR(__FUNCTION__))
Ragged r; 
r=m_dscope_header;
    const Mesh &Thtmp = *Thh;;
 int   nv = Thtmp.nv;
  int  nt = Thtmp.nt;
  int  nbe = Thtmp.neb;
// just write the mesh format
r.new_line();
r<<nv; r<<nt; r<<nbe; 
r.new_line();
    for (int ii = 0; ii < Thtmp.nv; ii++) {
      const Mesh::Vertex &vi(Thtmp(ii));
r<<vi.x; r<<vi.y; r<<vi.lab; 
r.new_line();

    }

int kkk=1;
    for (int ii = 0; ii < Thtmp.nt; ii++) {
      const Mesh::Triangle &vi(Thtmp.t(ii));
      int i0 = kkk+  Thtmp(ii, 0);
      int i1 = kkk+ Thtmp(ii, 1);
      int i2 = kkk+ Thtmp(ii, 2);
r<<i0; r<<i1; r<<i2; r<<vi.lab;
r.new_line();
    }
    for (int ii = 0; ii < Thtmp.neb; ii++) {
      const Mesh::BorderElement &vi(Thtmp.be(ii));    // 
      int i0 = kkk+ Thtmp.operator( )(vi[0]);
      int i1 =  kkk+Thtmp.operator( )(vi[1]);
r<<i0; r<<i1; r<<vi.lab;
r.new_line();
    }
dscope().send(r,DscopeBlock(),0);
// from marchywka@happy:/home/ubuntu/dev/freefem/FreeFem-sources-4.12/plugin/seq$ vi VTK_writer.cpp 
// TODO FIXME copy code for m_exp array here don't use evct
if (m_exp.size())
{
for(int ie=0; ie<m_exp.size(); ++ie)
{
int nbsol=nv;
  KN< double > valsol(nbsol);
  valsol = 0.;
  KN< int > takemesh(nbsol);
  takemesh = 0;
  MeshPoint *mp3(MeshPointStack(stack));

  for (int it = 0; it < nt; it++) {
    for (int iv = 0; iv < 3; iv++) {
      int i = (*Thh)(it, iv);
      mp3->setP(Thh, it, iv);
      valsol[i] = valsol[i] + GetAny< double >((*m_exp[ie])(stack));
      ++takemesh[i];
    }
  }

  for (int i = 0; i < nbsol; i++) {
    valsol[i] /= takemesh[i];
  }

} // ie 
r=m_evct_header;
dscope().send(r,DscopeBlock(),0);
} // have_evct


///////////////////////////////////////////////////////////

     }
   else {
MM_ERR(" je suis dans code mort  "<<MMPR(__FUNCTION__))
#if 0 
    MeshPoint *mp(MeshPointStack(stack)) , mps=*mp;
MM_ERR(" je suis dans "<<MMPR(__FUNCTION__))
     ofstream fp((*fn+".points").c_str());
     ofstream ff((*fn+".faces").c_str());
     fp.precision(12);
MM_ERR(" je suis dans "<<MMPR(__FUNCTION__))
     if (verbosity>1)
       cout << "  -- Opening files " << (*fn+".points") << " and " << (*fn+".faces") << endl;
    const   Fem2D::Mesh & Th=*Thh;
    long nbv=Thh->nv;
    long nbt=Thh->nt;
    ff << nbt << endl;
    fp << nbv << endl;
    KN<long> num(nbv);
    num=-1;
    int k=0;
    for (int it=0;it<nbt;it++)
      {
        const Fem2D::Triangle & K(Th[it]);

        num[Th(K[0])]=k++;
        num[Th(K[1])]=k++;
        num[Th(K[2])]=k++;

         ff << " 3 " << Th(K[0])+1 << ' ' << Th(K[1])+1 << ' ' << Th(K[2])+1 << ' '
         	<< " 0 0 0 " << K.lab << '\n';
      }
    if( verbosity>5)
     cout << "  - end writing faces  " << endl;
    for (int iv=0;iv<nbv;iv++)
      {
        const Fem2D::Vertex  & v(Th(iv));
        ffassert( iv == Th(num[iv]/3,num[iv]%3));
        mp->setP(Thh,num[iv]/3,num[iv]%3);

        fp << GetAny<double>((*xx)(stack)) << ' ';
        fp << GetAny<double>((*yy)(stack)) << ' ';
        if (zz)
         fp << GetAny<double>((*zz)(stack)) << ' ';
        else
         fp << " 0 ";
        fp << v.lab<< '\n';
      }
    if( verbosity>5)
     cout << "  - end writing points  " << endl;

   *mp= mps;
#endif
   }
  //  delete fn;   modif mars 2006 auto del ptr
   return SetAny<pmesh>(Thh);
}

////////////////////////////////////////////////////////////////

// this nominrally works but interace to edp code not stable
class MjmDscopeL :  public E_F0 { public:

// using at global or function scope works ... 
//using namespace mjm_dscope_globals;
   typedef pmesh  Result;
   Expression getmeshL;
   Expression names; // filename;
PMap  m_exp,m_str;
// TODO  need to hide this in dscope
Ragged m_dscope_header;
Ragged m_evct_header;
   MjmDscopeL(const basicAC_F0 & args)
    {
		args.SetNameParam();
      getmeshL=to<pmeshL>(args[0]);
      names=to<string*>(args[1]);

dscope_trace_header(m_evct_header,"");

dscope_start();

// exp[i0]= to< KN_<double> >(args[i0]);
//MM_ERR(" did the thing ")

dsacope_express(m_exp,m_str,args,2);

   } // compile time ctor
    static ArrayOfaType  typeargs() { return  ArrayOfaType(atype<pmeshL>(),atype<string*>(),true);}
    static  E_F0 * f(const basicAC_F0 & args){ return new MjmDscopeL(args);}
    AnyType operator()(Stack s) const ;

}; // MjmDscopeL
// home/ubuntu/dev/freefem/FreeFem-sources-4.12/plugin/seq$ vi vortextools.cpp 
// runtime code 
AnyType MjmDscopeL::operator()(Stack stack) const
{
//typedef unsigned int IdxTy;
typedef int IdxTy;
//MM_ERR(" je suis dans "<<MMPR(__FUNCTION__))
string * params =  GetAny<string*>((*names)(stack));
  using  Fem2D::MeshPointStack;
  const  MeshL & Th =* GetAny<pmeshL>((*getmeshL)(stack));
    typedef MeshL::Element Element;
    typedef MeshL::Vertex Vertex;
    int nt=Th.nt,nv=Th.nv;
// pv is the evaluated string in order but without position index 
DscopeParams pv;
dscope_params(pv,m_str,stack);

// send the mesh onnly 
if (m_exp.size()==0)
{
// TODO this does not work no header and no data points.. 
Ragged r;
// this params is also pv[0]
dscope_add_params(r,m_dscope_header,params);
    for(int k=0;k<nv;k++){
//        const Element & E(Th[k]);
      const Vertex &vi(Th(k));
//MM_ERR(MMPR4(k,vi.x,vi.y,vi.z))
}
dscope().send(r,DscopeBlock(),0);
} // mesh only 


if ((0!=m_exp.size()))
{

Ragged r;
if (pv.size()>1)
{
	StrTy id="freefem";
	BaseParams kvp(pv[1]);
	kvp.get(id,"id");
	//dscope_trace_header(m_evct_header,"");
	dscope_trace_header(r,id);

}
else
{
// r is specific to this call, can be modified etc. 
dscope_add_params(r,m_evct_header,params);
}

int nbsol=nv;
DscopeBlock b(nbsol,m_exp.size()+1);
 std::vector< KN< double > >  vvalsol(m_exp.size());
// made m_exp a map to retain the position info
auto ie=m_exp.begin();
MM_LOOP(ii,vvalsol) {(*ii).resize(nbsol); } 
MM_SZ_LOOP(ii,vvalsol,szvv)
{
  KN< double >& valsol=vvalsol[ii];
valsol=double(0);
  KN< int > takemesh(nbsol);
  takemesh = 0;
  MeshPoint *mp3(MeshPointStack(stack));

  for (int it = 0; it < nt; it++) {
    for (int iv = 0; iv < 2; iv++) {
      int i = Th(it, iv);
      mp3->setP(&Th, it, iv); 
      //valsol[i] = valsol[i] + GetAny< double >((*m_exp[ii])(stack));
      valsol[i] = valsol[i] + GetAny< double >((*((*ie).second))(stack));
      ++takemesh[i];
    }
  }

  	for (int i = 0; i < nbsol; i++) { valsol[i] /= takemesh[i]; }
	MM_SZ_LOOP(j,valsol,szv){ b(j,ii)=valsol[j];
}
++ie;
} // i 
    for(int k=0;k<nv;k++){
      	const Vertex &vi(Th(k));
		b(k,m_exp.size())=vi.x;
	} // for k 
dscope().send(r,b,0);
} // m_exp
   return SetAny<pmeshL>(&Th);
}


//////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////



static void init_mjmdscope() {
//  if(verbosity&&(mpirank==0) )  cout <<"datascope ";
// old 2 D mesh display 
  Global.Add("mjmdscope","(",new OneOperatorCode<MjmDscope>);
// display meshL oscilloscope format using 
  Global.Add("mjmdscopeL","(",new OneOperatorCode<MjmDscopeL>);
    
}
LOADFUNC(init_mjmdscope);

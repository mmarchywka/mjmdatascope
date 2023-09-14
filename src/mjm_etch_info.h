#ifndef MJM_ETCH_INFO_H__
#define MJM_ETCH_INFO_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"
#include "mjm_string_kvp.h"
//#include "mjm_dscope_rawfifo.h"
//#include "mjm_ff_mesh_adapt.h"
//#include "mjm_etchfront.h"
//#include "mjm_pawnoff.h"
//#include "mjm_blob.h"
//#include "mjm_affine_mesh_hood.h"
#include "mjm_collections.h"
//#include "mjm_mesh_tape.h"
#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>

/*
@software{,
  author = {Michael J Marchywka},
  title = {},
abstract=(),
institution={},
license={Knowledge sir should be free to all },
publisher={Mike Marchywka},
email={marchywka@hotmail.com},
authorid={orcid.org/0000-0001-9237-455X},
  filename = {},
  url = {},
  version = {0.0.0},
  date-started = {}
}
*/
/*
2023-07-06 updated mesa drivers to get medit to work better
and now making backup copy of this and take out junk code
to test the dof modifications. 
2023-06-26 Fixed IsIn to get stragglers and just reversed
fnew or fni with fwd flag depending on gradient although 
hard coded reversed for now. 
Now need to get boundary numbers for triangulated mesh.

2023-06-26 missing some triangles, etch curves not orented
so etching could be in wrong dir. Want to start with gradient 
to the right ( liquid to the left ). 

*/
// Sun 11 Jun 2023 08:52:22 AM EDT
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_ff_phase_bound   





mjm_global_credits::credit __credit__mjm_etch_info("mjm_etch_info" , "  ");

template <class Tr,class Te>
class mjm_etch_info 
{
// typedef mjm_ff_phase_bound Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;


// TYPEDEF 
typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_ragged_table Ragged;
typedef typename Ragged::Line Line ;
//typedef mjm_pawnoff<Tr>  Hand;
//typedef mjm_blob<Tr>  Blob;

//typedef mjm_affine_mesh_hood<Tr> MeshHood;

//typedef mjm_tape_set<Tr> TapeSet;
//typedef mjm_tape_vector<Tr> TapeVector;
// API
public:
std::vector<StrTy> dofs,fns;
mjm_etch_info(){m_stagger=0;}
mjm_etch_info(Te * pmom, const StrTy & msg)
{
MM_ERR(" setting up etch in fo "<<MMPR(msg))
m_mom=pmom;
m_eps=msg;
// this was a ssv line but use the query string format
BaseParams kvp(msg);
//StrTy 
m_dir="./play/";
m_stagger=0;
kvp.get(m_dir,"dir");
mesh_name="inputmesh";
mesh_fn="equsaven.msh";
kvp.get(mesh_fn,"mesh_fn");
m_base="equsaven";
kvp.get(m_base,"base");
kvp.get(m_stagger,"stagger");
mesh_fn=m_dir+mesh_fn;
dofs.push_back(StrTy("c"));
dofs.push_back(StrTy("nm"));
dofs.push_back(StrTy("pm"));
dofs.push_back(StrTy("v"));
dofs.push_back(StrTy("H"));
dofs.push_back(StrTy("Cl")) ;
MM_LOOP(ii,dofs){const StrTy fn=m_dir+m_base+(*ii)+".txt"; fns.push_back(fn);} 

} // ctor
const StrTy eps() const { return m_eps;}
const IdxTy stagger(const IdxTy n )  
{ MM_ERR(" setting from  "<<MMPR2(m_stagger,n)) m_stagger=n; return m_stagger;}
const IdxTy stagger() const {MM_ERR(MMPR(m_stagger))  return m_stagger;}
IdxTy read(const IdxTy flags)
{
MM_ERR(MMPR(mesh_fn));
m_mom->m_meshes.load(mesh_name,mesh_fn,0);
//IdxTy  read_dof( const StrTy  & dof, const StrTy & fn, const StrTy & nm, const IdxTy  & flags )
{MM_SZ_LOOP(i,dofs,dofsz) 
	{  
MM_ERR(MMPR(fns[i]));
m_mom->m_meshes.read_dof(dofs[i],fns[i],mesh_name,0); }} 
return 0; 
} // read

// find etch surface
// IdxTy isotrail(const StrTy &nm, const StrTy  & trailn, const StrTy & dofn, const D & v, const IdxTy flags)
// this is the OLD mesh 
IdxTy write(const IdxTy flags)
{
// write everything back out
MM_ERR(" stupid one " <<MMPR3(__FUNCTION__,mesh_name,flags))
m_mom->m_meshes.save(mesh_name,mesh_fn,0);
{ MM_SZ_LOOP(i,dofs,dofsz) 
{  m_mom->m_meshes.write_dof(dofs[i],fns[i],mesh_name,0);}} 
return 0; 
} // write
IdxTy write(const StrTy & nm, const StrTy & fn,  const IdxTy flags)
{
MM_ERR(MMPR4(__FUNCTION__,nm,fn,flags))
// write everything back out
m_mom->m_meshes.save(nm,m_dir+fn+".msh",0);
auto _fns=fns;
_fns.clear();
MM_LOOP(ii,dofs){const StrTy _fn=m_dir+fn+(*ii)+".txt"; _fns.push_back(_fn);} 
{ MM_SZ_LOOP(i,dofs,dofsz) 
{  m_mom->m_meshes.write_dof(dofs[i],_fns[i],nm,flags);}} 
return 0; 
} // write

// TODO  this was supposed to include problem specific junk 
IdxTy adopt(const StrTy & newmesh, const IdxTy flags)
{
MM_ERR(" calling stupid adopt "<<MMPR(newmesh))
m_mom->m_meshes.clear_adopt(newmesh);
 MM_SZ_LOOP(i,dofs,dofsz) 
{
//IdxTy  adopt_dof(const StrTy & nm , const StrTy & from, const StrTy & dof)
m_mom->m_meshes.adopt_dof(newmesh,mesh_name,dofs[i]);
} // i 
return 0;
} // adopt 

// IIRC inits no work on earlier compilers and not needed here 
StrTy mesh_name="inputmesh";
StrTy mesh_fn="equsaven.msh";
StrTy m_base="equsaven";
private:
Te * m_mom;
// string to initialize the etchfront params 
StrTy m_eps;
StrTy m_dir;
IdxTy m_stagger;
}; // mjm_etch_info


#endif // MJM_FF_PHASE_BOUND_H__ 

#ifndef MJM_BLOCKMATRIX_H__
#define MJM_BLOCKMATRIX_H__

#include "mjm_globals.h"
#include "mjm_generic_iterators.h"
// for column and eventual latex output 
#include "mjm_format_text.h"

#include <vector>
#include <sstream>
#include <string>
#include <cstring>

/*
Matrix with unlimited indexing power andpotential optimization
for power of 2 and predictable access patterns. Can be blocked
and maybe diestributed etc. 
Template parameters ex data type  ignored but could be used to fix nesting
levelts etc 
vi mjm_block_matrix.h 
now complains about gnu++11

g++ -DTEST_OUTER_MATRIX_MAIN__ -Wall -I.. -std=gnu++11  -Wl,-rpath,/home/marchywka/d/gcc/stage-4.9/usr/local/lib64/ -x c++ mjm_block_matrix.h 


 g++ -DTEST_OUTER_MATRIX_MAIN__ -Wall -I.. -std=gnu++11  -x c++ mjm_block_matrix.h


 2330  g++ -DTEST_BLOCK_MATRIX_MAIN__ -Wall -I.. -x c++ mjm_block_matrix.h 
 2331  ./a.out

*/
template <class Ty> class mjm_block_matrix 
{
typedef mjm_block_matrix Myt;
public:
typedef Ty D;
typedef unsigned int IdxTy;
typedef int Sint;
typedef std::string StrTy;
typedef std::stringstream SsTy;
typedef std::vector<D> PtTy;
typedef std::vector<IdxTy> LocBaseTy;
class LocationType : public LocBaseTy
{
typedef LocationType Myt;
typedef LocBaseTy Super;
public:
LocationType() : Super() {}
LocationType(const IdxTy n ) : Super(n) {}

StrTy to_string() const
{
SsTy ss;
ss<<"(";
StrTy sep="";
for (IdxTy i=0; i<size(); ++i) {ss<<sep<<(*this)[i]; sep=",";}
ss<<")";
return ss.str();
} 

}; // LocationType

//typedef std::vector<IdxTy> LocTy;
typedef LocationType LocTy;
typedef mjm_generic_itor<2> LoopItor2;
typedef LocTy location_type;

// API 

mjm_block_matrix(  )
: m_dims(dimwit(1)),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{
}
// is thesize of this is zero it will bomb 
mjm_block_matrix( const LocTy & dims)
: m_dims(dims),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{
}

mjm_block_matrix( const Myt & that)
: m_dims(that.m_dims),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(initalloc(that.m_ptr))
{

}

mjm_block_matrix( const D * data, const IdxTy dim )
: m_dims(dimwit(dim)),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{ for ( IdxTy i=0; i<dim; ++i) m_ptr[i]=data[i]; }

mjm_block_matrix( const IdxTy dim )
: m_dims(dimwit(dim)),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{ }



mjm_block_matrix( const IdxTy i,const IdxTy j )
: m_dims(dimwit(i,j)),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{}
mjm_block_matrix( const IdxTy i,const IdxTy j, const IdxTy k  )
: m_dims(dimwit(i,j,k)),m_sizes(elements(m_dims))
,m_top(m_dims.size()-1),m_ptr(valloc())
{}

LocTy dimwit(const IdxTy i )
{ LocTy x; x.push_back(i); return x; }

LocTy dimwit(const IdxTy i, const IdxTy j  )
{ LocTy x; x.push_back(i); x.push_back(j); return x; }

LocTy dimwit(const IdxTy i, const IdxTy j , const IdxTy k  )
{ LocTy x; x.push_back(i); x.push_back(j); x.push_back(k); return x; }

void diag(const D & v )
{
if (m_sizes.size()==0) return; 
const IdxTy sz=m_sizes[0];
MM_ILOOP(i,sz) m_ptr[i]=0;
IdxTy i=0;
IdxTy p=0;
while ((p=idx(i,i))<sz) { m_ptr[p]=v;++i; } 
} // diag
// should just exclude one but more general to define a plane 
//void rotation(const IdxTy exax, const D & theta )
void rotation(const IdxTy ex1,const IdxTy ex2, const D & theta )
{
diag(1); // this gets most of them.. 
const D ct=cos(theta);
const D st=sin(theta);
(*this)(ex1,ex1)=ct;
(*this)(ex2,ex2)=ct;
(*this)(ex2,ex1)=st;
(*this)(ex1,ex2)=-st;


//MM_ILOOP(i,ny())// row 
//{
//MM_ILOOP(j,nx()) // col 
//{
// }
// } // i 

} // diag





~mjm_block_matrix() { 
// this is calling dtor on elements wtf 
delete [] m_ptr; m_ptr=0; }

Myt & operator=(const std::vector<D> & that ) 
{
const IdxTy sz=that.size();
resize(sz);
for (IdxTy i=0; i<sz; ++i) (*this)(i)=that[i];

return *this;
}

Myt & operator=(const Myt & that) { this->from(that); return * this; } 
// copy from 
void from( const Myt& that)
{
 	m_dims=that.m_dims;
	m_sizes=elements(m_dims);
	m_top=(m_dims.size()-1);
	realloc(that.m_ptr);
}

template <class Td> void copy_to(Td & d)
{
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) d(i)=m_ptr[i];
}
template <class Td> void copy_to_libmesh(Td & d)
{
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) d.set(i,m_ptr[i]);
}


template <class Td> void copy_to2(Td & d)
{
LoopItor2 itor(m_dims);
while (itor) { d(itor[0],itor[1])=(*this)(itor[0],itor[1]); ++itor; }
}

template <class Td> void copy_from(const Td & d)
{
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) m_ptr[i]=d(i);
}

template <class Td> void copy_from2(const Td & d)
{
LoopItor2 itor(m_dims);
while (itor) { (*this)(itor[0],itor[1])=d(itor[0],itor[1]); ++itor; }
}
// these hokey things need a stride or iterator with stride or
// some way to slice etc. 
void assign_column(const Myt & data, const IdxTy col)
{
const IdxTy sz=data.size();
for (IdxTy i=0; i<sz; ++i ) (*this)(i,col)=data(i);
}
void axby_column(const Myt & data, const IdxTy col,const D & c1, const D & c2)
{
const IdxTy sz=data.size();
for (IdxTy i=0; i<sz; ++i ) (*this)(i,col)=data(i)*c1 + (*this)(i,col)*c2;
}

void copy_column( Myt & data, const IdxTy col) const
{
const IdxTy sz=data.size();
for (IdxTy i=0; i<sz; ++i ) data(i)=(*this)(i,col);
}

D dot_cols(const IdxTy i, const IdxTy j)
{
D x=0;
MM_ILOOP(k,ny()) { x+=(*this)(k,i)*(*this)(k,j); }
return x;
}
class _triple_type
{
typedef std::stringstream Ss; 
public: 
_triple_type() :i(0),j(0),v(0){}
_triple_type( const IdxTy _i,const IdxTy _j, const  D&  _v ) 
	:i(_i),j(_j),v(_v){}
bool operator<(const _triple_type & that) const { return v<that.v;}
StrTy dump() const { Ss ss; ss<<MMPR3(i,j,v); return ss.str(); } 
IdxTy i,j; D v;
}; // _triple_type

typedef _triple_type triple_type;
static bool Bit(const IdxTy f, const IdxTy b) { return (0!=(f&(1<<b))); } 
template <class Tt> void  to_triples( Tt & tt, const IdxTy flags)
{
const bool non_zed=Bit(flags,0);
const bool upper_tria=Bit(flags,1);
const bool lower_tria=Bit(flags,2);
const bool exclude_diag=Bit(flags,3);
const IdxTy i0=lower_tria?(exclude_diag?(1):0):0;
const IdxTy iif=ny(); //lower_tria?(exlude_diag?(i-1):i):ny();
for(IdxTy i=i0; i<iif; ++i) 
{
const IdxTy j0=upper_tria?(exclude_diag?(i+1):i):0;
const IdxTy jf=lower_tria?(exclude_diag?(i-1):i):nx();
for(IdxTy j=j0; j<jf; ++j) 
{
const D& v = (*this)(i,j);
if ((v!=0)||!non_zed) tt.push_back(triple_type(i,j,v));
} // j 
} // i 


} // to_triples 

Myt transpose() const
{
Myt x=*this;
// this should just reverse the dims and then have an iterator
// with variable length locations that traverses in fwd and reverse
// order but for now just do the simple thing.. 
LoopItor2 itor(m_dims);
LocTy rev(m_dims.size());
reverse(rev,m_dims);
x.reform(rev);

while (itor)
{
reverse(rev,itor.cursor());
x(rev)=(*this)(itor.cursor());
++itor;
} // itor 

return x;
}

template <class Tl >
//static void reverse(LocTy & d, const LocTy & s)
static void reverse(LocTy & d, const Tl & s)
{
const IdxTy sz=s.size();
const IdxTy sz1=sz-1;
for (IdxTy i=0; i<sz; ++i){
//MM_MSG(MMPR4(d[sz1-i],s[i],i,d.size()))
 d[sz1-i]=s[i];
}


}

void resize (const IdxTy size)
{ LocTy dim; dim.push_back(size); resize(dim); }

void resize (const IdxTy i,  const IdxTy j)
{ 
	LocTy dim;
 	dim.push_back(i);
 	dim.push_back(j);
 	resize(dim); 
}


void resize( const LocTy & dims)
{
	delete [] m_ptr; 
 	m_dims=dims;
	m_sizes=elements(m_dims);
	m_top=(m_dims.size()-1);
	m_ptr=(valloc());
}



bool form_is( const LocTy & dims) const { return (dims==m_dims);}
bool form_is (const IdxTy i,  const IdxTy j)
{ 
	LocTy dim;
 	dim.push_back(i);
 	dim.push_back(j);
	return form_is(dim); 
}




void reform( const LocTy & dims)
{
//	delete [] m_ptr; 
 	m_dims=dims;
	const IdxTy oldsz=m_sizes[0];
	m_sizes=elements(m_dims);
	const IdxTy _top=(m_dims.size()-1);
//	m_ptr=(valloc());
//	if (m_top!=_top) MM_ERR(" reofrming changes size from "<<m_top<<" t " <<_top)
	m_top=_top;
	if (oldsz!=m_sizes[0]) 
		MM_ERR(" reofrming changes size from "<<oldsz<<" t " <<m_sizes[0])
}

const IdxTy nx() const { return m_dims[m_top]; }
const IdxTy ny() const {if (m_top==0) return 0;  return m_dims[m_top-1]; }


const LocTy&  form() const { return m_dims; }
// I think this is the number of columns 
const IdxTy  dims() const { return m_dims[m_top]; }
const IdxTy  n_dims() const { return m_dims.size(); }
/// actually this should probably be size(i)
const IdxTy & dim(const IdxTy i) const { return m_dims[i]; }
const IdxTy & size(const IdxTy i) const { return m_dims[i]; }

const LocTy&  sizes() const { return m_sizes; }
const IdxTy  size() const { if (m_sizes.size()!=0) return m_sizes[0]; return 0;  }
D& operator[]( const LocTy & idx) { return m_ptr[pos(idx,m_sizes)]; }
// these things could be used to masquerade as vectors with small loss
// of index computation speed although just easier to make temps
//D& operator[]( const IdxTy & idx) { return m_virtual+m_stride*idx; }
D& operator[]( const IdxTy & idx) { return m_ptr[idx]; }
const D& operator[]( const IdxTy & idx) const { return m_ptr[idx]; }
D& operator()( const LocBaseTy & idx) { return m_ptr[pos(idx,m_sizes)]; }
const D& operator()( const LocBaseTy & idx)const { return m_ptr[pos(idx,m_sizes)]; }
D& operator()( const IdxTy & idx) { return m_ptr[idx]; }
const D& operator()( const IdxTy & idx)const { return m_ptr[idx]; }
IdxTy idx(const IdxTy & i, const IdxTy & j) const
{ return j+i*m_sizes[m_top]; }
//D& operator()( const IdxTy & i, const IdxTy & j){ return m_ptr[j+i*m_sizes[m_top]]; }
//const D& operator()( const IdxTy & i, const IdxTy & j)const { return m_ptr[j+i*m_sizes[m_top]]; }
D& operator()( const IdxTy & i, const IdxTy & j){ return m_ptr[idx(i,j)]; }
const D& operator()( const IdxTy & i, const IdxTy & j)const 
{ return m_ptr[idx(i,j)]; }

D& operator()( const IdxTy & i, const IdxTy & j, const IdxTy & k)
{ return m_ptr[k+j*m_sizes[m_top]+i*m_sizes[m_top-1]]; }
const D& operator()( const IdxTy & i, const IdxTy & j, const IdxTy &k )const 
{ return m_ptr[k+j*m_sizes[m_top]+i*m_sizes[m_top-1]]; }
const D& operator()( const IdxTy & i, const IdxTy & j, const IdxTy &k,const IdxTy &l  )const 
{ return m_ptr[l+k*m_sizes[m_top]+j*m_sizes[m_top-1]+i*m_sizes[m_top-2]]; }
D& operator()( const IdxTy & i, const IdxTy & j, const IdxTy &k,const IdxTy &l  ) 
{ return m_ptr[l+k*m_sizes[m_top]+j*m_sizes[m_top-1]+i*m_sizes[m_top-2]]; }

D& rel(const IdxTy & base, const IdxTy &dim, const IdxTy & offset)
{ // kluge, need a "1" pushed on top
if (dim==0) return m_ptr[base+offset]; 
return m_ptr[base+m_sizes[m_top-dim]*offset]; 
}
static LocTy loc(const IdxTy * p,const IdxTy sz)
{
LocTy x;
// const IdxTy * pptr=p;
//while (pptr!=0){ x.push_back(*pptr); ++pptr;}
for (IdxTy i=0; i<sz; ++i) x.push_back(p[i]);
return x;
}

//private:
LocTy elements(const LocTy & in) const
{
	const IdxTy sz=in.size();
	if (sz==0) return LocTy();
	LocTy szv(sz);
	szv[sz-1]=in[sz-1];
	for (Sint i=sz-2; i>=0; --i) szv[i]=szv[i+1]*in[i];	
return szv;
}
const IdxTy pos(const LocBaseTy & idx, const LocTy & szs) const
{
if (idx.size()<1) return 0;
IdxTy p=idx[m_top];
for (IdxTy i=0; i<m_top; ++i) p+=idx[i]*szs[i+1];

#ifdef PARANOID
if (true)
{
for (IdxTy i=0; i<szs.size(); ++i) 
	if (idx[i]>=szs[i]) 
MM_ERR(" out of range "<<i<<" "<<idx[i] <<" vs "<<szs[i])

}
#endif

return p;
}
// this need only cary from the fieldincremented 
Sint inc( LocTy & idx, const IdxTy del) const { return inc(idx,del,m_top); } 
Sint inc( LocTy & idx, const IdxTy del, const IdxTy imin) const
{
if (idx.size()<1) return -1;
idx[imin]+=del;
//carry(idx);
carry(idx,imin);
// this fails on last indexin rnage 
if (idx[0]>=m_dims[0]) return 1; 
return assert_range(idx,__LINE__);
}

void carry( LocTy & idx) const
{
IdxTy sz=idx.size();
// could fail on entry, not exit 
assert_range(idx,__LINE__);
// not suer what to do with incomplete index array 
if (sz==0 ) return; // fing unsigned fudd 
for (IdxTy i=m_top; i>0; --i)
{
// dangerous with unsigned lol i also zero 
while (idx[i]>=m_dims[i]) { ++idx[i-1]; idx[i]-=m_dims[i]; } 

}
assert_range(idx,__LINE__);

}
void carry( LocTy & idx, const IdxTy imin) const
{
IdxTy sz=idx.size();
// could fail on entry, not exit 
// and in fact it was lol wtf. 
//assert_range(idx,__LINE__);
// not suer what to do with incomplete index array 
if (sz==0 ) return; // fing unsigned fudd 
for (IdxTy i=imin; i>0; --i)
{
bool carry=0;
// dangerous with unsigned lol i also zero 
while (idx[i]>=m_dims[i]) { carry=true; ++idx[i-1]; idx[i]-=m_dims[i]; } 
if (!carry ) break; 
}
// this fails on last index of range. 
//assert_range(idx,__LINE__);

}




template <class Ta> Sint  assert_range( LocTy & idx, const Ta & lable="") const
{
Sint rc=0;
IdxTy sz=idx.size();
if (sz!=m_sizes.size()){
 MM_ERR(lable<<" location error "<<m_sizes.size()<<" vs "<<sz)
rc=1; }
for (IdxTy i=0; i<sz; ++i)
{
// dangerous with unsigned lol i also zero 
if (idx[i]>=m_dims[i])
{
MM_ERR(lable<<" index error for "<<i<<" value "<<idx[i]<<" exceeds "<<m_dims[i])
 rc+= (1<<i);
} 

}

return rc;
} //assert_range




const IdxTy postwo(const LocTy & idx, const LocTy & szs) const
{
if (idx.size()<1) return 0;
IdxTy p=idx[0];
for (IdxTy i=1; i<idx.size(); ++i) p+=idx[i]<<szs[i-1];
return p;
}
// nneds to have D comparable to zero 
bool assert_zero() const
{
const IdxTy sz=size();
for (IdxTy i=0; i<sz; ++i) if (m_ptr[i]!=0) 
{
	MM_ERR(" non zero at "<<i<<" val "<< m_ptr[i])
	return false;
}

return true;
}

D * valloc() const 
{
if ( m_dims.size()==0) return 0;
//MM_ERR(info_string())
//MM_ERR(" size is "<<m_sizes[0])
const IdxTy sz=m_sizes[0];
//D * p= new D[m_sizes[0]]; 
D * p= new D[sz]; 
// this was great until "D" was an object lol 
//::memset(p,0,sizeof(D)*m_sizes[0]);
for(IdxTy i=0; i<sz; ++i) p[i]=D();
return p;
}
D * initalloc(const D* that)  
{
const IdxTy sz=m_sizes[0];
//D * p= new D[m_sizes[0]]; 
D * p= new D[sz]; 
//::memcpy(p,that,sizeof(D)*m_sizes[0]);
//for(IdxTy i=0; i<sz; ++i) p[i]=D();
MM_ONCE(" the copy ctor for block matrix makes no sense",)
for(IdxTy i=0; i<sz; ++i) p[i]=that[i];
m_ptr=p;
return p;

}
D * realloc(const D* that)  
{
// TODO moving the delete after memcpy allow self copy
// but also increases peak memory usage.. 
delete [] m_ptr;
//MM_ERR(info_string())
//MM_ERR(" size is "<<m_sizes[0])
D * p= new D[m_sizes[0]]; 
::memcpy(p,that,sizeof(D)*m_sizes[0]);
m_ptr=p;
return p;
}

// not sure this does a COLUMN lol 
Myt & copy_column(const Myt & that,const IdxTy dest_row, const IdxTy src_row)
{
const IdxTy sz=m_dims[m_top];
const IdxTy thatsz=that.m_dims[that.m_top];
const IdxTy cnt=(sz<thatsz)?sz:thatsz;
MM_ERR(" copy_col "<<sz<<" "<<thatsz<<" "<<cnt<<" dest_row "<<dest_row<<" src "<<src_row)
copy_whole( m_ptr+(sz*dest_row), that.m_ptr+(thatsz*src_row), cnt);
return *this;
}
// should call copy ctor etc but ok most of time 
void copy_whole( D* dest, const D* src, const IdxTy cnt)
{ ::memcpy(dest,src,sizeof(D)*cnt); }

//D& operator()( const IdxTy & i, const IdxTy & j){ return m_ptr[j+i*m_sizes[m_top]]; }
void zero_row(const IdxTy row)
{
const IdxTy rs=m_sizes[m_top];
zero(row*rs,rs);

}
// this does not work with some target types like rationals
void zero() { ::memset(m_ptr,0,sizeof(D)*m_sizes[0]); }
void zero_robust() {for(IdxTy i=0; i<m_sizes[0]; ++i ) m_ptr[i]=0;  }
void set_robust(const D & v) {for(IdxTy i=0; i<m_sizes[0]; ++i ) m_ptr[i]=v;  }
void fillbyte(const int c ) { ::memset(m_ptr,c,sizeof(D)*m_sizes[0]); }
void zero(const IdxTy first, const IdxTy len) { ::memset(m_ptr+first,0,sizeof(D)*len); }
void set_to_index(const IdxTy idx)
{

LoopItor2 ii(m_dims);
while (ii)
{
(*this)(ii.cursor())=ii[idx];
++ii;

} //ii 

} // set_to_index



// /* FUDD = ot looks like a scoping problem NOT a bizarre operator pborlem fudd 
template <class Td>  explicit operator mjm_block_matrix<Td> ()
{
mjm_block_matrix<Td> x(m_dims);
LoopItor2 ii(m_dims);
while (ii)
{
x(ii.cursor())=Td((*this)(ii.cursor()));
++ii;
}
} //ii 
// */

StrTy info_string() const 
{
SsTy ss;
const IdxTy sz=m_dims.size();
for (IdxTy i=0; i<sz; ++i) ss<<"("<<i<<" "<<m_dims[i]<<" "<<m_sizes[i]<<")";

return ss.str();
}
StrTy loc_string(const LocTy & loc) const 
{
assert_range(loc);
SsTy ss;
const IdxTy sz=loc.size();
ss<<"( ";
for (IdxTy i=0; i<sz; ++i) ss<<loc[i]<<" ";
ss<<")";
return ss.str();
}

Myt & push_form( const LocTy form)
{ m_form_stack.push_back(m_dims); reform(form); return *this; }
Myt & push_form( const IdxTy i)
{ LocTy loc; loc.push_back(i); return push_form(loc); }
Myt & push_form( const IdxTy i, const IdxTy j)
{ LocTy loc; loc.push_back(i); loc.push_back(j);  return push_form(loc); }

Myt & pop_form( )
{ reform(m_form_stack.back()); m_form_stack.pop_back(); return *this; }

template <class Tother> StrTy  diffsq( const Tother & that, bool only_diff)
{
SsTy ss;
LocTy dimsx=m_dims;
LocTy dims2x2;
dims2x2.push_back(16);
dims2x2.push_back(16);
reform(dims2x2);
IdxTy i=0;
LocTy curse(2);
// this does nto handle zero size too well. 
do { 
auto x1=(*this)(curse);
auto x2=that(curse[0],curse[1]);
auto diff=x2-x1;
if ((diff!=0)||!only_diff)
{
ss<<i<<" "<<curse[0]<<" "<<curse[1]<<" "<<x1;
ss <<" "<<x2<<" "<<(x2-x1)<<CRLF;
}

 ++i; 
} while ((*this).inc(curse,1)==0);


reform(dimsx);
return ss.str();
}

// note that "that" is non-const for reformation which should be
// mutable... 
template <class Tys,class Tl> void A_times_x(Tys & d, Tl & that)
{
//Myt & left=*this;
LoopItor2 itor(dim(0),dim(1));
while (itor)
{
const IdxTy & r=itor[0];
const IdxTy & c=itor[1];
d(r)+=(*this)(r,c)*that(r);

++itor;
} // itor 

}
// 2024-04 added operator no idea why not before
Myt operator*(const Myt & that)const  { return times(that); } 
Myt times(const Myt & that) const 
{
const IdxTy dims=m_top+1;
const IdxTy tdims=that.m_top+1;

const IdxTy rd=dim(0);
const IdxTy cd=that.dim(1);
const IdxTy cx=dim(1);
const IdxTy rx=that.dim(0);
if (cx!=rx)
{
MM_ERR(" bad dimensions "<<rx<<" "<<cx<<" and "<<rd<<" "<<cd)
}
Myt x(rd,cd);
//LoopItor2 titor(that.dim(0),that.dim(1));
//LoopItor2 itor(dim(0),dim(1));
MM_ONCE("not implemente right maybe danger will robinson  ",)
for(IdxTy i=0; i<rd; ++i)
{
for(IdxTy j=0; j<cd; ++j)
{
D sum=0;
// row from this col from that 
for(IdxTy k=0; k<cx; ++k)
{
//MM_ERR(MMPR4(i,k,j,(*this)(i,k))<<MMPR(that(k,j)))
 sum+=(*this)(i,k)*that(k,j);
}
x(i,j)=sum;
} // j 
} // i 

return x; 
}


D sum() const
{
D x=0; 
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) x+=m_ptr[i];
return x;
}

D norm(const D & p=2.0) const
{
D x=0; 
const IdxTy sz=m_sizes[0];
// this could be much aster for l2 etc 
for (IdxTy i=0; i<sz; ++i) x+=::pow(m_ptr[i],p);
return ::pow(x,1.0/p);
}


Myt operator+(const D & val)
{
if (check_pathological()) { return Myt(); }
Myt x(m_dims);
//LoopItor2 itor(m_dims); whille (itor) { ++itor; }
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) x.m_ptr[i]=m_ptr[i]+val;
return x;
} // add scalar

Myt operator+(const Myt & that)
{
if (check_pathological()) { return Myt(); }
Myt x(m_dims);
//LoopItor2 itor(m_dims); whille (itor) { ++itor; }
const IdxTy sz=m_sizes[0];
if (sz!=that.m_sizes[0])
{ MM_ERR(" sizes do not match "<<sz<<" is not "<<that.m_sizes[0]) } 
for (IdxTy i=0; i<sz; ++i) x.m_ptr[i]=m_ptr[i]+that.m_ptr[i];
return x;
} // add scalar




Myt operator-(const Myt & that) const 
{
// if (check_pathological()) { return Myt(); }
Myt x(m_dims);
//LoopItor2 itor(m_dims); whille (itor) { ++itor; }
const IdxTy sz=(that.m_sizes[0]<m_sizes[0])?that.m_sizes[0]: m_sizes[0];
for (IdxTy i=0; i<sz; ++i) x.m_ptr[i]=m_ptr[i]-that.m_ptr[i];
return x;
} // add scalar




Myt operator*(const D & val)
{
if (check_pathological()) { return Myt(); }
Myt x(m_dims);
//LoopItor2 itor(m_dims); whille (itor) { ++itor; }
const IdxTy sz=m_sizes[0];
for (IdxTy i=0; i<sz; ++i) x.m_ptr[i]=m_ptr[i]*val;
return x;
} // mult scalar


bool check_pathological()
{
if (m_sizes.size()<1)
{ MM_ERR(" sizes is zero so irnoring operation") 
return true; }

return false;
}
Myt outer_product(const Myt & that) const
{ 
Myt x(outer_product(m_dims,that.m_dims));
LoopItor2 k(x.m_dims);
LoopItor2 i(m_dims);
while (i)
{
// ideally at least one of these would traverse in sequence allowing
// for popinter incements lol. 
LoopItor2 j(that.m_dims);
while (j)
{
// this does rely on itor order. 
x(k.cursor())=(*this)(i.cursor())*that(j.cursor());
++k;
++j;
} // j 
++i;
} // i 
return x; 
}
LocTy outer_product(const LocTy & x, const LocTy & y) const
{
const IdxTy sx=x.size();
const IdxTy sy=y.size();
const IdxTy sz=sx+sy;
LocTy z(sz);
for (IdxTy i=0; i<sx; ++i) z[i]=x[i];
for (IdxTy i=sx; i<sz; ++i) z[i]=y[i-sx];
return z; 
}
// for code generators, traverse and list all element values
void catalog()
{
const Myt & x=*this;
// the cursor LocTy differs, one has extra crap cBaseTy;
//std::map<StrTy, std::vector<LocTy> > values;
std::map<StrTy, std::vector<LocBaseTy> > values;
std::map<StrTy,StrTy> names;
const StrTy base="v";
IdxTy cnt=0;
LoopItor2 itor(m_dims);
while (itor) 
{
SsTy ss;
ss<<x(itor.cursor());
if (values.find(ss.str())==values.end())
{
SsTy sn;
sn<<base<<cnt;
++cnt;
names[ss.str()]=sn.str();
}
values[ss.str()].push_back(itor.cursor());
	++itor;
} // itor
std::cout<<generate(values,names)<<CRLF;


} //catalog


template <class Tn, class Tv>
StrTy generate(const Tv & values,Tn & names)
{
const IdxTy sz=m_dims[m_top];
const IdxTy szdim=m_top;
const Myt & x=*this;
const StrTy vtype="ValTy";
const StrTy casttype="ValTy";
std::stringstream ss;
for (auto ii=names.begin(); ii!=names.end(); ++ii)
{ 
auto fudd_const = values.find((*ii).first);
const IdxTy fudd=(*fudd_const).second.size();
ss<<" const "<<vtype<<" "<<(*ii).second<<" = "
<<casttype<<"("<<(*ii).first<<"); // n="<<fudd<<CRLF; }
MM_ERR("// dims "<<x.m_dims.size()<<" "<<x.m_top<<" "<<x.m_dims[0]);
LoopItor2 itor(x.m_dims);
while (itor)
{
for (IdxTy j=0; j<sz; ++j)
{
ss<<"mat(";
for (IdxTy i=0; i<szdim; ++i) { ss<<itor[i];  ss<<",";  }
SsTy sn;
sn<<x(itor.cursor());
const StrTy varname=names[sn.str()];
ss<<itor[x.m_top]<<")="<<varname<<"; ";
++itor; 
} // j
ss<<CRLF;
} //itor 
return ss.str();
} // generate

// move all elements cyclically by amount delta along dimension dim 
// dim zero woulc be row, dim 1 column etc. 
Myt cycle(const IdxTy dim, const int  delta) const
{
typedef int I;
Myt x(m_dims); // =*this;
LoopItor2 itor(m_dims);
if (itor.cursor().size()<=dim) 
{ MM_ERR( " cursor size seems wrong for dim "<<dim<<" "<<itor.cursor().size()) } 
const IdxTy modv=m_dims[dim];
while (itor)
{
// here D is the data type aarrg
const D val=(*this)(itor.cursor());
const IdxTy idx=itor[dim];
const IdxTy newv=(I(idx)+delta+I(modv))%modv;
itor[dim]=newv;
//MM_ERR(" fudd "<<idx<<" "<<modv<<" "<<newv<<" "<<val)
x(itor.cursor())=val;
itor[dim]=idx; 
++itor;
} // itor 

return x;

}  // cycle
StrTy to_string() const
{

SsTy ss;

LoopItor2 itor(m_dims);
while (itor)
{
ss<<itor.loc_string(itor.cursor())<<"="<<(*this)(itor.cursor());
++itor; 
if (!itor) ss<<CRLF;
else if (itor[m_top]==0) ss<<CRLF;
}

return ss.str(); 

} // to_string

// note also this requires a data tpe coercible into a double 
template <class Ts> // fudder does not convert quoted string 
StrTy to_string_cols(const StrTy & label, const IdxTy p, const IdxTy n
, const Ts sep=StrTy(" ")) const
{
typedef mjm_format_text FormTy;
SsTy ss;
ss.precision(p);
IdxTy cnt=0;
//const StrTy sep=" ";
LoopItor2 itor(m_dims);
while (itor)
{ // no sep here to save space 
if (cnt==0 ) { ss<<label<<""<<itor.loc_string(itor.cursor()); } 
//ss<<sep<<(*this)(itor.cursor());
ss<<sep<<FormTy::format_float((*this)(itor.cursor()),3,9);
++cnt;
++itor; 
if (!itor) ss<<CRLF;
else if (itor[m_top]==0){ cnt=0;  ss<<CRLF; }
else if (cnt==n) { cnt=0; ss<<CRLF; } 
}

return ss.str(); 

} // to_string_cols

StrTy dump(const IdxTy flags=0) const
{
SsTy ss;
if (m_sizes.size()==0) return ss.str();
const IdxTy sz=m_sizes[0];
//const IdxTy sz=m_dims[0];
MM_ILOOP(i,sz) { ss<<"["<<i<<"]="<<m_ptr[i];}
return ss.str();
}


// often want to combine mats for a given row. 
template <class Ts, class Tl > // fudder does not convert quoted string 
StrTy to_string_cols(const Tl & label, const IdxTy p, const IdxTy n
,const Myt & that, const Ts sep=StrTy(" "),const bool enable_format=true ) const
{
typedef mjm_format_text FormTy;
SsTy ss;
ss.precision(p);
IdxTy cnt=0;
//const StrTy sep=" ";
LoopItor2 itor(m_dims);
LoopItor2 itor_that(that.m_dims);
while (itor)
{ // no sep here to save space 
if (cnt==0 ) { ss<<label<<""<<itor.loc_string(itor.cursor()); } 
//ss<<sep<<(*this)(itor.cursor());
if (enable_format)
{ss<<sep<<FormTy::format_float((*this)(itor.cursor()),3,9);}
else { ss<<sep<<(*this)(itor.cursor()); }
++cnt;
++itor; 
const bool eol=((!itor)|| (itor[m_top]==0)|| (cnt==n));
if (eol)
{ cnt=0; 

while (itor_that)
{
if (enable_format)
{ss<<sep<<FormTy::format_float((that)(itor_that.cursor()),3,9);}
else {ss<<sep<<(that)(itor_that.cursor());}

++itor_that;
const bool eol_that=((!itor_that)|| (itor_that[that.m_top]==0));
if (eol_that) break;
} // itor_that

ss<<CRLF;

 } // eol 

} // itor 

return ss.str(); 
}


////////////////////////////////////////

template <class Ts, class Tl > // fudder does not convert quoted string 
StrTy to_string_cols_R(const Tl & label, const IdxTy p, const IdxTy n
,const Myt & that, const Ts sep=StrTy(" "),const bool enable_format=true ) const
{
typedef mjm_format_text FormTy;
SsTy ss;
ss.precision(p);
IdxTy cnt=0;
//const StrTy sep=" ";
LoopItor2 itor(m_dims);
LoopItor2 itor_that(that.m_dims);
while (itor)
{
if (cnt==0 ) { ss<<label<<" "<<itor.loc_string(itor.cursor(),false); } 
//ss<<sep<<(*this)(itor.cursor());
if (enable_format)
{ss<<sep<<FormTy::format_float((*this)(itor.cursor()),3,9);}
else { ss<<sep<<(*this)(itor.cursor()); }
++cnt;
++itor; 
const bool eol=((!itor)|| (itor[m_top]==0)|| (cnt==n));
if (eol)
{ cnt=0; 

while (itor_that)
{
if (enable_format)
{ss<<sep<<FormTy::format_float((that)(itor_that.cursor()),3,9);}
else {ss<<sep<<(that)(itor_that.cursor());}

++itor_that;
const bool eol_that=((!itor_that)|| (itor_that[that.m_top]==0));
if (eol_that) break;
} // itor_that

ss<<CRLF;

 } // eol 

} // itor 

return ss.str(); 
}


////////////////////////////////////////



//////////////////////////////////////////


//////////////////////////////////////////







////////////////////////

StrTy to_string_R(const StrTy & label ) const
{
SsTy ss;
const StrTy sep=" ";
LoopItor2 itor(m_dims);
while (itor)
{ // no sep here to save space 
 ss<<label<<" "<<itor.loc_string(itor.cursor());  
//ss<<sep<<(*this)(itor.cursor());
ss<<sep<<((*this)(itor.cursor()));
ss<<CRLF;
++itor; 
}

return ss.str(); 

} // to_string_cols

//////////////////////




StrTy to_min_max_string() const
{

SsTy ss;
LocBaseTy lmin, lmax;
D mmin=0,mmax=0;
bool first=true;
LoopItor2 itor(m_dims);
while (itor)
{
const D & v=(*this)(itor.cursor());
if (first||(v>mmax)) { mmax=v; lmax=itor.cursor(); }
if (first||(v<mmin)) { mmin=v; lmin=itor.cursor(); }
 first=false; 
++itor; 
}
ss<<" min "<<itor.loc_string(lmin)<<"="<<mmin;
ss<<" max "<<itor.loc_string(lmax)<<"="<<mmax;
 ss<<CRLF;

return ss.str(); 

} // to_string

// manual control needs to integrate with form change
void pad_cols_to(const IdxTy i)
{ while (m_column_names.size()<i) m_column_names.push_back(StrTy("."));
}
// stupid default is a ref to temp.. fudd 
const StrTy  col_name(const IdxTy i )  const
{
//	pad_cols_to(dims());
	if (i>=dims()) return "invalidcol";
	if (i>=m_column_names.size()) return ".";
	return m_column_names[i];
}
void col_name(const StrTy & nm, const IdxTy i)
{
	pad_cols_to(dims());
	if (i>=dims()) {MM_ERR( "invalidcol"<<i<<" versus "<<dims()) return;  } 
	m_column_names[i]=nm;

}
StrTy col_names() const
{
const StrTy sep=" ";
SsTy ss;
const IdxTy sz=m_column_names.size();
if (sz!=0) ss<<m_column_names[0];
for (IdxTy i=1; i<sz; ++i) ss<<sep<<m_column_names[i];
return ss.str();
}
//D& operator()( const IdxTy & i, const IdxTy & j){ return m_ptr[j+i*m_sizes[m_top]]; }
void index_elements(const  D & pfx, const IdxTy flags) 
{
const IdxTy sw=flags&15;
const IdxTy nc=m_dims[m_top];
const IdxTy nr=m_dims[m_top-1];
MM_ILOOP(i,nr) { 
MM_ILOOP(j,nc) { 
SsTy ssij;
switch (sw)
{
//various C codes... 
case 1: { ssij<<pfx<<"["<<i<<"]["<<j<<"]"; (*this)(i,j)=ssij.str(); break; } 
case 2: { ssij<<pfx<<"("<<i<<","<<j<<")"; (*this)(i,j)=ssij.str(); break; } 
// latex 
case 8: { ssij<<pfx<<"_{"<<i<<","<<j<<"}"; (*this)(i,j)=ssij.str(); break; } 

default:
ssij<<"_"<<i<<"_"<<j; (*this)(i,j)=pfx+ssij.str();
} // sw

} // j
} //i 
} // index_elements
Myt minor(const IdxTy _i, const IdxTy _j)
{
const IdxTy nc=m_dims[m_top]-1;
const IdxTy nr=m_dims[m_top-1]-1;
MM_ERR(MMPR2(nc,nr))
Myt x(nr,nc);
MM_ILOOP(i,nr) { 
//if (i==_i) continue;
const IdxTy ii=(i>=_i)?(i+1):i;
MM_ILOOP(j,nc) { 
//if (j==_j) continue;
const IdxTy jj=(j>=_j)?(j+1):j;
x(i,j)=(*this)(ii,jj);
} // j
} // i
return x;
} // minor
Myt invert3x3()
{
Myt x=Myt(3,3);
const D det=(*this)[0]*(*this)[4]*(*this)[8]
+(*this)[1]*(*this)[5]*(*this)[6]
+(*this)[2]*(*this)[3]*(*this)[7]
-(*this)[2]*(*this)[4]*(*this)[6]
-(*this)[1]*(*this)[3]*(*this)[8]
-(*this)[0]*(*this)[5]*(*this)[7];
if (det==0){  x.zero(); return x; } 
x(0,0)=(*this)[4]*(*this)[8]- (*this)[5]*(*this)[7]; 
x(0,1)=-((*this)[1]*(*this)[8]- (*this)[2]*(*this)[7]); 
x(0,2)=(*this)[1]*(*this)[5]- (*this)[2]*(*this)[4]; 

x(1,0)=-((*this)[3]*(*this)[8]- (*this)[5]*(*this)[6]); 
x(1,1)=(*this)[0]*(*this)[8]- (*this)[2]*(*this)[6]; 
x(1,2)=-((*this)[0]*(*this)[5]- (*this)[2]*(*this)[3]); 

x(2,0)=(*this)[3]*(*this)[7]- (*this)[6]*(*this)[4]; 
x(2,1)=-((*this)[0]*(*this)[7]- (*this)[1]*(*this)[6]); 
x(2,2)=(*this)[0]*(*this)[4]- (*this)[1]*(*this)[3]; 
const D di=1.0/det;
MM_ILOOP(i,9){  x[i]*=di; }

return x; 
}



// MEMBERS


// each index must be mod this entry
// index zero is mod m_dims[top] 1,m_dims[top-1] etc 
 LocTy m_dims;
 LocTy m_logs;
// these are the total sizes for multiply index -1
// - m_sizes[0] is the toal array size
// m_sizes top multpiples index 1 msizes_top*ind1+ind0 for example 
//  there is no const "1" on top however. 
 LocTy m_sizes;
IdxTy m_top;
D * m_ptr;

// this really needs a form OBJECT including header names 
typedef std::vector<LocTy> FormStack;
FormStack m_form_stack;
typedef std::vector<StrTy> ColNames;
ColNames m_column_names;
typedef std::vector<ColNames> NameStack;
NameStack m_name_stack;


} ; // mjm_block_matrix

#ifdef TEST_OUTER_MATRIX_MAIN__

int main(int argc, char ** argv)
{
typedef mjm_block_matrix<double> T;
//typedef T::location_type L;
const T::D  d1[]={1,2,3,4};

const T::D  d2[]={5,6,7,8};
T m1(d1,4), m2(d2,4);
MM_ERR(m1.to_string())
MM_ERR(m2.to_string())

MM_ERR("location")
MM_ERR((m1.outer_product(m2)).to_string())
return 0;
}  //main

#endif // TEST_OUTER_MATRIX_MAIN__


#ifdef TEST_BLOCK_MATRIX_MAIN__


int main(int argc, char ** argv)
{
typedef mjm_block_matrix<double> T;
typedef T::location_type L;
const T::IdxTy  szi[]={4,5,6,7};
const T::IdxTy  els=sizeof(szi)/sizeof(T::IdxTy);
MM_ERR("location")
L s=T::loc(szi,els);
MM_ERR("ctor")
T::IdxTy i=0;
T x(s);
MM_ERR("dtor")
L cursor(els);
do { x(cursor)=i; ++i; } while (x.inc(cursor,1)==0);
i=0;
L cursor2(els);
do { MM_ERR(i<<" "<<x(cursor2))  ++i; } while (x.inc(cursor2,1)==0);
x(1,2,3,4)=1234;
MM_ERR(" test "<<x(0,0,0,0)<<" "<<x(3,4,5,6)<<" "<<x(1,2,3,4))


return 0;
}


#endif // main



#endif


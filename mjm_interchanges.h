#ifndef MJM_INTERCHANGES_H__
#define MJM_INTERCHANGES_H__


#include <stdlib.h>
#include <stdio.h>
//#include <tcl.h>
#include <math.h>
#include <cmath>
#include <mjm_globals.h>
//#include <mjm_templates.h>
//#include <mjm_38spatial.h>
//#include <mjm_molecule.h>
#include <iostream>

//#include <cblas.h>
// not sure this works with linpack?
//#include <complex>

//#include "mjm_sample_spectrum.h"
#include "mjm_io.h"

class simple_complex
{
typedef simple_complex Myt;
typedef double DataTy ;
public:
simple_complex() :m_r(0),m_i(0) { }
simple_complex(const DataTy & r) :m_r(r),m_i(0) { }
simple_complex(const DataTy & r, const DataTy & i) :m_r(r),m_i(i) { }
simple_complex(const double * x) :m_r(x[0]),m_i(x[1]) { }
const DataTy & r() const { return m_r; }
const DataTy & i() const { return m_i; }
// really bad name?
const DataTy  abs() const { return m_r*m_r+m_i*m_i; }
const DataTy  mag() const { return ::sqrt(m_r*m_r+m_i*m_i); }
Myt operator*( const Myt & x) const
{ return Myt(x.r()*r()-x.i()*i(), x.i()*r()+x.r()*i()); }
// multi by conj of x
Myt conj( const Myt & x) const
{ return Myt(x.r()*r()+x.i()*i(), -x.i()*r()+x.r()*i()); }

Myt operator-( const Myt & x) const { return Myt(-x.r()+r(), -x.i()+i()); }
Myt operator-( ) const { return Myt(-r(),-i()); }

Myt & operator+=( const Myt & x)
{ m_r+=x.m_r; m_i+=x.m_i; return *this;  }

template <class Ty > Ty & pr(Ty & os) { os<<m_r<<" "<<m_i; return os; }
private:
DataTy m_r,m_i;

};





template <class Ty=double, class Td=int> 
class mjm_interchange_data
{
public:
typedef mjm_generic_traits Tr;
typedef Tr::StrTy StrTy;
typedef Tr::IsTy IsTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::OsTy OsTy;
typedef Td Derived;
//typedef double D;
typedef Ty  D;
typedef std::vector<D> P;
// just for ease for now  
// argh......
typedef std::vector<P> V;
static OsTy & info() { return std::cout; } 
protected:

V m_data;
//template <class Ty> Ty print(const Ty & x ) { return x;} 
public:
typedef D data_type;
typedef P array_type;
typedef V matrix_type;
//template <class Os> void print(const D & x) { Derived::print(x); } 
template <class Os> void print(Os & os, const D & x) { ((Derived*)this)->print(os,x); } 
//template <class Td, class Ts> void assign(Td& d, const Ts & s) { ((Derived*)this)->assign(d,s); } 
P & operator[](const IdxTy i) { return m_data[i];} 
const P & operator[](const IdxTy i) const { return m_data[i];} 

D & operator()(const IdxTy i, const IdxTy j) { return m_data[i][j]; }
const D & operator() (const IdxTy i, const IdxTy j) const { return m_data[i][j]; }


P & operator()(const IdxTy i) { return m_data[i]; }
const P & operator()(const IdxTy i) const { return m_data[i]; }

// probably do not want to do this.. 
Derived & push_back( const P & row ) 
{
m_data.push_back(row);
return *this;
} 

/*
void assignWhole(const D * const  d, const IdxTy m, const IdxTy n)
{
m_data.clear();
D  const * dptr=d;
for ( IdxTy i=0; i<n; ++i)
{
// could just keep using the old one...
P p;
for ( IdxTy j=0; j<m; ++j)
{
p.push_back(*dptr);
++dptr;
}
m_data.push_back(p);
}

}

IdxTy size() const { return m_data.size(); } 
*/


void load(const D * const  d, const IdxTy m, const IdxTy n)
{
m_data.clear();
D  const * dptr=d;
for ( IdxTy i=0; i<n; ++i)
{
// could just keep using the old one...
P p;
for ( IdxTy j=0; j<m; ++j)
{
p.push_back(*dptr);
++dptr;
}
m_data.push_back(p);
}

}

void load_lines(IsTy * is)
{
m_data.clear();
IdxTy line_size=0;
enum { MAXLINE=1<<20, MAXLEN=MAXLINE-2};
typedef char Chr;
typedef std::stringstream SS;
Chr line[MAXLINE];
while ((*is).good()&&!(*is).eof())
{
(*is).getline(line,MAXLEN);
// this was pickig up junk line, should also have uniform size
if (!((*is).good()&&!(*is).eof()) ) break;
SS ss(line);
// could just keep using the old one...
P p;
while ((ss).good()&&!(ss).eof())
{
D d;
ss>>d;
//if((ss).good()&&!(ss).eof()) 
p.push_back(d);
}
info()<<MM_MARK<<" pushing back size "<<p.size()<<" for line "<<m_data.size()<<CRLF;
if ( line_size==0) line_size=p.size();
if ( line_size!=p.size()) info()<<MM_MARK<<" size is not uniform "<<p.size()<<" vs "<<line_size<<" at "<<m_data.size()<<CRLF;

m_data.push_back(p);
}

}







IdxTy size() const { return m_data.size(); } 

////////////////////////////////////////////////////////////////////
void zero(P & x, const IdxTy sz, const D & zed)
{
x.clear();
for (IdxTy i=0; i<sz; ++i) x.push_back(zed);
}

void zero(V & x, const P & p, const IdxTy sz)
{
x.clear();
for (IdxTy i=0; i<sz; ++i) x.push_back(p);

}

D * raw() const 
{
// really this is dumb but this is a SLOW struct
const IdxTy n=m_data.size();
IdxTy sz=0;
for(IdxTy i=0; i<n; ++i) sz+=m_data[i].size();
std::cout<<MM_MARK<<" export to raw size "<<sz<<" "<<n<<CRLF;
std::cout.flush();
D * a= new D[sz];
//std::cout<<MM_MARK<<" new fails "<<sz<<" "<<n<<CRLF;

IdxTy ptr=0;
for(IdxTy i=0; i<n; ++i)
{
const P & p=m_data[i];
IdxTy m=p.size();
for(IdxTy j=0; j<m; ++j){ a[ptr++]=p[j];
//std::cout<<MM_MARK<<" "<<ptr<<CRLF; std::cout.flush();
} 
}
return a;

}

template <class Ts> void element_multiply(const Ts & that)
{
const IdxTy n=m_data.size();
if ( n==0) return; 
const IdxTy m=m_data[0].size();
for(IdxTy i=0; i<n; ++i)
{
P & p=m_data[i];
const P & q=that.m_data[i];
//IdxTy m=p.size();
for(IdxTy j=0; j<m; ++j){ D x= q[j];  p[j]=p[j]*x; } 

}
}

template <class Os> void dump_sparse(Os & os)
{
const IdxTy n=m_data.size();
for(IdxTy i=0; i<n; ++i)
{
const P & p=m_data[i];
const IdxTy m=p.size();

for(IdxTy j=0; j<m; ++j)
{
if ( p[j]!=0){ os<<MM_MARK<<" "<<i<<" "<<j<<" ";//<<p[j]<<CRLF;
print(os,p[j]);
os<<CRLF;

} // if
}
}
}

template <class Os,class Ts> void dump_matrix(Os & os, const Ts & lbl)
{
const IdxTy n=m_data.size();
for(IdxTy i=0; i<n; ++i)
{
const P & p=m_data[i];
const IdxTy m=p.size();
os<<MM_MARK<<" "<<lbl<<" ";
for(IdxTy j=0; j<m; ++j)
{
//if ( p[j]!=0) os<<MM_MARK<<" "<<i<<" "<<j<<" "<<p[j]<<CRLF;
//os<<" "<<p[j];
D xx=p[j];
//if ( xx>0) xx=::log(p[j]*9)/::log(10);
//else if ( xx<0) xx=-::log10(-p[j]*9)/::log(10);
//os<<" "<<(int)xx;
print(os,xx);
//os<<" "<<xx;
os<<" ";
}
os<<CRLF;

}
}


};



//template <class Ty=double> 
 class mjm_2D_data : public mjm_interchange_data<double, mjm_2D_data> // <Ty>
{
//typedef mjm_2D_data Myt;
typedef mjm_2D_data Myt;
public:

/*
*/
mjm_2D_data(){}
mjm_2D_data(const IdxTy m, const IdxTy n)
{
set_size(m,n); 
} 
template <class Ty> mjm_2D_data(const Ty & that )
{
const IdxTy n=that.size();
const IdxTy m=that[0].size();
set_size(n,m);
for ( IdxTy i=0; i<n; ++i)
{ IdxTy p=i*m;
for ( IdxTy j=0; j<m; ++j)
{
// this is another fudding vector 
m_data[i][j]=that[i][j];
++p;
}
}

}

void set_size(IdxTy m,IdxTy n)
{
P r;
zero(r,n,0.0);
zero(m_data,r,m);
}

template <class Os> void print(Os & os, const D & x) { os<<x;  } 

};

 class mjm_2D_data_uint : public mjm_interchange_data<unsigned int, mjm_2D_data_uint> // <Ty>
{
//typedef mjm_2D_data Myt;
typedef mjm_2D_data_uint Myt;
public:

/*
*/
mjm_2D_data_uint(){}
mjm_2D_data_uint(const IdxTy m, const IdxTy n)
{
set_size(m,n); 
} 
void set_size(IdxTy m,IdxTy n)
{
P r;
zero(r,n,0);
zero(m_data,r,m);
}

template <class Os> void print(Os & os, const D & x) { os<<x;  } 

};






 class mjm_2D_data_complex : public mjm_interchange_data<simple_complex,mjm_2D_data_complex> // <Ty>
{

typedef mjm_2D_data_complex Myt;
typedef D CoTy;
public:

void set_size(IdxTy m,IdxTy n)
{
P r;
zero(r,n,CoTy(0,0));
zero(m_data,r,m);
}
mjm_2D_data_complex() {}
mjm_2D_data_complex(const mjm_2D_data & that)
{
assign(that);
}
Myt operator=(const mjm_2D_data & that)
{
assign(that);
return *this;

}
void assign (const mjm_2D_data & that)
{
m_data.clear();
const IdxTy sz=that.size();
for (IdxTy i=0; i<sz; ++i)
{
P y;
const IdxTy szj=that[i].size();

for (IdxTy j=0; j<szj; ++j)
{
y.push_back(CoTy(that(i,j),0));
}
m_data.push_back(y);

}
//return *this;
}

template <class Os> void print(Os & os, const D & x) 
{ os<<x.r()<<" "<<x.i();  } 

 
};







//class mjm_2D_data : public mjm_2D_data_template<double> { typedef mjm_2D_data Myt; };

class mjm_xy_data : public mjm_interchange_data<double>
{
typedef mjm_xy_data Myt;

V m_data;
public:


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
void load(IsTy * is)
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




};


#endif


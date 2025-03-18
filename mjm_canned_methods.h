#ifndef MJM_CANNED_METHODS_H__
#define MJM_CANNED_METHODS_H__
 
#include "mjm_globals.h"


// #define USE_RATIONALS
#ifdef USE_RATIONALS
#define MJM_RATIONAL_BIGINT mjm_bigint_mmpr
#include "mjm_bigint_mmpr.h"
//#include "mjm_cursing_list.h"
// some of these pickup the math libs... 
// add for the multiplicity tests
#include "mjm_integrals.h"
//#include "../copied_code/gmp_factorize.h"
#include "mjm_rational.h"
#endif // USE_RATIONALS

// for the presence absence vector 
//#include "mjm_ordering.h"
#include "mjm_char_mat.h"
#include "mjm_part_iterators.h"
#include "mjm_iterator_base.h"

#include "mjm_data_model_error_log.h"
// need this for Alfredo lol
//#include "mjm_rational.h"
// TODO FIXME move there bin search etc 
//#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
// add day number to notes 
//#include "mjm_calendar.h"
#include "mjm_strings.h"
#include "mjm_svg_writer.h"


//#include "mjm_string_index.h"

#include "mjm_cli_ui.h"
#include "mjm_fasta_ii.h"

//#include "mjm_collections.h"
//#include "mjm_svg_writer.h"

#ifdef USE_RATIONALS
#include "mjm_rational_things.h"
#include "mjm_table_managers.h"
#endif



//3245  echo parse-biom-json /home/marchywka/d/zymo/in682.171203.zymo/demo.Bac16Sv34/otus/otu_table.biom 5 | ./mjm_biom_hdf5.out 2>xxx
//#include "mjm_biom_hdf5.h"
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
// rand()
#include <stdlib.h>
#include <stdint.h>

/*


2018-08-02 copy from mjm_linc_graph
*/
// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
//#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<x); 
#ifndef MM_DMEL
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
#endif
/*

*/


////////////////////////////////////////////////////////////////


namespace canned_methods_traits
{
class Tr {
public:
typedef unsigned int IdxTy;
typedef double D;
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef mjm_block_matrix<D> MyBlock;
typedef  data_model_error_log Dmel; 
//typedef unsigned int KeyCode;
//typedef unsigned __int128 KeyCode;
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
//typedef std::map<IdxTy, IdxTy > Locations;
//typedef std::set<IdxTy > Locations;
//typedef std::vector<IdxTy > Locations;
//typedef mjm_string_ordering Ordering;

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }
static int myatoix(const char * c) { 
if (*c=='#') 
{
return myatoi(StrTy("0x")+StrTy(c+1)); 
}
return ::strtol(c,0,0); 
}

//typedef mjm_sparse_matrix<D> MySparse;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


#if 0 

#endif

}; // trees_and_tables_traits
///////////////////////////////////////////////////////////////


class layout_blocks 
{
typedef  canned_methods_traits::Tr  Tr;
//typedef  linc_graph_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef layout_blocks Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock MyBlock;
//typedef Tr::Ordering Ordering;
typedef mjm_ragged_table Ragged;
typedef std::vector<StrTy> Strings;
typedef std::vector<D> Numbers;

class group_inverse
{
public:
group_inverse(): m_j(~0), m_seq(~0) {}
void add_location(const IdxTy j) { m_j=j; } 
void set_seq(const IdxTy j) { m_seq=j; } 
IdxTy index() const { return m_seq;}
IdxTy j() const { return m_j;}

IdxTy m_j, m_seq;
} ; // group_inverse

typedef group_inverse Jinv;

typedef std::map<StrTy, Jinv> JinvMap;

public:

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }

void setup_geometry(ReadWriteMap&  rwm)
{

m_tm=50+30; rwm.get("tm",m_tm);
m_bm=330;rwm.get("bm",m_bm);
m_lm=200;rwm.get("lm",m_lm);
m_rm=50;rwm.get("rm",m_rm);

m_catpitch=30;rwm.get("catpitch",m_catpitch);
m_ysz=600;rwm.get("ysz",m_ysz);
//m_ynames=.5*D(m_tm); rwm.get("ynames",m_ynames);
m_ynames=0; rwm.get("ynames",m_ynames);
m_namesz=10; rwm.get("namesz",m_namesz);
m_xstarts=10+m_lm*.5; rwm.get("xstarts",m_xstarts);
m_xend=m_rm*.5; rwm.get("xend",m_xend);
// ; y+=yinc; // m_namesz+2;
//D ymin=tm;
//D ymax=ymin+ysz;
m_yleg=m_tm+m_ysz+30;
m_szleg=m_catpitch*.9; rwm.get("szleg",m_szleg);
m_colorleg="red";
m_xs=IdxTy( m_catpitch*m_cats+.5)+m_lm+m_rm;
m_ys=m_ysz+m_tm+m_bm;
m_wr=m_catpitch;
m_hr=.01*m_ysz;
m_xmin=1e-4;
m_vmin=4;
m_vmax=m_vmin-1;

m_y_date_loc=m_ys-m_bm-20; 
D delyd=0;
rwm.get("deldatey",delyd);
m_y_date_loc+=delyd;
rwm.get("datey",m_y_date_loc);
MM_ERR(" ASSFUDD " <<MMPR4(m_y_date_loc,delyd,m_ys,m_bm))
m_ygrid=9;
m_xz=m_lm;
m_xf=m_lm+m_cats*m_catpitch;
m_yz=m_tm;
m_yf=m_tm+m_ysz;
m_xleg=.9*m_lm;

}
void default_xy_lab(ReadWriteMap & rwm)
{

//m_ylab= StrTy(" Range limited log ratio to mean");
m_ylab= StrTy("Limited log ratio to median"); rwm.get("ylab",m_ylab);
m_ylab_sz=.2*m_lm; rwm.get("ylab_sz",m_ylab_sz);
m_ylab_x=.2*m_lm; rwm.get("ylab_x",m_ylab_x);
m_ylab_y=.5*m_ysz+m_tm+m_ylab_sz*m_ylab.length()*.5*.5; rwm.get("ylab_y",m_ylab_y);

}
void default_xy_day_lab(ReadWriteMap & rwm)
{

//m_ylab= StrTy(" Range limited log ratio to mean");
m_ylab= StrTy("Normalized amounts "); rwm.get("ylab",m_ylab);
m_ylab_sz=.2*m_lm; rwm.get("ylab_sz",m_ylab_sz);
m_ylab_x=.2*m_lm; rwm.get("ylab_x",m_ylab_x);
m_ylab_y=.5*m_ysz+m_tm+m_ylab_sz*m_ylab.length()*.5*.5; rwm.get("ylab_y",m_ylab_y);



}




void default_rainbod(Strings & rainbow)
{
rainbow.push_back("red");
rainbow.push_back("blue");
rainbow.push_back("green");
rainbow.push_back("orange");
rainbow.push_back("black");

}
void colors_and_day_names(ReadWriteMap & rwm)
{
m_sample_order.clear();
IdxTy i=0;
MM_LOOP(ii,m_samples) { m_sample_order[(*ii).first]=i; ++i; } 

Strings rainbow;
default_rainbod( rainbow);
MM_SZ_LOOP(i,m_xlabels,xxx)
{
Ss ss;
StrTy xlabel="xlabel";
ss<<m_serial_to_date[(i+m_day_min)];
StrTy & n=m_xlabels[i];
if (n.length()==0) xlabel=ss.str();
else xlabel=n;
//MM_ERR(MMPR2(rwm.size(),name))
rwm.get(xlabel,xlabel);
m_xlabels[i]=xlabel;
}
if (false) { 
{
MM_SZ_LOOP(i,m_colors,xxx)
{
Ss sc;
StrTy col="color";
sc<<col<<i;
StrTy c=rainbow[i%rainbow.size()];
rwm.get(sc.str(),c);
m_colors[i]=c;
}
} // names 
IdxTy i=0;
MM_LOOP(ii,m_samples)
{
const StrTy & nameauto=(*ii).first;
m_names[i]=nameauto;
Ss sc;
StrTy col="color";
sc<<col<<":"<<nameauto;
rwm.get(sc.str(),m_names[i]);
++i;

} // samples

} // false
else
{
MM_LOOP(ii,m_samples)
{
const StrTy & nameauto=(*ii).first;
IdxTy i=m_sample_order[nameauto];
{
Ss sc;
StrTy col="color";
sc<<col<<i;
StrTy c=rainbow[i%rainbow.size()];
rwm.get(sc.str(),c);
m_colors[i]=c;
}
Ss sc;
m_names[i]=nameauto;
StrTy col="name";
sc<<col<<":"<<nameauto;
rwm.get(sc.str(),m_names[i]);
++i;


} // ii  samples 



} // true



} // colors_and_day_names


void colors_and_names(ReadWriteMap & rwm)
{
Strings rainbow;
default_rainbod( rainbow);
MM_SZ_LOOP(i,m_names,xxx)
{
Ss ss,sc;
StrTy name="name";
StrTy col="color";
ss<<name<<i;
sc<<col<<i;
StrTy & n=m_names[i];
if (n.length()==0) name=ss.str();
else name=n;
//MM_ERR(MMPR2(rwm.size(),name))
rwm.get(name,name);
m_names[i]=name;
StrTy c=rainbow[i%rainbow.size()];
rwm.get(sc.str(),c);
m_colors[i]=c;
} // names 


}

// xy plots with headers and limits derived from data scan
// 5855 17567 17567 2018-02-05 1  Greta B-6  m_n=1 m_ntot=1 m_tot=12.5 m_max=12.5 m_avg=12.5 m_avg_dose=0 m_freq=1 m_invalid=0
// ./mjm_linc_graph.out -cmd "read-ragged in xxx2a" -cmd "add-ragged p period 3" -cmd "read-ragged p lys.txt"  -cmd "snacks-txt-svg srxxx 1 in p" -quit

// ss<<dayfirst<<" "<<daylast<<" "<<lbl<<" "<<pereff<<" ";
// StrTy  binlab=ss.str();
// os<<bin<<sep<<binlab<<sep<<dog<<sep<<food<<sep<<ss.str()<<CRLF;
//StrTy lbl=sf.serial_to_date[daylast]


void setup_day_plot(Ragged & r, Ragged & pr)
{

MM_ERR("setups day plot ")
ReadWriteMap&  rwm= m_rwm;
pr.to_map(rwm);
// the thing needs a day serial string date and line info
//m_bin_idx=1;rwm.get("bin_idx",m_serial_idx);
m_bin_idx=1;rwm.get("bin_idx",m_bin_idx);
// TODO the serial index has to match the date used in the liner entry
// or else it will create offset date labels in the plot 
m_serial_idx=2;rwm.get("serial_idx",m_serial_idx);
m_date_idx=3;rwm.get("date_idx",m_date_idx);
m_sample_idx=5;rwm.get("sample_idx",m_sample_idx);
m_food_idx=6;rwm.get("food_idx",m_sample_idx);
m_value_idx=16;rwm.get("value_idx",m_value_idx);

m_plot_thick=5;rwm.get("plot_thick",m_plot_thick);
m_dates_label=" All dates approximate";
rwm.get("dates_label",m_dates_label);

const IdxTy prsz=pr.size();
for(IdxTy i=0; i<prsz; ++i)
{
const Ragged::Line& line=pr.line(i);
const IdxTy szl=line.size();
if (szl<1) continue;
const StrTy & cmd=line[0];
//MM_ERR(MMPR(cmd))
if (cmd=="foods") {
for(IdxTy j=1; j<szl; ++j){

 ++m_foods_to_plot[line[j]]; 
MM_ERR(" adding food "<<MMPR2(line[j],m_foods_to_plot.size()))
}
}

} // i<prsz
MM_LOOP(ss,m_foods_to_plot) { MM_ERR(" FUDDL SHOT "<<MMPR2((*ss).first,(*ss).second)) }

m_day_max=0;
m_day_min=0;
const IdxTy rsz=r.size();
for(IdxTy i=0; i<rsz; ++i)
{
const Ragged::Line& line=r.line(i);
const IdxTy szl=line.size();
if (false) { Ss ss; for(IdxTy k=0; k<szl; ++k) { ss<<MMPR2(k,line[k]); } MM_ERR(ss.str()) }
const StrTy& bin=(line[m_bin_idx]);
const IdxTy serial=myatoi(line[m_serial_idx]);
const StrTy& date=(line[m_date_idx]);
const StrTy& sample=(line[m_sample_idx]);
const StrTy& food=(line[m_food_idx]);
if (m_day_max<serial) { m_day_max=serial; if (m_day_min==0) m_day_min=serial; }
if (m_day_min>serial) m_day_min=serial;
m_date_to_serial[date]=serial;
m_serial_to_date[serial]=date;
if (m_foods_to_plot.find(food)==m_foods_to_plot.end()) continue;


const D value=atof(line[m_value_idx].c_str());
if (i==0) { m_y_max=value; m_y_min=value; }
if (value>m_y_max) m_y_max=value;
if (value<m_y_min) m_y_min=value;
//MM_ERR(MMPR2(serial,date))
++m_samples[sample];
++m_bins[bin];
} 
// this is the number of data points. 
m_cats=m_day_max-m_day_min+1;

setup_geometry(rwm);
default_xy_day_lab(rwm);
m_jmax=m_samples.size();
m_colors=Strings(m_samples.size());
m_names=Strings(m_samples.size());
m_xlabels=Strings(m_cats);
colors_and_day_names(rwm);
} // setup_day_plot

// the log ratio to median plots 
void setup(Ragged & r, Ragged & pr)
{
MM_ERR("setups")
//Ragged & r=m_r;
ReadWriteMap&  rwm= m_rwm;
pr.to_map(rwm);
//mjm_svg_writer sw;
//std::ofstream  os(fn);
m_text_idx=0;rwm.get("text_idx",m_text_idx);
m_jmax=0;rwm.get("jmax",m_jmax);
m_jfirst=1;rwm.get("jfirst",m_jfirst);
m_first_line=1; rwm.get("first_line",m_first_line);

if ( m_jmax==0) if (m_first_line!=0)
{
if (r.size()==0) return; 
const Ragged::Line& line=r.line(0);
const IdxTy szl=line.size();
if (szl==0){ MM_ERR(" header is zer length danger will robinson") }
m_jmax=line.size()-1;
//for (IdxTy j=0; j<m_jmax; ++j)
}


m_min_line=m_jfirst+m_jmax;
r.min_line_size(m_min_line);

m_cats=r.size()-m_first_line;
const bool use_geom=true;
if (use_geom)
{
setup_geometry(rwm);
default_xy_lab(rwm);
}
else{ } //use_geom 

// angroups 
m_colors=Strings(m_jmax);
m_names=Strings(m_jmax);

// otu's
m_xlabels=Strings(m_cats);
// if there is a header start with that 
if (m_first_line!=0)
{
const Ragged::Line& line=r.line(0);
const IdxTy sz=line.size();
for (IdxTy j=0; j<m_jmax; ++j)
{
const IdxTy idx=j+m_jfirst;
if (idx>=sz) break;
const StrTy & vs= line[idx];
m_names[j]=vs;
}
} // 


m_x_start=Numbers(m_jmax);
m_x_frac=Numbers(m_jmax);
m_x_height=Numbers(m_jmax);
D xstart_def=0; rwm.get("xstart",xstart_def);
D xfrac_def=1; rwm.get("xfrac",xfrac_def);
for (IdxTy j=0; j<m_jmax; ++j) { 
D xi=xstart_def;
D xf=xfrac_def;
D xh=1;
StrTy nm=m_names[j];
if (nm.c_str()[0]==0) { Ss ss; ss<<j; nm=ss.str(); }
StrTy xin="xin-"+nm;
StrTy xsn="xsn-"+nm;
StrTy xhn="xhn-"+nm;
 rwm.get(xin,xi);
 rwm.get(xsn,xf);
 rwm.get(xhn,xh);
m_x_start[j]=xi; m_x_frac[j]=xf; 
m_x_height[j]=xh;
} 

for (IdxTy j=0; j<m_jmax; ++j) {m_jinv[m_names[j]].add_location(j);  }
//for (IdxTy j=0; j<m_jmax; ++j) 
IdxTy i=0; MM_LOOP(ji,m_jinv) { (*ji).second.set_seq(i); ++i; } 
//IdxTy i=0; MM_LOOP(ji,jinv) { (*ji).second=i; ++i; } 

IdxTy disperse=0; rwm.get("disperse",disperse);
if (disperse==1)
{
//std::map<StrTy,IdxTy> jinv; 
const IdxTy nuniq=m_jinv.size();
const D del=1.0/nuniq;
for (IdxTy j=0; j<m_jmax; ++j) {
const StrTy & nm=m_names[j];
const IdxTy idx=m_jinv[nm].index();
D xi=del*idx; 
D xf=del; 
D xh=1;
// still allow over ride 
StrTy xin="xin-"+nm;
StrTy xsn="xsn-"+nm;
StrTy xhn="xhn-"+nm;
 rwm.get(xin,xi);
 rwm.get(xsn,xf);
 rwm.get(xhn,xh);
m_x_start[j]=xi; m_x_frac[j]=xf; 
m_x_height[j]=xh;

//m_colors[j]=
}



}


MM_ERR("setups")
colors_and_names(rwm);
MM_ERR("setups")
} // setup


void find_day_statistics(Ragged & r)
{
// this needs to order the samples and locat the days. 
m_so=r.ordering();
const IdxTy sz=m_samples.size();
const IdxTy rsz=r.size();
m_dense_data.resize(m_cats,sz);
MM_LOOP(ii,m_foods_to_plot)
{ m_dense_map[(*ii).first].resize(m_cats,sz); }
//m_sample_order.clear();
IdxTy i=0;
//MM_LOOP(ii,m_samples) { m_sample_order[(*ii).first]=i; ++i; } 

i=m_day_min;
for(; i<=m_day_max; ++i)  { m_point_order[i]=i-m_day_min;  } 

//MM_LOOP(ss,m_foods_to_plot) { MM_ERR(" FUDDL SHOT "<<MMPR2((*ss).first,(*ss).second)) }

//std::vector<IdxTy>  som(m_jmax);
//for(IdxTy k=0; k<sz; ++k) som[k]=k;
m_vmin=0;
m_vmax=0;
for(IdxTy i=0; i<rsz; ++i)
{

const Ragged::Line& line=r.line(i);
//const IdxTy szl=line.size();
const IdxTy serial=myatoi(line[m_serial_idx]);
//const StrTy& date=(line[m_date_idx]);
const StrTy& sample=(line[m_sample_idx]);
const StrTy& date=(line[m_date_idx]);
const StrTy& food=(line[m_food_idx]);
//_if (food==StrTy("lysinehcl")) 
//{MM_ERR(" ASDFASFASFASDF "<<MMPR(food))}
if (m_foods_to_plot.find(food)==m_foods_to_plot.end()) continue;
//MM_ERR("FUDDING FOO D"<<MMPR3(food,date,m_foods_to_plot.size()))
const D value=atof(line[m_value_idx].c_str());
if (i==0) { m_vmax=value; m_vmin=value; }
if (value>m_vmax) m_vmax=value;
if (value<m_vmin) m_vmin=value;
if (m_sample_order.find(sample)==m_sample_order.end())
{
MM_ERR(" not found sample "<<MMPR4(date,food,sample,m_sample_order.size())) 
}
const IdxTy n=m_sample_order[sample];
if (m_point_order.find(serial)==m_point_order.end())
{
MM_ERR(" mising point order "<<MMPR4(serial,sample,date,food)<<MMPR3(m_day_min,m_day_max,m_point_order.size()))
}
const IdxTy x=m_point_order[serial];
// this needs to be sparse or find invalid like negative 
if (false) {
m_dense_data(x,n)+=value; // atof(line[m_value_idx]);
m_sparse_data[n][x]+=value;
}
else{ 
m_dense_map[food](x,n)+=value; // atof(line[m_value_idx]);
m_sparse_map[food][n][x]+=value;
}


}
MM_ERR(" ASSFUDD  shot "<<MMPR4(m_vmin,m_vmax,m_day_min,m_day_max))
} // find_day_statistics


void find_statistics(Ragged & r)
{
MM_ERR("stats")
m_so=r.ordering();
const IdxTy sosz=m_so.size();
m_merit=Numbers(sosz+0*m_cats);
m_avg=Numbers(sosz+0*m_cats);
m_dense_data.resize(m_cats,m_jmax);

const IdxTy ll=m_cats+m_first_line;
for(IdxTy i=m_first_line; i<ll; ++i)
{
const IdxTy ip=i-m_first_line;
D sum=0;
D n=0;
D mmax=0;
D mmin=0;
const Ragged::Line line=r.line(i);
if (line.size()>=4) {
for (IdxTy j=0; j<m_jmax; ++j)
{
//MM_ERR(" SHOT FURCJK RRRR FUDD "<<MMPR3(i,j,line.size()))
const StrTy & vs= line[m_jfirst+j];
const D dv=atof(vs.c_str());
//MM_ERR(MMPR4(vs,dv,i,j))
m_dense_data(ip,j)=dv;
if (dv==0) continue;
if (n==0) { mmax=dv; mmin=dv; }
if (dv>mmax) mmax=dv;
if (dv<mmin) mmin=dv;
sum+=dv;
n+=1;
} //j 
} // if 
//if ( n!=0) m_avg[ip]=sum/n; else m_avg[ip]=0;
//if ( n==m_jmax) m_merit[ip]=mmax/mmin;
//else m_merit[ip]=-mmax;
std::vector<IdxTy>  som(m_jmax);
for(IdxTy k=0; k<m_jmax; ++k) som[k]=k;
std::sort(som.begin(),som.end(),[&](const IdxTy& a, const IdxTy& b) 
{
return m_dense_data(ip,a)<m_dense_data(ip,b);
});

// this is actually now the MEDIAN 
if ( n!=0) m_avg[i]=sum/n; else m_avg[i]=0;
IdxTy loc=m_jmax>>1;
 m_avg[i]= m_dense_data(ip,som[loc]);
if (m_avg[i]==0)
{
MM_ERR(" moving median to non zero ")
while ((loc<m_jmax)&&(m_avg[i]==0))
{
 m_avg[i]= m_dense_data(ip,som[loc]);
++loc;
}

} // avg==0  

if ( n==m_jmax) m_merit[i]=mmax/mmin;
else m_merit[i]=-mmax;

} //i 


//const IdxTy ll=m_cats+m_first_line;
for(IdxTy i=m_first_line; i<ll; ++i)
{
const Ragged::Line line=r.line(i);
if (line.size()<=m_text_idx) continue; // fudd 
const StrTy & text= line[m_text_idx]; // line[text_idx];
m_xlabels[i-m_first_line]=text;
} // i 

} // find_statistics 


void sort(Ragged & r)
{
MM_ERR("sort")
r.sort(m_first_line,m_so,[&](const IdxTy& a, const IdxTy& b) -> bool
{
if (m_merit[a]<0)
{
if (m_merit[b]<0) return m_merit[a]<m_merit[b];
return !false;
}
if (m_merit[b]<0)
{
//if (merit[a]<0) return merit[a]>merit[b]; 
return false;
}


return m_merit[a]>m_merit[b];

});

} // sort

void write_svg_day_best_labels(std::ostream & os, mjm_svg_writer& sw, const IdxTy flags)
{
std::vector<IdxTy> bestpos(m_jmax);
std::vector<D> bestd(m_jmax);
if ( m_foods_to_plot.size()==0){
MM_ERR(" starved for foods lol")
 return;
}
const StrTy & food=(*(m_foods_to_plot.begin())).first;
//auto & dd=m_dense_data;
auto & dd=m_dense_map[food];
//auto & dd=m_dense_data;
auto & sd=m_sparse_map[food];

for(IdxTy i=0; i<m_cats; ++i)
{

for(IdxTy j=0; j<m_jmax; ++j)
{

//const D& vthis=m_dense_data(i,j);
const D& vthis=dd(i,j);
if (vthis==m_vmin) continue;
D  closest=2*(m_vmax-m_vmin)+m_vmin;
for(IdxTy k=0; k<m_jmax; ++k)
{
if (j==k) continue;
//const D& vthat=m_dense_data(i,k);
const D& vthat=dd(i,k);
const D diff=(vthat>vthis)?(vthat-vthis):(vthis-vthat);
if (diff<closest) closest=diff;
} //k 
if (closest>bestd[j]) {bestd[j]=closest; bestpos[j]=i;}
} // j
} /// i
MM_LOOP(is,m_samples)
{
int  i=m_sample_order[(*is).first];
 int x=bestpos[i];
//auto & sdi=m_sparse_data[i] ;
auto & sdi=sd[i] ;
while (sdi.find(x)==sdi.end()) { if (x==0) break ; x-=1; }



int del=1;
if (int(x+del)>=int(m_cats)) del=-1;
if (del>0) { 
while (sdi.find(x+del)==sdi.end()) { if (int(x+del+1)==int(m_cats)) break ; del+=1; }}
if (sdi.find(x+del)==sdi.end()) { del=-1; }
if (del<0) { while (sdi.find(x+del)==sdi.end()) { if ((x+del)==0) break ; del-=1; }}


//D v = m_dense_data((x),i);
D v = dd((x),i);
//D vdel = m_dense_data((x)+del,i);
D vdel = dd((x)+del,i);
v=(v-m_vmin)/(m_vmax-m_vmin);
vdel=(vdel-m_vmin)/(m_vmax-m_vmin);
MM_ERR(MMPR4((*is).first,v,vdel,x)<<MMPR(del))
const D xi=m_lm+x*m_catpitch;
const D xdel=m_lm+(x+del)*m_catpitch;
const D y=m_tm+(1.0-v)*m_ysz;
const D ydel=m_tm+(1.0-vdel)*m_ysz;
const D & x0=xi;
const D & y0=y;
const D & x1=xdel;
const D & y1=ydel;
MM_ERR(MMPR4(x0,y0,x1,y1)<<MMPR4(m_tm,m_ysz,m_lm,m_catpitch))
const D sz=m_plot_thick*2;
const StrTy & text=m_names[i];
const StrTy & col=m_colors[i];
//os<<sw.label_line_text_upoff(text+(*is).first, x0,y0,  x1,  y1, sz,  m_plot_thick*2, col, "center",.5) ; 
os<<sw.label_line_text_upoff(text, x0,y0,  x1,  y1, sz,  m_plot_thick*2, col, "center",.5) ; 

} // is

}
// pass first 10 then 101 aka 2 and  5
void write_svg_day_data(std::ostream & os, mjm_svg_writer& sw, const IdxTy flags)
{
Ss ss;
const bool label_foods=true;
const bool histogram_like=!true;
//const bool write_name=((flags&(1<<0))!=0);
//const bool write_rect=((flags&(1<<1))!=0);
//const bool skip_first_pass=((flags&(1<<2))!=0);
//const IdxTy first=skip_first_pass?1:0;
//std::map<IdxTy, IdxTy> olaps;
//const IdxTy sz=m_samples.size();
//const IdxTy nx=m_cats;
MM_LOOP(ff,m_foods_to_plot)
{
const StrTy & food=(*ff).first; // (*(m_foods_to_plot.begin())).first;
MM_ERR(" writing food "<<MMPR2(food, m_foods_to_plot.size()))
//auto & dd=m_dense_data;
//auto & dd=m_dense_map[food];
//auto & dd=m_dense_data;
auto & sd=m_sparse_map[food];


IdxTy i=0; 
IdxTy ifudd=0;
//MM_LOOP(ii,m_sparse_data)
MM_LOOP(is,m_samples)
//for(IdxTy i=0; i<sz; ++i)
{
//auto ii=m_sparse_data.find(i);
//if (ii==m_sparse_data.end()) continue;
if (m_sample_order.find((*is).first)==m_sample_order.end())
{
MM_ERR(" not found "<<MMPR3(i,(*is).first, m_sample_order.size()))
}
i=m_sample_order[(*is).first];
MM_ERR(MMPR4((*is).first,i,ifudd,m_names[i]))
++ifudd;
//auto  ii=m_sparse_data.find(i);
auto  ii=sd.find(i);
//if (ii==m_sparse_data.end()) { ++i; continue;}
if (ii==sd.end()) { ++i; continue;}
auto & s=(*ii).second;
std::vector<D> xv,yv;
MM_LOOP(jj,s)
//for(IdxTy j=0; j<nx; ++j)
{
//auto jj=s.find(j);
//if (jj==s.end()) continue;
const D ip=(*jj).first-0; // m_day_min;
const D x=m_lm+ip*m_catpitch;
//D v=m_sparse_data[(*ii).first][(*jj).first];
D v=sd[(*ii).first][(*jj).first];
v=(v-m_vmin)/(m_vmax-m_vmin);
if ((v>1)||(v<0))
{
//D vfudd=m_sparse_data[(*ii).first][(*jj).first];
D vfudd=sd[(*ii).first][(*jj).first];
MM_ERR(" v oor "<<MMPR4((*ii).first, (*jj).first, v,vfudd)<<MMPR2(m_vmax,m_vmin))
}
const D y=m_tm+(1.0-v)*m_ysz;
if (xv.size()>0){ 
D xlast=xv.back();
D ylast=yv.back();
if ((x!=xlast)||(y!=ylast)) { xv.push_back(x); yv.push_back(y);// }
if ((label_foods) &&((int(ip)%5)==3)){
ss<<sw.label_line_text(food,xlast,ylast,x,y,2,"white","center",0) ; 
ss<<CRLF;
}}}


if (ip==0) { xv.push_back(x); yv.push_back(y);}

D ipnext=ip+1; // m_day_min;
if ( histogram_like) { 
auto jjn=jj; ++jjn;
if ( (jjn)!=s.end())  ipnext=(*(jjn)).first;
}
const D xbar=m_lm+(ipnext)*m_catpitch;
xv.push_back(xbar);
yv.push_back(y);

} // j 
const StrTy col=m_colors[i];



// write the path at once. 
//StrTy vector_text( const Tvec& x, const Tvec& y,  const StrTy & color="#00ffffff",const D & sw=10, const D & opaq=.8) const
// TODO FIXME size could be 1, need to accomodate 
if (xv.size()!=0) { os<<sw.vector_text(xv,yv,col,m_plot_thick,.7);  
const bool label_mid_segs=!true;
if ( label_mid_segs) { 
const StrTy & text=m_names[i];
const IdxTy ix=(xv.size()-1)>>1;
const D & x0=xv[ix];
const D & y0=yv[ix];
const D & x1=xv[ix+1];
const D & y1=yv[ix+1];
const D sz=m_plot_thick*2;
os<<sw.label_line_text(text+(*is).first, x0,y0,  x1,  y1, sz,  col, "center",3) ; }
}
++i;
} // i 
} // ff 
if (label_foods) os<<ss.str();
#if 0 
#endif
//////////////////////////////////////////


 write_svg_day_best_labels( os, sw, flags);
}








// pass first 10 then 101 aka 2 and  5
void write_svg_data(std::ostream & os, mjm_svg_writer& sw, const IdxTy flags)
{
const bool write_name=((flags&(1<<0))!=0);
const bool write_rect=((flags&(1<<1))!=0);
const bool skip_first_pass=((flags&(1<<2))!=0);
const bool write_density=((flags&(1<<3))!=0);
const bool  write_name_on_plot =((flags&(1<<4))!=0); //  true; 
const IdxTy first=skip_first_pass?1:0;
// group_name, otu/etc, div, density 
std::map< StrTy, std::map<IdxTy,  std::map<int,D> > > rho;
IdxTy ndiv=100;
//std::map<IdxTy, IdxTy> olaps;
for(IdxTy pass=first; pass<2; ++pass)
{
for (IdxTy j=0; j<m_jmax; ++j)
{
const StrTy & cj=m_colors[j];
const StrTy & name=m_names[j];
const IdxTy ll= m_cats+m_first_line;
for(IdxTy i=m_first_line; i<ll; ++i)
{
const IdxTy ip= i-m_first_line;
const IdxTy is=m_so[i];
const D x=m_lm+ip*m_catpitch;
const D xrs=m_x_start[j];
const D xf=m_x_frac[j];

//const Ragged::Line line=r.line(is);
//MM_ERR(" SHOT FURCJK RRRR FUDD "<<MMPR3(i,j,line.size()))
//if (line.size()<=(m_jfirst+j)) continue; 
//const StrTy & vs= line[jfirst+j];
const D vraw=m_dense_data(is-m_first_line,j); // atof(vs.c_str());
D v=m_xmin;
D avgg=m_avg[is-0*m_first_line];
if (avgg<=0)
{
const StrTy & text=m_names[j];
MM_ERR(" bad average ( or MEDIAN )  "<<MMPR4(is,avgg,j,text)<<MMPR(m_xlabels[is-m_first_line]))
avgg=1;
}
if (vraw!=0) v=log(vraw/avgg)/log(10);
//if (vraw!=0) v=log(vraw/m_avg[is-0*m_first_line])/log(10);

//else  v=log(v/m_avg[is-m_first_line])/log(10);
//MM_ERR(MMPR4(i,is,vraw,v))
if (pass==0)
{
if (vraw!=0) { 
if ( m_vmax<m_vmin) { m_vmax=v; m_vmin=v;
//MM_ERR(MMPR4(i,is,vraw,v))

 }
else {
if ( v>m_vmax) m_vmax=v;
if (v<m_vmin) m_vmin=v;
//MM_ERR(MMPR4(i,is,vraw,v)<<MMPR2(m_vmax,m_vmin))
} // max<min
 } // vraw
 
} // pass==0
if (pass==1) { 
v=(v-m_vmin)/(m_vmax-m_vmin);
// v is between zero and 1 
const D y=m_tm+(1.0-v)*m_ysz;
const D xr=x+xrs*m_wr;
const D xw=m_wr*xf;
const D txr=xr+.1*m_wr;
const D tyr=y+.75*m_hr;
if (write_density)
{
	int div=v*ndiv;
	if ( vraw==0) div=~0;	
	++rho[name][ip][div];

} // write_density 

// TODO overlaps ???? 
//if ( vraw!=0) os<<sw.shade_rect_text(x,y,wr,hr,cj,.5);
//if ( vraw!=0){ if (write_rect) { os<<sw.shade_rect_text(x+xrs*m_wr,y,m_wr*xf,m_hr,cj,.8); }
if ( vraw!=0){ if (write_rect) { os<<sw.shade_rect_text(xr,y,xw,m_hr,cj,.8);

os<<CRLF;
 }
if (write_name) {
 StrTy  text=m_names[j];
Ss ss; ss<<" "<<vraw<<" "<<m_avg[is]<<" ";
text=text+ ss.str();
//os<<sw.gtext_text(text,x+.1*m_wr,y+.75*m_hr,2, StrTy("white"),  StrTy("start") ,0); 
os<<sw.gtext_text(text,txr,tyr,2, StrTy("white"),  StrTy("start") ,0); os<<CRLF;
}
}
else
{
if (write_rect) { 
const D y=m_tm+(1.0-0)*m_ysz;
 //os<<sw.shade_rect_text(x,y,wr,hr,cj,.5);
 os<<sw.shade_rect_text(x,y,m_wr,m_hr,cj,.8);
os<<sw.gtext_text(StrTy("*"),x+(D(j+1)/(m_jmax+2))*m_wr , y, m_szleg, cj,  StrTy("middle") ,0); 
os<<CRLF;
}
}
} // pass ==0 

} // 
} // 
if (pass==1)
{
if (write_density)
{
D dy=m_ysz/ndiv;
//std::map<StrTy, IdxTy> jinv;
//for (IdxTy j=0; j<m_jmax; ++j)
//{
//jinv[m_names[j]]=j;
//}
MM_LOOP(rr,m_jinv)
{
IdxTy j=(*rr).second.j();
const StrTy & name=(*rr).first; // m_names[j];

const StrTy & cj=m_colors[j];
//const IdxTy ll= m_cats+m_first_line;
const D xrs=m_x_start[j];
const D xf=m_x_frac[j];
const D xh=m_x_height[j];
auto & m=rho[name];
for(IdxTy i=0; i<m_cats; ++i)
{
const IdxTy ip= i;
//const IdxTy is=m_so[i];
const D x=m_lm+ip*m_catpitch;
const D xr=x+xrs*m_wr;
const D xw=m_wr*xf;
auto & m2=m[ip];
IdxTy cmax=0;
IdxTy cmin=0;
if ( write_name_on_plot)
{
D tyr=m_tm+(1.0-0)*m_ysz-4;
const D fudd_kluge=1;
os<<sw.gtext_text(name,xr+fudd_kluge,tyr,2, StrTy("black"),  StrTy("end") ,90); 
os<<sw.gtext_text(name,xr+fudd_kluge,m_tm+2,2, StrTy("black"),  StrTy("start") ,90); 
os<<CRLF;
}

MM_LOOP(mm,m2){ 

if ((*mm).second>cmax) cmax=(*mm).second;
if ((*mm).second<cmin) cmin=(*mm).second;
if (cmin==0) cmin=(*mm).second;

 }
if (cmin<1)
{
MM_ERR( " bad histogram "<<MMPR4(cmin,cmax,name,ip))

}
const D lcmax=log(cmax);
const D lcmin=log(cmin);
const D lrange=(lcmax>lcmin)?(1.0/(lcmax-lcmin)):0;
const D del=1.0/ndiv;
MM_LOOP(mm,m2)
{
const int div=(*mm).first;
const IdxTy cnt=(*mm).second;
const D v=del*div; // 1.0*div/ndiv;
D y=m_tm+(1.0-v)*m_ysz;
D yend=y-dy*xh;
if (div<0) {y=m_tm+m_ysz+2; yend=y+.01*m_ysz*xh; } 
D vopq=1.0*cnt/cmax;
if (lrange!=0)  vopq=(log(cnt)-lcmin)*lrange;
if (vopq<.1) vopq=.1;
if (vopq>1) vopq=1;
os<<sw.shade_rect_text(xr,y,xw,yend-y,cj,vopq); os<<CRLF; 
} // mm 

} // i

} // j

} // write_density 

} // pass==1

if (pass==0)
{
//Ss ss,sr;
//ss<<std::setprecision(2)<<(m_vmax*1.1);
//ss>>vmax;
// this is negative need to make it fudding more negative fudd 
//sr<<std::setprecision(2)<<(m_vmin*1.1);
//sr>>vmin;
// expand the range by about 10 percent
//D r=m_vmax-m_vmin;
// want vmax=vmax+.1r;
D e= 1.1; // (m_vmax!=0)?(1.0+.1*r/m_vmax):0;
MM_ERR(" STARTING RANGE "<<MMPR2(m_vmax,m_vmin)) 
exp_range(m_vmax,true,e);
//e= (m_vmin!=0)?(1.0+.1*r/m_vmin):0;
exp_range(m_vmin,!true,e);
MM_ERR(" ENDING RANGE "<<MMPR2(m_vmax,m_vmin)) 

//r=m_vmax-m_vmin;
// t=m_vmin+i*r/ygrid
//vmin*ygrid=i*r
// i=-min/r*ygrid
//if ( 
// f=m'/r' 
// f*(M-m)=m; f*M-m(f+1)=0i; n=f/(f+1)*M
}

} // pass

}

void  exp_range( D & x, const bool top, const D&  e )
{
const bool f=(top&&(x>0))||!(top&&(x<0));
D y=x*((f)?e:(1.0/e));
Ss ss;
ss<<std::setprecision(2)<<y;
ss>>x;
}

D  trunc( const D & x, const IdxTy n  )
{
Ss ss;
ss<<std::setprecision(n)<<x;
D y;
ss>>y;
return y;
}





void write_svg_grid(std::ostream & os, mjm_svg_writer& sw)
{


os<<sw.gtext_text(m_ylab, m_ylab_x, m_ylab_y, m_ylab_sz, m_colorleg,  StrTy("end") ,90);
//StrTy hstripes(const D & xz,const D & xf, const D & yz, const D & yf,const IdxTy & n, const StrTy& col, const D & thick)
os<<sw.hstripes( m_xz, m_xf, m_yz,  m_yf, m_ygrid, StrTy("black"), 1);
os<<sw.vstripes( m_xz, m_xf, m_yz,  m_yf, m_cats, StrTy("black"), 1);

}
void write_svg_ylabels(std::ostream & os, mjm_svg_writer& sw)
{

for(IdxTy i=0; i<=m_ygrid; ++i)
{
const D f=D(i)/m_ygrid;
const D y=m_tm+m_ysz*(1.0-f);
Ss ss;
ss<<std::setprecision(3);
D q=((m_vmin+(m_vmax-m_vmin)*f));
//if ((q>=1e3)&&(q<1e5)) ss<<std::fixed<<std::setprecition(0);
if ((q>=1e3)&&(q<1e5)) ss<<int(q);
else ss<<q;
const StrTy & text=ss.str();
// this ASSFUDD does not anchor on fudding right or fudding levt 
os<<sw.gtext_text(text, m_xleg, y, m_szleg, m_colorleg,  StrTy("end") ,0);
os<<CRLF;

}



}
void write_svg_xlabels(std::ostream & os, mjm_svg_writer& sw)
{

const IdxTy ll=m_cats+m_first_line;
for(IdxTy i=m_first_line; i<ll; ++i)
{
const IdxTy ip=i-m_first_line;
const IdxTy is=m_so[i];
//const Ragged::Line line=r.line(is);
//if (line.size()<=text_idx) continue; // fudd 
const StrTy & text= m_xlabels[is-m_first_line]; // line[text_idx];
const D x=m_lm+ip*m_catpitch+m_szleg;
os<<sw.gtext_text(text, x, m_yleg, m_szleg, m_colorleg,  StrTy("end") ,-45);
os<<CRLF;
} // i 
}
void write_svg_day_xlabels(std::ostream & os, mjm_svg_writer& sw)
{

const IdxTy ll=m_xlabels.size();
for(IdxTy i=0; i<ll; ++i)
{
const IdxTy ip=i; // -m_first_line;
const IdxTy is=m_so[i];
//const Ragged::Line line=r.line(is);
//if (line.size()<=text_idx) continue; // fudd 
const StrTy & text= m_xlabels[is]; // line[text_idx];
//MM_ERR(MMPR3(text,is,i))
const D x=m_lm+ip*m_catpitch+m_szleg;
os<<sw.gtext_text(text, x, m_yleg, m_szleg, m_colorleg,  StrTy("end") ,-45);
os<<CRLF;
} // i 
// TODO WTF???
//os<<sw.gtext_text(m_dates_label,m_lm+m_catpitch*(m_cats>>1) , m_yleg+10*m_szleg, m_szleg, m_colorleg,  StrTy("center") ,0);
MM_ERR(" ASSFUDD " <<MMPR3(m_y_date_loc,m_ys,m_dates_label))
os<<sw.gtext_text(m_dates_label,m_lm+m_catpitch*(m_cats>>1) , m_y_date_loc, m_szleg, m_colorleg,  StrTy("center") ,0);

os<<CRLF;

}


void write_svg_legend(std::ostream & os, mjm_svg_writer& sw, const IdxTy flags)
{
const D xnames=m_lm+10;
D y=m_ynames;
D xstarts=m_xstarts; // xnames-m_lm*.5; // ; y+=yinc; // m_namesz+2;
D xfn=m_xf+m_xend; // m_rm*.5;
D yinc=m_namesz+2;
D x=xstarts; // xnames-m_lm*.5;
D wr=.5*m_wr;
D pn=0;
const bool use_jinv=(flags==1);
const IdxTy njinv=m_jinv.size();
const IdxTy nlabels=(!use_jinv)?m_jmax:njinv;
for (IdxTy j=0; j<m_jmax; ++j)
{
const StrTy & text=m_names[j];
 const D p=10+text.length()*m_namesz;
if  ( p>pn) pn=p;
}
if (y==0)
{
if (pn==0) { y=.6*m_tm;}
else
{
IdxTy dy=IdxTy(.5+D(wr+pn)/D(xfn-xnames)*nlabels);
y=m_tm-(m_namesz+2)*dy;
if (y>=m_tm) y=.9*m_tm;
if (y<0) y=.1*m_tm;
}

}

if (!use_jinv) {
for (IdxTy j=0; j<m_jmax; ++j)
{
if (x>=xfn) { x=xstarts; y+=yinc; }
const StrTy & cj=m_colors[j];
const StrTy & text=m_names[j];
os<<sw.shade_rect_text(x,y,wr,m_hr,cj,.8);
os<<sw.gtext_text(text, x+wr, y+.5*m_hr , m_namesz, cj,  StrTy("start") ,0);
os<<CRLF;

if ( pn==0)  x+=wr+10+text.length()*m_namesz;
else  x+=wr+pn; // 10+text.length()*m_namesz;
} // j

} else // jinv
{
MM_LOOP(jj,m_jinv)
{
const IdxTy j=(*jj).second.j();
if (x>=xfn) { x=xstarts; y+=yinc; }
const StrTy & cj=m_colors[j];
const StrTy & text=m_names[j];
os<<sw.shade_rect_text(x,y,wr,m_hr,cj,.8);
os<<sw.gtext_text(text, x+wr, y+.5*m_hr , m_namesz, cj,  StrTy("start") ,0);
os<<CRLF;

if ( pn==0)  x+=wr+10+text.length()*m_namesz;
else  x+=wr+pn; // 10+text.length()*m_namesz;
} // j




} // jinv



}


//TODO FIXME want to force a grid line at zero sometimes 
void adjust_grid()
{
const D vmaxi=m_vmax;
const D vmini=m_vmin;
const IdxTy ygridi=m_ygrid;
D e= 1.1;  
exp_range(m_vmax,true,e);
exp_range(m_vmin,!true,e);
D del=1;
IdxTy ibest=ygridi;
for(IdxTy i=5; i<16; ++i)
{

D diff=m_vmax-m_vmin;
D idelg=D(i-0)/diff;
idelg=idelg-IdxTy(idelg);
if (idelg<del) {del=idelg; ibest=i; }
}
m_ygrid=ibest;
MM_ERR(MMPR4(m_vmax,vmaxi,m_vmin,vmini)<<MMPR2(m_ygrid,ygridi))
}

void write_day_svg(std::ostream & os, mjm_svg_writer& sw)
{
MM_ERR("write_svg")
os<<sw.start_text(" test foo",m_xs,m_ys);
os<<sw.frame_text("#FFFFFF",m_xs,m_ys);
os<<CRLF;
// need this to get the max and min ut shold draw laer 
//write_svg_day_data(os, sw,2);
MM_ERR(" writing grid")
adjust_grid();
write_svg_grid(os,sw);
MM_ERR(" writing xlabels")
write_svg_day_xlabels(os,sw);
MM_ERR(" writing ylabels")
write_svg_ylabels(os,sw);
write_svg_legend(os,sw,0);
MM_ERR(" writing data")
write_svg_day_data(os, sw,5);

MM_ERR(" writing done")

os<<sw.end_text();
os<<CRLF;


} // write_day_svg


void write_svg(std::ostream & os, mjm_svg_writer& sw,const IdxTy flags=0)
{
MM_ERR("write_svg")
os<<sw.start_text(" test foo",m_xs,m_ys);
os<<sw.frame_text("#FFFFFF",m_xs,m_ys);
os<<CRLF;
// need this to get the max and min ut shold draw laer 
const bool write_rect=((flags&(1<<0))==0);
const bool write_density=!write_rect;
const bool write_names_on_plot=((flags&(1<<1))!=0);
const bool adjust=true; 
if  ( write_rect) write_svg_data(os, sw,2);
if (write_density) write_svg_data(os, sw,8+(write_names_on_plot?16:0));
if (adjust) adjust_grid();
write_svg_grid(os,sw);
write_svg_xlabels(os,sw);
write_svg_ylabels(os,sw);
if ( write_rect) write_svg_legend(os,sw,0);
if ( write_density) write_svg_legend(os,sw,1);
if (write_rect) write_svg_data(os, sw,5);
//if (write_density) write_svg_data(os, sw,4);


os<<sw.end_text();
os<<CRLF;


} // write_svg

protected:

ReadWriteMap   m_rwm;
//Ragged m_r;

IdxTy m_text_idx,m_jmax,m_jfirst,m_min_line,m_cats;
IdxTy m_tm,m_bm, m_lm,m_rm;
D m_catpitch, m_ysz,m_ynames,m_namesz,m_yleg,m_szleg;
D m_xstarts,m_xend;
StrTy m_colorleg;
D m_xs,m_ys,m_wr,m_hr,m_xmin,m_vmin,m_vmax;
D m_y_date_loc;
IdxTy m_ygrid, m_xz, m_xf,m_yz,m_yf,m_xleg;
StrTy m_ylab;
D m_ylab_x, m_ylab_y, m_ylab_sz;
Strings m_colors, m_names,m_xlabels;

Numbers m_merit, m_avg;
Numbers m_x_frac, m_x_start,m_x_height;
IdxTy m_first_line;

MyBlock m_dense_data;
typedef std::map< IdxTy , std::map< IdxTy, D > > SparseData;
SparseData  m_sparse_data;

std::map<StrTy, MyBlock> m_dense_map;
std::map<StrTy, SparseData> m_sparse_map;


Ragged::sort_order m_so;

JinvMap m_jinv;

// other plot thypes
IdxTy m_bin_idx,m_serial_idx, m_date_idx, m_sample_idx, m_food_idx,m_value_idx;
StrTy m_dates_label;
//IdxTy m_samples;
D m_y_min, m_y_max,m_plot_thick;
IdxTy m_day_min,m_day_max;
std::map<StrTy,IdxTy> m_bins;
std::map<IdxTy,StrTy> m_serial_to_date;
std::map<StrTy,IdxTy> m_date_to_serial;
std::map<StrTy,IdxTy> m_samples;
std::map<StrTy,IdxTy> m_sample_order;
std::map<IdxTy,IdxTy> m_point_order;
std::map<StrTy,IdxTy> m_foods_to_plot;
std::map<StrTy,IdxTy> m_dogs_to_plot;
};  //layout_blocks







class latex_blocks 
{
typedef  canned_methods_traits::Tr  Tr;
//typedef  linc_graph_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef latex_blocks Myt;
//typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef mjm_ragged_table Ragged;
//typedef Tr::Ordering Ordering;
public:
static bool bb(const IdxTy bit, const IdxTy flags)
{
return ((1<<bit)&flags)==0;
}
template<class Ty>
static void table_header( OsTy & os, const IdxTy nval, const Ty & val_names, const IdxTy flags, const StrTy nm="Name")
{
if (bb(0,flags)) { os<<"\\begin{table}[H]"<<CRLF; os<<"\\centering"<<CRLF; }
if (!bb(3,flags)) { os<<"\\begin{minipage}{.45\\linewidth}"<<CRLF; }

if (bb(1,flags)) { 
os<<"\\begin{tabular}{|"; os<<"l|"; for(IdxTy p=1; p<=(nval); ++p) os<<"r|";
//for(IdxTy p=1; p<(nval); ++p) os<<"r|";
os<<"}"<<CRLF; }
if (bb(2,flags)) { 
os<<"\\hline"<<CRLF;
//os<<"Name";
if ( bb(4,flags)) { os<<nm;  
for(IdxTy p=0; p<(nval); ++p) os<<"&"<<val_names[p]; }
else
{
if ( nval!=0) os<<val_names[0];
for(IdxTy p=1; p<=(nval); ++p) os<<"&"<<val_names[p]; }

//}
os<<"\\\\"<<CRLF;
os<<"\\hline"<<CRLF;
}

} 

template<class Tx, class Ty>
static void table_line( OsTy & os, const IdxTy nval, const Tx & name, const Ty & v, const IdxTy prec=3)
{

os<<name;
for(IdxTy q=0; q<nval; ++q)
{
Ss ss;
ss<<std::setprecision(prec)<<v[q];
os<<"&"<<ss.str(); // v[q];
}
os<<"\\\\"<<CRLF;
}

template< class Ty>
static void table_line( OsTy & os, const Ty & v, const IdxTy prec=3)
{
const IdxTy nval=v.size();
for(IdxTy q=0; q<nval; ++q)
{
Ss ss;
ss<<std::setprecision(prec)<<v[q];
if (q!=0) os<<"&";
os<<ss.str(); // v[q];
}
os<<"\\\\"<<CRLF;
}
static void nonzed(OsTy &os,const char * nm, const StrTy & v)
{ nonzed(os,StrTy(nm),v); } 
static void nonzed(OsTy &os,const StrTy& nm, const StrTy & v)
{
if (v.c_str()[0]!=0) os<<"\\"+nm+"{"+v+"}"<<CRLF;
}
static StrTy segged(const char * base, const IdxTy part, const IdxTy parts,
const IdxTy flags)
{
return segged(StrTy(base),part,parts,flags);
}
static StrTy segged(const StrTy& base, const IdxTy part, const IdxTy parts,
const IdxTy flags)
{
const bool human=((flags&1)==1);
if (parts<2) return base;
if (base.c_str()[0]==0) return base;
Ss ss;
if (human) { ss<<"Part "<<part<<" of "<<parts<<". "<<base; }
else { ss<<base<<part; }
return ss.str();
}

static void table_end( OsTy & os, const IdxTy flags,
     const StrTy caption="",const StrTy & label="" )
{
os<<"\\hline"<<CRLF;
os<<"\\end{tabular}"<<CRLF;
if (! bb(3,flags) ) { os<<"\\end{minipage}"<<CRLF;} 
nonzed(os,"caption",caption);
nonzed(os,"label",label);
os<<"\\end{table}"<<CRLF;
}

template <class Ty> 
static void add_superscripts(StrTy & s, const Ty v, const IdxTy flags=0)
{
const bool always=(flags&1);
const IdxTy sz=v.size();
if (!always) if (sz==0) return;
//s="$\\textrm{"+s+StrTy("}^{\\left(");
if (sz) s="$\\textrm{"+s+StrTy("}^{\\left(");
else s="$\\textrm{"+s+StrTy("}");
if (sz) s=s+v[0];
for (IdxTy i=1; i<sz; ++i) s=s+StrTy(",")+v[i];
if (sz) s=s+StrTy("\\right)}$");
else  s=s+StrTy("$");

}
static void variable(OsTy & os, const StrTy & nm, const StrTy & v)
{ os<<"\\newcommand{\\"<<nm<<"}{"<<v<<"}"<<CRLF; }

static void comment(OsTy & os, StrTy & text, const bool block)
{
const StrTy fudd;
comment(os,text,block,fudd); 
}
static void comment(OsTy & os, StrTy & text, const bool block, const StrTy & ty
)
{ 
if (!block) { os<<"%"<<text<<CRLF; return ; }

os<<"\\begin{"<<ty<<"comment}"<<CRLF;
os<<text<<CRLF;
os<<"\\end{"<<ty<<"comment}"<<CRLF;

}


// add capability to make table more rectangular for
// essentially lists - side by side etc 
template <class Ty>
static void ragged_table(OsTy & os, Ragged & r, Ty & m,const IdxTy inflags=0)
{
const IdxTy rsz=r.size();
const IdxTy headers=Tr::myatoi(m["headers"]);
IdxTy maxsz=Tr::myatoi(m["maxsz"]);
IdxTy minsz=Tr::myatoi(m["minsz"]);
IdxTy side_by_side=Tr::myatoi(m["side-by-side"]);
const StrTy caption=m["caption"];
const StrTy label=m["label"];
const bool sbs=(side_by_side>0);
const IdxTy hshift=((inflags&1)!=0)?16:0;
Ss hs,hs2,os2;
const Ragged::Line & hline=(rsz>=headers)?r.line(0):Ragged::Line();
//const IdxTy hsz=hline.size();
const IdxTy hsz=hline.size()-1;
if (hsz== (~0)) {
MM_ERR(" hline size is zero ")
return; 
}
const IdxTy flags=sbs?8:0;
table_header(hs,hsz,hline,flags+hshift);
os<<hs.str();
if (sbs)
{
table_header(hs2,hsz,hline,9+hshift);
}

StrTy cap=caption;
StrTy lbl=label;
IdxTy part=0;
IdxTy crsz=sbs?( rsz/2):rsz;
IdxTy parts=(maxsz==0)?1:(crsz+maxsz-1-headers)/maxsz;
IdxTy j=0;
for(IdxTy i=headers; i<rsz; ++i)
{
if ((maxsz!=0)&&((j%maxsz)==0)&&(i!=headers))
{++part;
cap=segged(caption,part,parts,1);
lbl=segged(label,part,parts,0);
if ( sbs ) { 
os<<"\\end{tabular} \\end{minipage}"<<CRLF;
os<<hs2.str();
os<<os2.str();
MM_ERR(MMPR(hs2.str()))
hs2.str("");
os2.str("");
}
table_end(os,flags, cap,lbl); os<<hs.str();

} // maxsz

table_line(os,r.line(i),0);
++j;
if ( sbs ){ ++i;  if ( i<rsz){   table_line(os2,r.line(i),0);} } 
}

if ( parts>1) 
{ ++part;
cap=segged(caption,part,parts,1);
lbl=segged(label,part,parts,0);
}
if ( sbs ) { 
os<<"\\hline"<<CRLF;
os<<"\\end{tabular}"<<CRLF;
os<<"\\end{minipage}"<<CRLF;
os<<hs2.str();
os<<os2.str();
MM_ERR(MMPR(hs2.str()))
}
table_end(os,flags,cap,lbl);
} // ragged_table




}; // lastex_blocks

// needs a param glod  type 
#if 0 
class mjm_digraph
{
//typedef  violin_box_etc_traits::Tr  Tr;
typedef  canned_methods_traits::Tr  Tr;
//typedef  linc_graph_traits::util  Util;
//typedef  Traits::Tr  Tr;
typedef mjm_digraph Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
//typedef Tr::Ordering Ordering;
typedef std::vector<StrTy> Words;

typedef Tr::MyBlock  MyBlock;
typedef violin_box_etc_params ParamGlob;

typedef string_tokenizer St;
typedef std::map<StrTy,StrTy> SeqMap;
class node_entry
{
public:
typedef node_entry Myt;
node_entry(): m_v(0) {}
node_entry(const D & x): m_v(x) {}
operator double () const   { return m_v; } 
Myt & operator+=( const D & x) { m_v+=x; return *this; } 
D m_v;
}; // node_entry

typedef IdxTy Nono;
typedef node_entry Novo;
typedef std::map<Nono,Novo> Dmap;
typedef std::map<Nono,Dmap> Gr;

typedef std::map<Nono,IdxTy> Hm;
typedef std::map<Nono,D> Score;

//typedef std::map<IdxTy,IdxTy> TermMap;
//typedef std::map<IdxTy,TermMap> SpeciesMap;
//typedef std::map<IdxTy,SpeciesMap> GenusMap;

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); } 
int myatoi(const char * c) const { return ::strtol(c,0,0); }

public:
mjm_digraph(): m_read_only(false)
,m_min_words(bad()),m_max_words(0),m_fields(0)

 {}
static const IdxTy &  bad() { static const IdxTy b=~0U;  return b; }
const IdxTy  size() const { return m_gr.size(); }

void load(const StrTy & fn )
{
std::ifstream  ifs(fn);  
load(ifs); 
}
void load(std::istream & is )
{
 IdxTy lines=0;
CommandInterpretter li(&is);
Setup(li,m_flp);
    while (li.nextok())
    {
        const IdxTy sz=li.size();
		if (sz<3) 
		{
			MM_ERR(" bad line "<<li.line())
			continue;
		}
		const Words & w=li.words();
		const D score=atof(w[0].c_str());
		const IdxTy dest=m_st(w[1]);
		const IdxTy src=m_st(w[2]);
		m_gr[src][dest]+=score;
	//	m_table.add(li.words());
		++lines;
    } // nextok()
//MM_ERR(" done loading "<<MMPR3(m_min_words,m_max_words,m_fields)<< MMPR3(lines,size(),m_order.dump(0)))
//MM_ERR(MMPR(lines)<< stats())
MM_ERR(MMPR(lines))
}
void subs(const IdxTy flags) 
{
const bool dump_walk=((flags&(1<<0))==0);
const bool max2=((flags&(1<<1))!=0);
//IdxTy ml=0;
//IdxTy lev_lim=10;
//const IdxTy max_level=ml;
MM_LOOP(ii,m_gr)
{
Hm  hits;
Hm  level;
Score score;
const IdxTy & i=(*ii).first;
const StrTy & nm=m_st(i);
++hits[i];
const Dmap & dm= (*ii).second;
walk(m_gr,dm,hits,level,1,score,1);
if (dump_walk) { 
IdxTy max=0;
if (max2)
{
MM_LOOP(jj, hits) {
	IdxTy j=(*jj).first; 
	const IdxTy & lvl=level[j];
if ( lvl>max) max=lvl;
//if (max>lev_lim) break; 	
 } // j 
//if (max>lev_lim) continue; 	
} // max2

MM_LOOP(jj, hits) 
{	
	IdxTy j=(*jj).first; 
	const StrTy & nmj=m_st(j); 
	const IdxTy & lvl=level[j];
	const D  & s=score[j];
//	if ((max_level==0)||(lvl<max_level)) 
	{ MM_ERR(MMPR4(nm, lvl,nmj,s)) }

} 
}

} // ii 


} // subs

void walk(const Gr & gr, const Dmap & dm, Hm & hits, Hm & level, const IdxTy d, Score & s,  const D & scorein )
{

MM_LOOP(ii, dm)
{
const IdxTy & i=(*ii).first;
//if (hits[i]!=0) continue; // return; 
// this needs to get to the end... 
//if (hits.find(i)!=hits.end()) continue; // return; 
const StrTy & nm=m_st(i);
if (false) {MM_ERR(MMPR(nm)) } 
++hits[i];
IdxTy&  dh=level[i];
const D scoreout=scorein*(*ii).second;
if (( dh==0)||( dh>d)) {  dh=d;  s[i]=scoreout;  } 
if (hits[i]>1) continue; // return; 
const auto ig= gr.find(i);
const Dmap & dmi= (*ig).second;
walk(gr,dmi,hits,level,d+1,s,scoreout);
} // ii 

} // walk




private:
void Setup( CommandInterpretter&  li, const ParamGlob & flp)
{
//li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,';');
li.readline_ok(false); li.use_alt_eol('\r',false); li.set_split(1,' ');
}

ParamGlob m_flp;
//SeqMap m_map;
Gr m_gr;
//Ordering m_order;
St m_st;
bool m_read_only;
IdxTy m_min_words,m_max_words,m_fields;

}; // mjm_digraph end

#endif






class mjm_canned_methods 
{
typedef canned_methods_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_canned_methods Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;


//typedef linc_graph_params ParamGlob;
typedef mjm_logic_base VariableStore;


typedef mjm_ragged_table Ragged;
typedef std::map<StrTy, Ragged> RaggedMap;

typedef mjm_fasta::fasta_file Fasta;
typedef std::map<StrTy, Fasta> FastaMap;

//typedef mjm_clustalw_aln Aln;
//typedef std::map<StrTy, Aln> AlnMap;

//typedef mjm_digraph Dig;
//typedef std::map<StrTy, Dig> DigMap;


//typedef mjm_char_mat Vec;
typedef mjm_char_mat CharMat;


//typedef std::vector<StrTy> Words;
typedef Ragged::Line Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
typedef mjm_tax_tree TaxTree;
typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;

#ifdef USE_RATIONALS
typedef mjm_rational_things  MulRat;
typedef  MulRat::big_coef_type big_coef_type;
typedef mjm_table_base<big_coef_type> TableMan;
#endif



public:

static int myatoi(const StrTy & s ) { return myatoi(s.c_str()); } 
static int myatoi(const char * c) { return ::strtol(c,0,0); }
static int myatoix(const char * c) { 
if (*c=='#') 
{
return myatoi(StrTy("0x")+StrTy(c+1)); 
}
return ::strtol(c,0,0); 
}


public :
//mjm_linc_graph():m_dmel(new Dmel()) {Init();}
#if 0
mjm_linc_graph(int argc,char **_args) : m_dmel(new Dmel())
{
// not sure this should be done first, user can invoke it 
Init();
// kluge to allow defaults lol 
const IdxTy ikluge=argc+10;
char * args[ikluge];
char dummy[2]; dummy[0]=0;
for (IdxTy i=0; i<IdxTy(argc); ++i) args[i]=_args[i];
for (IdxTy i=argc; i<ikluge; ++i) args[i]=&dummy[0];
int i=1;
// yeah probably a map or something better but this is easy 
while (i<argc)
{
const int istart=i;
//m_tree.config("-tree",i,argc,args);
//m_flp.config("-params",i,argc,args);
//configi(m_points,"-points",i,argc,args);
//m_flp.config_set("-set-param",  i,  argc, args);
//m_tree.config_set("-set-branch",  i,  argc, args);
cmdlcmd( i, argc, args);
if (i==istart) {++i; MM_ERR(" did nothing with "<<args[i]) } 

}
}
~mjm_linc_graph()
{
//clear_handlers();
delete m_dmel;
}
#endif

////////////////////////////////////////////////////////
// command block


static void stream_edit_fasta(const StrTy & d,const StrTy & s,const IdxTy flags)
{
std::ifstream is(s);
std::ofstream os(d);
MM_ERR(" stream_edit_fasta "<<MMPR2(d,s))
std::map<StrTy,IdxTy> komap;
std::map<StrTy,IdxTy> igmap;
igmap[StrTy("Uncultured")]=1;
IdxTy entries=0;
IdxTy wrote=0;
bool writing=true;
CommandInterpretter li(&is);
//ConmfigureReader(  li,debug);
    while (li.nextok())
    {
    // TODO FIXME ASSFUDD putting ok() in nextok() IS FUDDED ASDFASDFASD
        if (!li.ok()) break;
        const IdxTy sz=li.size();
		if (sz==0) continue;
		const char c=li.line()[0];
		if (c=='>')
		{
			++entries;
			if (sz<2) break; // just keep whatever
			IdxTy kloc=1;
			StrTy key=li.words()[kloc];
			while (igmap.find(key)!=igmap.end()) 
			{
				++kloc;
			 	if ( kloc>=sz) break;
				key=li.words()[kloc];
			}
			++komap[key]; 
			writing= (komap[key]==1);
			if (writing) ++wrote;
		}
		if (writing) { os<<li.line()<<CRLF;  }
	} // nextok()

MM_ERR(MMPR2(entries,wrote))
}


static void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv )
{
const StrTy d=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
stream_edit_fasta(d,s,flags);
MM_ERR("stream edit fasta fn into name "<<MMPR3(d,s,flags))
}



static void cmd_read_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
{
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
m[name].load(fn);
MM_ERR("reading fasta fn into name "<<MMPR3(fn,m[name].size(),name))
}



static void cmd_write_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
{
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
MM_ERR("writing fasta fn into name "<<MMPR3(fn,m[name].size(),name))
Fasta & fasta=m[name];
std::ofstream os(fn);
MM_SZ_LOOP(i,fasta,fsz)
{
os<<fasta.name(i)<<CRLF;
os<<fasta.seq(i)<<CRLF;
}
}


static void cmd_add_to_fasta(Cip & cip , LocalVar & lv, FastaMap & m  )
{
const StrTy name=cip.p1;
const StrTy seq=cip.p2;
Fasta & f=m[name];
StrTy sn=cip.wif(3);
if (sn.length()==0)
{ Ss ss ; ss<<"seq"<<(f.size()+1); sn=ss.str(); }
sn=StrTy(">") + sn;
f.add(sn,seq); // load(fn);
MM_ERR("adding to  "<<MMPR4(name,m[name].size(),sn,seq))
}

static void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv, FastaMap & m, RaggedMap & rm  )
{
//const StrTy name=cip.p1;
const StrTy dfasta=cip.p1;
// file with sequence numbers only 
const StrTy sfasta=cip.p2; // wif(3);
// annotation in zymo format
const StrTy srag=cip.wif(3);
const StrTy qiime=cip.wif(4);
//CommandInterpretter li;
Fasta & df=m[dfasta];
const Fasta & sf=m[sfasta];
const Ragged & sr=rm[srag];
const Ragged & qiiq=rm[qiime];
// xlate_map xlate_field_map(const IdxTy key, const IdxTy val)
Ragged::xlate_double_map qii_map=qiiq.xlate_field_max_map(8,3);
MM_ERR(" SHTFUDD "<<MMPR(qii_map.size()))
std::map<StrTy, StrTy> kv;
const IdxTy szsr=sr.size();
for (IdxTy i=0; i<szsr; ++i)
{
//MM_ERR(" ASSFUDD "<<MMPR2(i,sr.line(i).size()))
const IdxTy fudd=sr.line(i).size();
if (fudd!=1)
{
MM_ERR(" the ragged file needs to be loaded with flag bit 0 set (1) "<<MMPR3(i,fudd,srag))
return;
}
const char * c=sr.line(i)[0].c_str();
StrTy k="";
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; k=StrTy(d);
while (c[j]!=0) { if ((c[j]!='\t') &&(c[j]!=' ')) break; ++j; }
kv[k]=StrTy(c+j);

}

MM_ERR("zumo merge fasta fn into name "<<MMPR3(dfasta,sfasta,srag)<<MMPR3(df.size(),sf.size(),sr.size()))
const IdxTy sz=sf.size();
StrTy sep=" ";
for(IdxTy i=0; i<sz; ++i)
{
const StrTy & seq=sf.seq(i);
StrTy sn=sf.name(i);
StrTy sqname="";
const char * c=sn.c_str();
if (c[0]!=0) { 
IdxTy j=0;
while (c[j]!=0) { if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
char d[j+1]; memcpy(d,c,j); d[j]=0; 
StrTy k=StrTy(d+1);

StrTy annot=kv[k];
const char * shotfudd=sn.c_str();
sqname=StrTy(shotfudd+1); // sn; // k;
//MM_ERR(MMPR2(k,annot))
StrTy sfx=""; 
auto im=qii_map.find(sqname);
//MM_ERR(" looking up "<<MMPR(sqname))
if (im!=qii_map.end()) { Ss sfudd; sfudd<<" qzscore="; sfudd<<(*im).second; sfx=sfudd.str(); } 
else { 
// MM_ERR(" ASSFUCL "<<MMPR2(sqname,(*(qii_map.begin())).first))

 } 
sn=sn+sep+zymo_condense(annot)+sfx;
}
df.add(sn,seq); // load(fn);
}



}


static StrTy zymo_condense(const StrTy & s ) // const
{
const char * c=s.c_str();
const IdxTy sz=s.length();
//IdxTy start=0;
StrTy genus="";
StrTy species="";
if (sz<3) return StrTy();
for (IdxTy i=0; i<(sz-3); ++i)
{
if ((c[i]=='g')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
// fudding buserror? wtf memcpy fudded up? 
// the bus error got me, I think the length j is just wrong doh... 
char d[j+1-i]; //memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
genus=StrTy(d+3);
continue;
}
if ((c[i]=='s')&&(c[i+1]=='_')&&(c[i+2]=='_'))
{
IdxTy j=i+3;
while (c[j]!=0) {if (c[j]==';') break;  if (c[j]=='\t') break; if (c[j]==' ') break; ++j; }
// memcpy probably ok just had wrong length.. 
char d[j+1-i]; // memcpy(d,c+i,j); 
d[j-i]=0; 
for (IdxTy k=0; k<(j-i); ++k) d[k]=c[k+i];
species=StrTy(d+3);
continue;
}



} // i 
return genus+StrTy(" ")+species;
} 

//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;




template <class Trm> 
static void cmd_transpose_ragged(Cip & cip , LocalVar & lv, Trm & rm  ) 
//static void cmd_transpose_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
{
const StrTy cmd=cip.cmd();
const StrTy dname=cip.p1;
const StrTy sname=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
typename Trm::mapped_type & r=rm[sname];
//Ragged & r=rm[sname];
typename Trm::mapped_type  & d=rm[dname];
//Ragged & d=rm[dname];
MM_ERR(" cmd_transpose_ragged "<<MMPR3(dname,sname,flags)<<MMPR2(r.size(),d.size()))
r.transpose(d);
MM_ERR(" cmd_transpose_ragged "<<MMPR3(dname,sname,flags)<<MMPR2(r.size(),d.size()))
} // transpose

template <class Trm> 
static void cmd_read_ragged(Cip & cip , LocalVar & lv, Trm & rm  ) 
//static void cmd_read_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
{
const StrTy cmd=cip.cmd();
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const IdxTy nwords=myatoi(cip.wif(4));
const bool load_string=(cmd=="string-ragged");
const bool read_code=(cmd=="read-code");
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
const bool csv=((flags&16)!=0);
const bool pdftext=((flags&32)!=0);
const bool ignore_2delims=((flags&64)!=0);
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_read_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR4(debug,csv,pdftext,ignore_2delims)<<MMPR4(load_string, load_parsed,tabsep,ignore_hash))
typename Trm::mapped_type  & r=rm[name];
//Ragged & r=rm[name];
r.sep(" ");
if (tabsep) r.sep("\t");
if (csv) r.sep(",");
if (ignore_2delims) r.splitter(6);
if (ignore_hash) r.ignore_hash(ignore_hash);
if (pdftext) r.pdf_columns();
if (read_code) r.code_splitter(nwords);
if (load_string) r.load_from_string(fn,debug);
else 
if (load_parsed) r.load(fn,debug);
else r.load_lines(fn,debug);
MM_ERR(MMPR2(r.size(),name))
}

/*


StrTy dump( const IdxTy flags=0,const StrTy &  s=" ") const
{
Ss ss;
IdxTy i=0;
const bool use_space=((flags&1)!=0);
const bool add_seq=((flags&2)==0);
const bool add_quote=((flags&4)!=0);
const bool add_escapes=((flags&8)!=0);
const bool use_s=((flags&16)!=0);


*/

template <class Trm> 
static void cmd_write_ragged(Cip & cip , LocalVar & lv, Trm & rm  ) 
//static void cmd_write_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR( "cmd_write_ragged changed recently fudded up flags ")
/*const bool load_string=(cmd=="string-ragged");
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
const bool csv=((flags&16)!=0);
const IdxTy wflags=0;
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_write_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR2(debug,csv)<<MMPR4(load_string, load_parsed,tabsep,ignore_hash)<<MMPR(wflags))
//Ragged & r=rm[name];
*/
typename Trm::mapped_type  & r=rm[name];
/*
r.sep(" ");
if (tabsep) r.sep("\t");
if (csv) r.sep(",");
if (ignore_hash) r.ignore_hash(ignore_hash);
*/
//if (load_string) r.load_from_string(fn,debug);
//else 
//if (load_parsed) r.load(fn,debug);
//else r.load_lines(fn,debug);
const IdxTy wflags=flags;
r.write_file(fn,wflags);

MM_ERR(MMPR2(r.size(),name))
}


template <class Trm,class Os > 
static void cmd_dump_ragged(Os & os, Cip & cip , LocalVar & lv, Trm & rm  ) 
//static void cmd_write_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn="ostream"; // cip.p1;
const StrTy name=cip.p1;
const IdxTy wflags=myatoi(cip.wif(3));
const IdxTy flags=myatoi(cip.p2);
const bool load_string=(cmd=="string-ragged");
const bool load_parsed=((flags&1)==0);
const bool tabsep=((flags&2)!=0);
const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
const bool csv=((flags&16)!=0);
//const IdxTy wflags=0;
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_dump_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR2(debug,csv)<<MMPR4(load_string, load_parsed,tabsep,ignore_hash)<<MMPR(wflags))
//Ragged & r=rm[name];
typename Trm::mapped_type  & r=rm[name];
r.sep(" ");
if (tabsep) r.sep("\t");
if (csv) r.sep(",");
if (ignore_hash) r.ignore_hash(ignore_hash);
//if (load_string) r.load_from_string(fn,debug);
//else 
//if (load_parsed) r.load(fn,debug);
//else r.load_lines(fn,debug);
r.dump_os(os,wflags);

MM_ERR(MMPR2(r.size(),name))
}





template <class Trm> 
static void cmd_add_ragged(Cip & cip , LocalVar & lv, Trm & rm  ) 

//static void cmd_add_ragged(Cip & cip , LocalVar & lv, RaggedMap & rm  ) 
{
const StrTy cmd=cip.cmd();
const StrTy name=cip.p1;
const StrTy fn=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
//const bool load_string=(cmd=="string-ragged");
//const bool load_parsed=((flags&1)==0);
//const bool tabsep=((flags&2)!=0);
//const bool ignore_hash=((flags&4)!=0);
const bool debug=((flags&8)!=0);
//MM_ERR(" cmd_read_ragged from "<<fn<<" to "<<name<<" flags "<<flags)
MM_ERR(" cmd_add_ragged from "<<MMPR4(cmd,fn,name,flags)<<MMPR(debug))
//Ragged & r=rm[name];
typename Trm::mapped_type  & r=rm[name];
r.add(cip.words2());

MM_ERR(MMPR2(r.size(),name))
}

template <class Ttm > static void cmd_tt(Cip & cip , LocalVar & lv, Ttm & ttm ) 
{
const StrTy cmd=cip.p1;
const StrTy name=cip.p2;
const StrTy p1=cip.wif(3);
TaxTree & tt=ttm[name];
MM_ERR(" cmd_tt "<<MMPR3(cmd,name,p1))
do { 
if (cmd=="load-tax-tree") { tt.read_ncbi_dmp(p1);  continue; }
if (cmd=="dump-tax-tree") { tt.dump(p1);  continue; }
if (cmd=="write-tax-single") { tt.write_single(p1);  continue; }
if (cmd=="save-tax") { tt.write_composite(p1);  continue; }
if (cmd=="load-tax") { tt.read_composite(p1);  continue; }
if (cmd=="dump-lineages") { tt.dump_lineages(p1,2);  continue; }
if (cmd=="dump-normal-lineages") { tt.dump_lineages(p1,3);  continue; }
if (cmd=="check-normal-lineages") { tt.dump_lineages(p1,4);  continue; }
if (cmd=="traverse-full-tree") { tt.traverse_full_tree(p1,4);  continue; }
} while (false); 


}
// this needs a rational version 
static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////


}; //mjm_canned_methods



/////////////////////////////////////////////////////////

#ifdef  TEST_CANNED_METHODS__
int main(int argc,char **args)
{
typedef mjm_linc_graph  Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif

#undef MM_DMEL 

#endif

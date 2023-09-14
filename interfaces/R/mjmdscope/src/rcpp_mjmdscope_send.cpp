
#include <Rcpp.h>
//#include "../../mjm_dscope_interface.h"
#include "mjm_dscope_interface.h"

using namespace Rcpp;

typedef mjm_dscope_interface<> Dscope;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef mjm_string_base_params<Dscope::traits_type> BaseParams;
typedef double ValTy;
typedef mjm_block_matrix<ValTy> Block;

class _fix_globs
{
public:
_fix_globs() { mjm_global_flags::mm_err_enable=false; }
}; // _fix_globs
_fix_globs glob_fix;

Dscope m_dscope;
Ragged m_dscope_header;


void dscope_load_header(Ragged & r, const StrTy & id )
{
Dscope::traits_type::Ss ss;
//ss<<"# id test pattern chunks 0 "<<CRLF;
//ss<<"# oscilloscope 2D "<<CRLF;
ss<<"# id test pattern chunks "<<id<<CRLF;
ss<<"# time-snap 1D "<<CRLF;
ss<<"# trace "<<CRLF;
//ss<<"# trace "<<CRLF;
//ss<<""<<CRLF;
r.load(ss,false);
} // dscope_load_header

// [[Rcpp::export]]
List mjmdscopesend( SEXP ind ) {
 	Rcpp::DataFrame dfin(ind);
static bool init=false;
if (!init)
{
std::string  dscope_init="rawfifo launch  2"; // std::string(to<string*>(args[2]));
/*
Dscope::traits_type::Ss ss;
//ss<<"# id test pattern chunks 0 "<<CRLF;
//ss<<"# oscilloscope 2D "<<CRLF;
ss<<"# id test pattern chunks R "<<CRLF;
ss<<"# time-snap 1D "<<CRLF;
ss<<"# trace "<<CRLF;
//ss<<"# trace "<<CRLF;
//ss<<""<<CRLF;
m_dscope_header.load(ss,false);
*/
dscope_load_header(m_dscope_header,"R");
m_dscope.load(dscope_init,0);

 { mjm_global_flags::mm_err_enable=false; }
init=true;
} // init 
int len=dfin.length();   
int sz=dfin.nrows();
Block dblock(sz,len);
 for (int i=0;i<len;++i) {
// there is a fail_fast interface if there is no listenre... 
NumericVector fu=dfin[i]; 
 for (int j=0;j<sz;++j) {
// first is column
double v=fu[j];
// f-ing onje based f 
dblock(j,i)=v;
//MM_ERR(MMPR4(i,j,v,dblock(j,i)))
} // j 
} // i 
Ragged r =m_dscope_header;
r.new_line();
r<<"#";
r<<"names";
const CharacterVector sh=dfin.names();
 for (int i=0;i<len;++i) { r<<sh[i]; }
//r.new_line();
m_dscope.send(r,dblock,0);

   // CharacterVector x = CharacterVector::create( "foo", "bar" )  ;
   // NumericVector y   = NumericVector::create( 0.0, 1.0 ) ;
    List z  ; //           = List::create( x, y ) ;

    return z ;
}

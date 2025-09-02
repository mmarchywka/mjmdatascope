#ifndef MJM_GLOBALS_H__
#define MJM_GLOBALS_H__


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace mjm_global_flags
{

bool mm_err_enable=true; 
bool mm_delete_temps=true; 
int mm_err_io=0;
std::vector<int> State;
void push()
{
int s=mm_err_enable?1:0;
 s|=mm_delete_temps?2:0;
State.push_back(s);
} // push
void pop()
{
if (State.size())
{
int x=State.back();
mm_err_enable=(x&1);
mm_delete_temps=(x&2);
State.pop_back();
} // size

} // pop 

void push_verbose(const bool x )
{
push();
mm_err_enable=x; 
}


}; // mjm_global_flags

namespace mjm_global_credits
{
std::map<std::string, int> mm_credited;
std::vector<std::string> mm_credit_vector;
std::string credits()
{
std::string s;
auto ii=mm_credited.begin(); while (ii!=mm_credited.end())  
{ const auto & x=mm_credit_vector[(*ii).second]; if (x.length()>2) s=s+"\n"+x ; ++ii; }
return s;
}
void mm_credit(const std::string & src, const std::string & v)
{
if (mm_credited.find(src)==mm_credited.end()) { mm_credited[src]=mm_credit_vector.size();
mm_credit_vector.push_back(v); 
}
}
std::string about() 
{
std::string s;
for(unsigned int  i=0; i<mm_credit_vector.size(); ++i) s=s+mm_credit_vector[i]+"\n";

return s;
}
class credit
{
public:
credit(const std::string & src, const std::string & v) { mm_credit(src,v); } 
} ; // credit

}; // mjm_global_credits


#define MM_GETSET(t,x,y) const t & x()const  { return y; } void  x(const t & asdf )  { y=asdf; } 
#define MM_GETSETP(t,x,y)  t * x()  { return y; } void  x( t * asdf )  { y=asdf; } 

//#define SETGET(nm,mem,ty) void nm(const ty & x) { mem=x; } const ty & nm() const { return mem; }
#define MM_SETGETSTATE(nm,flag,bit) void nm(const bool x) { MM_ERR("ASSFUCK "<<x)  Set(flag,bit,x); MM_ERR("shita "<<flag ) } const bool nm() const { return Bit(flag,bit); }




#define MM_WAIT {std::cerr<<" wiat for input "<<__FILE__<<__LINE__<<"\n"; std::cerr.flush(); StrTy x;  std::cin >> x; }
#define MM_MARK __FILE__<<__LINE__
#define MM_MARKF __FUNCTION__<<__FILE__<<__LINE__
#define MM_TRACE std::cout<< __FILE__<<__LINE__<<"\n"; std::cout.flush();
#define MM_MSG(msg) {  std::cout<< __FILE__<<__LINE__<<" "<<msg<<"\n"; std::cout.flush(); } 
#define MM_MSG_CRLF {  std::cout<<"\n"; std::cout.flush(); } 
//#define MM_ERR(msg) std::cerr<< __FILE__<<__LINE__<<" "<<msg<<"\n"; std::cout.flush(); std::cerr.flush();

#define MM_TYPE typedef typename
#ifdef LINE_LOCK_MM_ERR
#warning line buffering stderr for mt code output do not hold locks in parameters 
// so where does the mutex live? these are all headers....
#include <pthread.h>

//void enter_serial() { pthread_mutex_lock( mutex1() ); }
//void exit_serial() { pthread_mutex_unlock( mutex1() ); }

///////////////////////////////////////////

#define MJM_EXCEPT_THROW(x) {Ss xxxss;xxxss<<MM_MARK<<" "<<x; throw std::logic_error(xxxss.str());}

#ifdef MJM_DIE_ON_DOH
#define MM_DOH(x) {Ss xxxss;xxxss<<MM_MARK<<" "<<x; throw std::logic_error(xxxss.str());}
#else
#define MM_DOH(x) {MM_ERR(x) }
#endif
/////////////////////////////////////////

// for only one compliation unit, otherwise ifdef the main doh 
namespace mjm_globals {

typedef pthread_mutex_t Mutex; 
Mutex* mutex() 
{
static  Mutex m=  PTHREAD_MUTEX_INITIALIZER;
return &m;
}
void global_io_lock_mutex() {  pthread_mutex_lock(mutex()); } 
void global_io_unlock_mutex() {  pthread_mutex_unlock(mutex()); } 
}; 

#define MM_ERR(msg) { { if ( mjm_global_flags::mm_err_enable) { std::cout.flush(); try { mjm_globals::global_io_lock_mutex(); std::cerr<< __FILE__<<__LINE__<<" "<<msg<<"\n";std::cerr.flush(); } catch (...) { mjm_globals::global_io_unlock_mutex(); throw ; } mjm_globals::global_io_unlock_mutex(); }}}
#define MM_ERRF(msg) { { if ( mjm_global_flags::mm_err_enable) { std::cout.flush(); try { mjm_globals::global_io_lock_mutex(); std::cerr<<__FUNCTION__<<":"<< __FILE__<<__LINE__<<" "<<msg<<"\n";std::cerr.flush(); } catch (...) { mjm_globals::global_io_unlock_mutex(); throw ; } mjm_globals::global_io_unlock_mutex(); }}}
#else
#ifdef MIXED_CERR_STREAMS
#define MM_ERR(msg) {{ if ( mjm_global_flags::mm_err_enable) { std::stringstream ____ss; ____ss<< __FILE__<<__LINE__<<" "<<msg<<"\n"; fprintf(stderr,"%s",____ss.str().c_str()); std::cerr.flush();  }}} 
#define MM_ERRF(msg) {{ if ( mjm_global_flags::mm_err_enable) { std::stringstream ____ss; ____ss<<__FUNCTION__<<":"<< __FILE__<<__LINE__<<" "<<msg<<"\n"; fprintf(stderr,"%s",____ss.str().c_str()); std::cerr.flush();  }}} 
#else
#define MM_ERR(msg) {if (mjm_global_flags::mm_err_enable) { std::cout.flush(); std::cerr<< __FILE__<<__LINE__<<" "<<msg<<"\n";std::cerr.flush(); }} 
#define MM_ERRF(msg) {if (mjm_global_flags::mm_err_enable) { std::cout.flush(); std::cerr<< __FILE__<<__LINE__<<" "<<__FUNCTION__<<" "<<msg<<"\n";std::cerr.flush();} } 
#endif //  MIXED_CERR_STREAMS

#endif // LINE_LOCK_MM_ERR

#define MM_SS(msg)  Ss ss;  ss<< __FILE__<<__LINE__<<" "<<msg;  
#define MM_SSF(msg)  Ss ss;  ss<< __FILE__<<__LINE__<<" "<<__FUNCTION__<<" "<<msg;  

#define MM_STATUS(msg) { std::cout.flush();  std::cerr<< __FILE__<<__LINE__<<" "<<msg<<"              \r"; std::cerr.flush(); } 
#define MM_ONCE(msg, action ) { static bool dun=false; if (!dun&&(mjm_global_flags::mm_err_enable)) {  std::cerr<< __FILE__<<__LINE__<<" ONCE "<<msg<<"\n"; std::cout.flush(); dun=true; }  action  } 
#define MM_DUNCE(msg, interval ) { static IdxTy i=0; ; if ((i%interval)==0) {  std::cerr<< __FILE__<<__LINE__<<" DUNCE "<<i<<" "<<msg<<"\n"; std::cout.flush(); std::cerr.flush();  }  ++i;  } 
// #define CRLF "\r\n"
// this doesn't work since getling put cr on last field.....
#define CRLF "\n"
#define CR "\r"
#define CRLFS 1 

//http://stackoverflow.com/questions/6671698/adding-quotes-to-argument-in-c-preprocessor

#define QFUDDBOOST(x) #x
#define QUOTE(x) QFUDDBOOST(x)
#define MM_STR_LOC StrTy( __FILE__)+StrTy(QUOTE(__LINE__))
#define MM_STR_TIMESTAMP StrTy( __DATE__)+StrTy(QUOTE(__TIME__))
//#include QUOTE(FILE_ARG)
// this should be exapanded in the text formatt classes lol 
//#define MMZ( val,ty) " "<<((val==ty())?"":QUOTE( val=))<<((val==ty())?"":val)
#define MMZS( val)  ((StrTy(val)==StrTy(""))?"":" ")<<((StrTy(val)==StrTy(""))?"":QUOTE( val=))<<((StrTy(val)==StrTy(""))?StrTy(""):StrTy(val))

std::string MJM_BLANK(double d) { std::stringstream ss; if (d!=0) ss<<d; return ss.str(); } 
#define MMN( val) " "<<((val==0)?"":QUOTE( val=))<<MJM_BLANK(val)
#define MMPR( val) " "<<QUOTE( val=)<<val
#define MMPRS( val) " "<<QUOTE( val )<<" "<<val
#define MMPRSL( val) " "<<QUOTE( val )<<" "<<val<<CRLF
#define MMPR2( val1,val2) " "<<QUOTE( val1=)<<val1<<" "<<QUOTE(val2=)<<val2
#define MMPR3( val1,val2,val3) " "<<QUOTE( val1=)<<val1<<" "<<QUOTE( val2=)<<val2<<" "<<QUOTE( val3=)<<val3
#define MMPR4( val1,val2,val3,val4) " "<<QUOTE( val1=)<<val1<<" "<<QUOTE( val2=)<<val2<<" "<<QUOTE( val3=)<<val3<<" "<<QUOTE( val4=)<<val4
#define MMPR4S( val1,val2,val3,val4) " "<<QUOTE( val1 )<<" "<<val1<<" "<<QUOTE( val2 )<<" "<<val2<<" "<<QUOTE( val3 )<<" "<<val3<<" "<<QUOTE( val4 )<<" "<<val4
#define MMPR4SL( val1,val2,val3,val4) " "<<QUOTE( val1 )<<" "<<val1<<CRLF<<QUOTE( val2 )<<" "<<val2<<CRLF<<QUOTE( val3 )<<" "<<val3<<CRLF<<QUOTE( val4 )<<" "<<val4<<CRLF

#define MMPR5( val1,val2,val3,val4,val5) " "<<QUOTE( val1=)<<val1<<" "<<QUOTE( val2=)<<val2<<" "<<QUOTE( val3=)<<val3 <<" "<<QUOTE( val4=)<<val4 <<" "<<QUOTE( val5=)<<val5

#define MM_FAIL {int * i=0; *i=0; } 
#define MM_FAULT {int * i=0; *i=0; } 
#define MM_DIE(val)  { { mjm_global_flags::push_verbose(true);  MM_ERR(val) ; }  int * i=0; *i=0; } 
#define MYASSERT( val) {if (!val) MM_ERR(MMPR((val)))    }  

#define MM_LOG(exp,label,msg) { if (exp) { MM_MSG(label<<" "<<msg) ; } }

// too confusing for syntax highlighter.. 
//#define MM_ILOOP(v,n)  for( IdxTy v=0; v<(n) ; ++v) { 
#define MM_ILOOP(v,n)  for( IdxTy v=0; v<(n) ; ++v)  
#define MM_LOOP(v,obj)  for( auto v=obj.begin(); v!=obj.end(); ++v)  
#define MM_SZ_LOOP(v,obj,sz)  const IdxTy sz= obj.size();  for( IdxTy v=0; v<sz;  ++v)  
#define MM_NCSZ_LOOP(v,obj,sz)  IdxTy sz= obj.size();  for( IdxTy v=0; v<sz;  ++v)  

//MM_FLAG(show_seg,0,1)
#define MM_FLAG(v,bit,p) const bool v=((flags&(1<<bit))==(p<<bit));


// yes, ths makes exe larger and make duplicates
// and confusion maybe but temp ontil settle on fianl solution

namespace {
typedef std::ostream OsTy;
typedef std::istream IsTy;
typedef std::string StrTy;


OsTy & osx() { return std::cout; }
OsTy & esx() { return std::cerr; }
//const bool PATRONIZING =true; 
const bool PATRONIZING =false; 
//const bool PEDANTIC_AND_CONDESCENDING =true; 
const bool PEDANTIC_AND_CONDESCENDING =false; 



}; // anon for flags egtc. 

class mjm_generic_traits
{

public:

typedef char ChTy;
typedef std::string StrTy;
typedef std::istream IsTy;
typedef std::ostream OsTy;
typedef std::ofstream Ofs;
typedef std::stringstream SsTy;
typedef std::stringstream Ss;
typedef unsigned int IdxTy;
typedef  int IntTy;
typedef  int DelTy;
typedef std::runtime_error ErTy;
typedef double D;
};


#endif // guard

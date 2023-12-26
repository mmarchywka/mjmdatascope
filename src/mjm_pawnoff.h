
#ifndef MJM_pawnoff_H__
#define MJM_pawnoff_H__
// moved from mikemail project 2019-09-13

#include <dirent.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <iostream>
#include <stdexcept>
#include <exception>
#include <typeinfo>
#include <stdio.h>
#include <cstdio>
#include <string>
#include <cstdio>


#include "mjm_globals.h"
#include <mjm_thread_util.h>
#include "mjm_data_model_error_log.h"
//#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_blob.h"
//#include "mjm_canned_methods.h"
//#include "mjm_cli_ui.h"
//#include "mjm_tokenized_collections.h"


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

System calls to handle things we don't do in code now. Mostly a crutch
for testing stubs but also useful for mail handling. 

*/



////////////////////////////////////////////////////////////////

class pawnoff_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D; 
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;
public:
pawnoff_params( const StrTy & nm) : Super(nm) {}
pawnoff_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
//IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
//IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
// samples, 
//D time_step() const { return m_map.get_double("time_step",1e-13); }
//int bulk_layers() const { return m_map.get_int("bulk_layers",6); }
//StrTy end_date() const { return m_map.get_string("start_date","9999-99-99"); }
//IdxTy accmode() const { return m_map.get_uint("accmode",0); } // // 100;
//bool print_counts() const { return m_map.get_bool("print_counts",!true); }
// should be geneated code, do not edit  to make every entry for one name  
StrTy to_string(const StrTy & sep=" ") const
{
Ss ss;
ss<<"log_commands="<<log_commands()<<sep;
ss<<"exit_on_err="<<exit_on_err()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
return ss.str();
}
}; // pawnoff_params





template <class Tr> 
class mjm_pawnoff
{
//typedef mikemail_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_pawnoff Myt;
//typedef mjm_cli_ui<Myt> CliTy;
MM_TYPE Tr::IdxTy IdxTy;
MM_TYPE Tr::D D;
MM_TYPE Tr::Ss Ss;
MM_TYPE Tr::IsTy IsTy;
MM_TYPE Tr::OsTy OsTy;
MM_TYPE Tr::Ofs Ofs;
//typedef Tr::MyBlock  MyBlock;
//typedef pthread_mutex_t Mutex;
//typedef Tr::Mutex Mutex;
// slow but it works 
// date  -d "Sat, 13 Apr 2019 08:01:28 +0000" "+%Y-%m-%d:%k:%M:%S"
// echo "Sat, 13 Apr 2019 08:01:28 +0000"| date -f -  "+%Y-%m-%d:%k:%M:%S"
// TODO FIXME verify that this can transfer buffers... 
//enum {  OUTBUFSZ=(1<<24), MAX_BASH_TEXT=(1<<12) }; 
// the cat version may be interretting dollar signs etc 
//enum {  OUTBUFSZ=(1<<10), MAX_BASH_TEXT=(1<<12) }; 
enum {  OUTBUFSZ=(1<<10), MAX_BASH_TEXT=(0) }; 
enum { STACK_LIMIT=1<<16,STRING_GUARD=2};


typedef mjm_blob<Tr> Blob;

typedef mjm_thread_util<Tr> ThreadUtil;
//typedef typename  mjm_thread_util<Tr>::mutex_vector MutexVector;
typedef typename  ThreadUtil::mutex_vector MutexVector;

class bg_param
{
public:
bg_param() : parent(0),mode(0) {}
bg_param(Myt * p, const StrTy & p1, const StrTy & p2, const IdxTy m) 
: parent(p),s1(p1),s2(p2),mode(m) {}
StrTy dump()
{
Ss ss;
ss<<MMPR4((parent==0), s1,s2,mode);
return ss.str(); 
}
Myt * parent;
StrTy s1,s2;
IdxTy mode;
};
typedef bg_param BgParam;

enum { MU_SERIAL, MU_WRITE,MU_SZ};

static MutexVector & mv()
{
static MutexVector m = MutexVector(MU_SZ);
return m;
}
static void SEnterSerial(const IdxTy i) { mv().enter_serial(i); }
static void SExitSerial(const IdxTy i)  { mv().exit_serial(i); }



static IdxTy  Serial( const IdxTy d)
{
SEnterSerial(MU_SERIAL);
static IdxTy s=0;
IdxTy x=s;
s+=d;
SExitSerial(MU_SERIAL);
return x;

}

bool Bit(const IdxTy f, const IdxTy b) { return ( (f&(1<<b))!=0); } 



// API

	public:
typedef Blob blob;
mjm_pawnoff() { Init(); }
~mjm_pawnoff() { Cleanup(); }
void clean() { Cleanup(); }
void verify_temps(const bool x) {  m_verify_temps_gone=x;}
void debug(const IdxTy n ) { m_debug=n; }  
IdxTy  debug()const  { return  m_debug; }  
IdxTy view(const StrTy & h)
{ return More(h); } 
IdxTy view(const StrTy & h, const StrTy & b) { return More(h); }
StrTy dates(const StrTy & dates) { return Dates(dates); }
IdxTy raw(const StrTy & s) { return Raw(s); }
//StrTy LexiDays(const StrTy & dates)
StrTy lexi_day() { return LexiDays(StrTy()); }
StrTy lexi_days(const StrTy & dates) { return LexiDays(dates); }
StrTy serial_days(const StrTy & dates) { return SerialDays(dates); }
StrTy date_default(const StrTy & dates) { return DatesDefault(dates); }
// echo $(( $( date -d `date -d "0 days ago "  +%Y-%m-%d` +%s) /86400 )) 
StrTy fgen(const StrTy & t=StrTy()) //const
{
char c[1+((L_tmpnam+1)|255)] ; 
const StrTy fin=tmpnam(c)+m_ext;
m_files.push_back(fin);
return fin;
} // fgen 

StrTy get_clipboard() // const
{
//onst StrTy x=m_test_hand.get_output("","myclip -echo");
// if this is not a url on the clipboard, it prints a bunch of junk to stderr
// needs to trp... a
return get_output("","myclip -echo");
}
// the decoded output should not be a string at least a C string ??? 
StrTy base64( const StrTy & d)
{
char * result;
IdxTy len;
Base64(&result,&len,d);
result[len]=0;
StrTy x= StrTy(result);
return x;
} // base64

StrTy html( const StrTy & d)
{
char * result;
IdxTy len;
Html(&result,&len,d);
result[len]=0;
StrTy x= StrTy(result);
// mjm 2019-09-13
delete [] result; 
return x;
} // html 
StrTy qp( const StrTy & d)
{
char * result;
IdxTy len;
Qp(&result,&len,d);
result[len]=0;
StrTy x= StrTy(result);
// mjm 2019-09-13
delete [] result; 
return x;
} // qp 
///////////////////////////////////////
// hard coded and hard to change as static with limite
// params, no In thing... 
IdxTy cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, const IdxTy flags) //  const 
{
const bool always=true; // Bit(flags,0);
const bool file= Bit(flags,1);
const bool notrailingcrlf= Bit(flags,2); // 4 
const bool dfile= Bit(flags,3); /// 8
Blob d,err,cout;
StrTy _cmd;
//if (file) { _cmd="cat \""+s+"\" |"+cmd; }
if (file) { _cmd=cmd; d.load(s); }
else { d=s; _cmd=cmd; }
//IdxTy rc=out.hand().fileio(cout,err,d,_cmd,3);
IdxTy rc=fileio(cout,err,d,_cmd,3);
if (always||(rc!=0)) { MM_ERR(MMPR3(dest,s,dfile)<< MMPR4(file,notrailingcrlf,flags,rc)<<MMPR4(cmd,StrTy(cout),StrTy(err),StrTy(d))) }
if (notrailingcrlf) {
IdxTy ncrlf= cout.no_trailing_crlf();
if (always) { MM_ERR(MMPR(ncrlf)) }
 }
if (!dfile) dest=StrTy(cout);
else { cout.save(dest); }
return rc;
}// cmd_exec
StrTy sexc(IdxTy & rc, const StrTy &cmd)
{
StrTy rv;
const IdxTy _rc=simple_cmd_exec(rv,StrTy(),cmd,4);
rc=_rc;
return rv;
}
IdxTy simple_cmd_exec(StrTy & dest, const StrTy & s, const StrTy &cmd, const IdxTy flags)
{
const bool always= Bit(flags,0);
const bool file= Bit(flags,1);
const bool notrailingcrlf= Bit(flags,2); // 4 
const bool dfile= Bit(flags,3); /// 8
const bool errnonzed= !Bit(flags,4); /// 16
Blob d,err,cout;
StrTy _cmd;
//if (file) { _cmd="cat \""+s+"\" |"+cmd; }
if (file) { _cmd=cmd; d.load(s); }
else { d=s; _cmd=cmd; }
IdxTy rc=fileio(cout,err,d,_cmd,3);
const bool wrong=(rc!=0)||(err.size());
if (always||(wrong)) { MM_ERR(MMPR3(dest,s,dfile)<< MMPR4(file,notrailingcrlf,flags,rc)<<MMPR4(cmd,StrTy(cout),StrTy(err),StrTy(d))) }
if (notrailingcrlf) {
IdxTy ncrlf= cout.no_trailing_crlf();
if (always) { MM_ERR(MMPR(ncrlf)) }
 }
const bool reterr=errnonzed&&wrong;
Blob & retblob=reterr?err:cout;
if (!dfile) dest=StrTy(retblob);
else { retblob.save(dest); }
return rc+err.size();
}// simple_cmd_exec







///////////////////////////////////////
StrTy today( const IdxTy flags=0) { return get_output("" , m_today); } 
StrTy now( const IdxTy flags=0) { return get_output("" , m_now_fn); } 
StrTy get_output( const StrTy & d, const StrTy & cmd)
{
char * result;
IdxTy len;
GetOutput(&result,&len,d,cmd);
result[len]=0;
StrTy x= StrTy(result);
// mjm 2019-09-13
delete [] result; 
return x;
} // get_output
StrTy get_output_with_input_bg( const StrTy & d, const StrTy & cmd)
{
// this needs to create a dummy instance that
// is only used by the thread as this is not thread safe
// or I guess it could be thread safe... 
//BgParam*  bg =  ( new BgParam(this,d,cmd,0)) ;
BgParam*  bg =  ( new BgParam(new Myt(),d,cmd,0)) ;
if (Bit(m_debug,0)) {  MM_ERR(" XXXX "<<MMPR(d)) } 
//bg->parent=this; bg->s1=d; bg->s2=cmd; bg->mode=0;
ThreadUtil::fire_and_forget(1,Myt::Bg, bg);
return StrTy();
}


IdxTy get_output_with_input(Blob & dest, const Blob & data, const StrTy & cmd)
{
char * result;
IdxTy len;
GetOutputWithInput(&result,&len,data,cmd);
//result[len]=0;
dest.adopt(result,len,len+1);
//StrTy x= StrTy(result);
// mjm 2019-09-13
//delete [] result; 
return 0;
} // get_output



StrTy get_output_with_input( const StrTy & d, const StrTy & cmd)
{
char * result;
IdxTy len;
if (Bit(m_debug,0)) { MM_ERR(" XXXX "<<MMPR2(d,cmd)) }
GetOutputWithInput(&result,&len,d,cmd);
result[len]=0;
StrTy x= StrTy(result);
// mjm 2019-09-13
delete [] result; 
return x;
} // get_output



// return stderr in return vlue 
template <class Ty> StrTy  get_output(Ty & buf,  const StrTy & d, const StrTy & cmd)
{
char * result;
IdxTy len;
GetOutput(&result,&len,d,cmd);
//result[len]=0;
buf.accept(result,len);
StrTy x= ""; // StrTy(result);
// mjm 2019-09-13
delete [] result; 
return x;
} // get_output






void base64( char ** result, IdxTy * len, const StrTy & d)
{ Base64(  result, len,  d); }

IdxTy cmd(const StrTy & c)
{
IdxTy rc=std::system(c.c_str());
if (rc!=0) MM_ERR(MMPR(rc))
return rc; 
}

IdxTy fileio(Blob & cout, Blob & cerr, const Blob & cin, const StrTy & cmd, const IdxTy flags=0)
{ return Fileio(cout,cerr,cin,cmd,flags); }

IdxTy lineio(Blob & cout, Blob & cerr, const Blob & cin, const StrTy & cmd)
{ return Lineio(cout,cerr,cin,cmd); }



	private:
IdxTy Fileio(Blob & cout, Blob & cerr, const Blob & cin, const StrTy & _cmd,const IdxTy flags)
{
const bool trail_nl=Bit(flags,0);
const bool add_parens=Bit(flags,1);
const bool no_ins=Bit(flags,2); // 4
const bool flag_debug=!true||Bit(flags,8);
// vi leaves in a trailing nl doh 
StrTy cmd=_cmd;
if (flag_debug){  MM_ERR(MMPR4(flags,trail_nl,add_parens,cmd))}
if (trail_nl||add_parens)
{
const IdxTy inlen=_cmd.length();
if (inlen ==0 ) return 0;
char c[inlen+10];
IdxTy base=add_parens?1:0;
if (add_parens) c[0]='(';
memcpy(c+base,_cmd.c_str(),inlen);
const IdxTy dun=int(inlen)-((trail_nl&&(c[base+inlen-1]=='\n'))?1:0)
+(add_parens?2:0); 
if (add_parens) c[dun-1]=')';
c[dun]=0;
cmd=StrTy(c);
}
MM_ONCE(" Fileio is not thread of process safe doh ", ) 
//char c[1+((L_tmpnam+1)|255)] ; 
const StrTy fin=fgen(); // tmpnam(c);
const StrTy fout=fgen(); // tmpnam(c);
const StrTy ferr=fgen(); // tmpnam(c);
const StrTy frc=fgen(); // tmpnam(c);
// copy cin to fin
cin.save(fin);
Blob blank,brc;
blank.save(fout);
blank.save(ferr);
Ss ss;
//ss<<"cat "<<fin<<" | " <<cmd <<" >"<<fout<<" 2>"<<ferr<<" ; echo $? >> "<<frc;
 if (no_ins) ss<<"echo -n | "<<cmd <<"" << " >"<<fout<<" 2>"<<ferr<<" ; echo $? >> "<<frc;
else  ss<<cmd <<" < "<< fin  << " >"<<fout<<" 2>"<<ferr<<" ; echo $? >> "<<frc;
//if (flag_debug) 
{MM_ERR("XXXX "<<ss.str()) }

IdxTy rc=std::system(ss.str().c_str());
// copy cout and cerr from their files 
cout.load(fout);
cerr.load(ferr);
brc.load(frc);
if (brc.ptr()) rc=atoi(brc.ptr());
else { rc=~0;  MM_ERR(" no return code for "<<MMPR(ss.str())) } 
return rc;

}

IdxTy Lineio(Blob & cout, Blob & cerr, const Blob & cin, const StrTy & cmd)
{
MM_ONCE(" Lineio is not thread of process safe doh ", ) 
//char c[1+((L_tmpnam+1)|255)] ; 
const StrTy fin=fgen(); // tmpnam(c);
const StrTy fout=fgen(); // tmpnam(c);
const StrTy ferr=fgen(); // tmpnam(c);
const StrTy frc=fgen(); // tmpnam(c);
// copy cin to fin
cin.save(fin);
Blob blank,brc;
blank.save(fout);
blank.save(ferr);
std::vector<StrTy> cmdv;
const char * sp=cmd.c_str();
const IdxTy len=cmd.length();
IdxTy start=0;
for(IdxTy i=0; i<len; ++i)
{
if (sp[i]=='\n'){
 if (i>=start)  cmdv.push_back(StrTy(cmd,start,i-start));
start=i+1;
} 
} // i 
 if (len>start)  cmdv.push_back(StrTy(cmd,start,len-start));
IdxTy rc=0;
for(IdxTy i=0; i<cmdv.size(); ++i)
{
Ss ss;
//ss<<"cat "<<fin<<" | " <<cmd <<" >"<<fout<<" 2>"<<ferr<<" ; echo $? >> "<<frc;
if (i==0) ss<<cmdv[i] <<" < "<< fin  << " >"<<fout<<" 2>"<<ferr<<" ; echo $? >> "<<frc;
else ss<<cmdv[i] <<" < "<< fin  << " >>"<<fout<<" 2>>"<<ferr<<" ; echo $? >> "<<frc;
rc|=std::system(ss.str().c_str());

} 

// copy cout and cerr from their files 
cout.load(fout);
cerr.load(ferr);
brc.load(frc);
rc=atoi(brc.ptr());
return rc;

}


std::vector<StrTy> m_files;
pawnoff_params m_params;
StrTy m_ext;
StrTy m_nl, m_bashsep,m_more,m_date,m_today,m_now_fn,m_lexi_day,m_date_default,m_base64,m_html,m_qp;
StrTy m_temp_base;
IdxTy m_width;
IdxTy m_debug;
void Cleanup()
{
const bool norem= (mjm_global_flags::mm_delete_temps );
if (m_verify_temps_gone) {
MM_ERR(" pawnoff cleanup "<<MMPR2(norem,m_files.size())) } 
if (!norem  ) return;

MM_LOOP(ii,m_files)
{
if ( m_verify_temps_gone ) {  MM_ERR(" removing "<<MMPR((*ii))) } 
std::remove((*ii).c_str());

} // files

} // Cleanup 
void Init()
{
m_verify_temps_gone=true;
m_ext=".pawnoff";
m_width=120;
m_nl="\n";
m_bashsep="---XXX---";
m_more="less -X";
//m_date="date -f - \"+%Y-%m-%d:%k:%M:%S\"";
m_date="echo -n `date -f - \"+%Y-%m-%d:%H:%M:%S\"`";
m_now_fn="echo -n `date  \"+%Y-%m-%d-%H-%M-%S\"`";
//m_today="echo -n `date  \"+%Y-%m-%d:%H:%M:%S\" `";
m_today="echo -n `date  \"+%Y-%m-%d\" `";
//MM_ONCE(" fuxed m_today to exclude time wtf",)
m_lexi_day="echo -n `date -f -  +%Y-%m-%d 2>/dev/null `";


m_date_default="date -f - ";
m_base64="base64 -i -d ";
Ss ss;
ss<<"lynx -dump -stdin  --force-html -width="<<m_width<<" ";
// doh
//m_html="lynx -dump -stdin  --force-html -width="+m_width;
m_html=ss.str();
m_qp="tcucodec quote -d ";
m_temp_base="/tmp/mikemail";
m_debug=0;
}
Ss bash_string_input(const StrTy& s )
{
Ss ss; 
const StrTy& bashsep="---XXX---";
const StrTy& nl=m_nl; // "\n";
//const StrTy& more=m_more;
ss<<"tempfunc()"<<nl<<"{"<<nl;
ss<<"cat << "<<bashsep<<nl;
ss<<s;
ss<<nl;
ss<<bashsep<<nl<<"}"<<nl;
return ss;
}



/*

// https://stackoverflow.com/questions/2746168/how-to-construct-a-c-fstream-from-a-posix-file-descriptor

include <ext/stdio_filebuf.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main()
{
    ofstream ofs("test.txt");
    ofs << "Writing to a basic_ofstream object..." << endl;
    ofs.close();

    int posix_handle = fileno(::fopen("test.txt", "r"));

    __gnu_cxx::stdio_filebuf<char> filebuf(posix_handle, std::ios::in); // 1
    istream is(&filebuf); // 2

    string line;
    getline(is, line);
    cout << "line: " << line << std::endl;
    return 0;
}

https://en.cppreference.com/w/cpp/io/c/tmpfile

include <iostream>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
namespace fs = std::filesystem;
 
int main()
{
    std::FILE* tmpf = std::tmpfile();
    std::fputs("Hello, world", tmpf);
    std::rewind(tmpf);
    char buf[6];
    std::fgets(buf, sizeof buf, tmpf);
    std::cout << buf << '\n';
 
    // Linux-specific method to display the tmpfile name
    std::cout << fs::read_symlink(
                     fs::path("/proc/self/fd") / std::to_string(fileno(tmpf))
                 ) << '\n';
}

*/
 
// the cat version may be interretting dollar signs etc 
Ss string_input(const StrTy& s )
{
//if (s.length()<(MAX_BASH_TEXT)) return bash_string_input(s); 
const IdxTy serial=Serial(1);
Ss st; 
// TODO FIXME this needs to verify unique name
//st<<"/tmp/mikemail"<<serial<<".tmptxt";
st<<m_temp_base<<serial<<".tmptxt";

StrTy fn=st.str();
{std::ofstream os(fn);
os<<s; } 
Ss ss;
//const StrTy& bashsep="---XXX---";
const StrTy& nl=m_nl; // "\n";
//const StrTy& more=m_more;
ss<<"tempfunc()"<<nl<<"{"<<nl;
ss<<"cat "<<fn<<nl;
// this removes it without insuring the following command completes- ??? 
if (strncmp(fn.c_str(),m_temp_base.c_str(),m_temp_base.length())==0) 
{// MM_ONCE(" removing the temp files programmatcally ",)
{ ss<<"rm "<<fn<<nl; } }
ss<<"}"<<nl;

return ss;
}
 
// the cat version may be interretting dollar signs etc 
Ss bin_input(const Blob& s )
{
//if (s.length()<(MAX_BASH_TEXT)) return bash_string_input(s); 
const IdxTy serial=Serial(1);
Ss st; 
// TODO FIXME this needs to verify unique name
//st<<"/tmp/mikemail"<<serial<<".tmptxt";
st<<m_temp_base<<serial<<".tmptxt";

StrTy fn=st.str();
{std::ofstream os(fn);
//os<<s; 
os.write(s.ptr(),s.size());
} 

Ss ss;
//const StrTy& bashsep="---XXX---";
const StrTy& nl=m_nl; // "\n";
//const StrTy& more=m_more;
ss<<"tempfunc()"<<nl<<"{"<<nl;
ss<<"cat "<<fn<<nl;
// this removes it without insuring the following command completes- ??? 
if (strncmp(fn.c_str(),m_temp_base.c_str(),m_temp_base.length())==0) 
{// MM_ONCE(" removing the temp files programmatcally ",)
{ ss<<"rm "<<fn<<nl; } }
ss<<"}"<<nl;

return ss;
}






static void*  Bg(void * p)
{
BgParam * _bg= (BgParam*) p;
BgParam bg=*_bg;
delete _bg;
switch (bg.mode)
{
case 0 : {
bg.parent->get_output_with_input( bg.s1,bg.s2);

break; 
} 
default: MM_ERR(" bad mode "<<MMPR(bg.dump())) 

} // mode 
delete bg.parent;
return 0; 
} // Bg


IdxTy Raw(const StrTy & s)
{
Ss ss;
ss<<s;
IdxTy rc=std::system(ss.str().c_str());
return rc;


}
IdxTy  More(const StrTy & s,const IdxTy flags=0)
{
Ss ss=string_input(s);
//m_more="less -X";
const StrTy& more=m_more;
// echo to attempt to stop less from finishing. 
ss<<"tempfunc | "+more+StrTy(";echo");  // more";
//MM_ERR(ss.str());
//MM_ERR(" trying system call now ")
IdxTy rc=std::system(ss.str().c_str());
return rc;
}

IdxTy  FromDoTo(const StrTy & fin, const StrTy & proc, const StrTy & dest)
{
Ss ss;
ss<< StrTy("cat ");
ss<<fin;
ss<<StrTy(" | ");
ss<<proc;
ss<<StrTy(">");
ss<<dest;
return std::system(ss.str().c_str());
}
// https://stackoverflow.com/questions/478898/how-do-i-execute-a-command-and-get-output-of-command-within-c-using-posix
StrTy Dates(const StrTy & dates)
{
Ss ss= string_input(dates);
StrTy f= ss.str()+" tempfunc | "+m_date;
return  DoSomething(f);

}
StrTy LexiDays(const StrTy & dates)
{
StrTy _dates=dates;
if (dates.length()==0) _dates=StrTy("0 days ago");
Ss ss= string_input(_dates);
//MM_ERR("LEXIDAYS "<<MMPR(ss.str()))
StrTy f= ss.str()+" tempfunc | "+m_lexi_day;
return  DoSomething(f);
}
StrTy SerialDays(const StrTy & dates)
{
StrTy _dates=dates;
// echo $(( $( date -d `date -d "0 days ago "  +%Y-%m-%d` +%s) /86400 )) 
if (dates.length()==0) _dates=StrTy("0 days ago");
Ss ss= string_input(dates);
//if (dates.length()==0) ss<<CRLF;
//MM_ERR("SERIALAYS "<<MMPR(ss.str()))
//afck - fudkdate sends stderr to fking null  
//StrTy m_serial_day="echo -n $(( `date -f - +%s ` / 86400 )) ";
//StrTy m_serial_day="echo -n $(( `date -f - +%s ` / 86400 )) ";
/*
# autogenrated with mkscriptu at Thu Jul 30 17:16:17 EDT 2020 
sht=`date -f - +%s 2>/dev/null  `  2>/dev/null  # / 86400 
#echo -n $(( `date -f - +%s ` / 86400 ))
echo "$sht/86400"  | bc -l 2>/dev/null
#date $* 2>/dev/null
*/

const bool use_fdate=false;
if (use_fdate)
{
StrTy m_serial_day="fudddate";
StrTy f= ss.str()+" tempfunc | "+m_serial_day;
return  DoSomething(f);
}
StrTy m_serial_day="date -f - +%s 2>/dev/null";
StrTy f= ss.str()+" tempfunc | "+m_serial_day;
StrTy sd= DoSomething(f);
int x=0; 
{ Ss ss; ss<<sd; ss>>x; } 
{ Ss ss; ss<<int(x/86400);
return ss.str(); } 
} // SerialDays




StrTy DatesDefault(const StrTy & dates)
{
Ss ss= string_input(dates);
StrTy f= ss.str()+" tempfunc | "+m_date_default;
return  DoSomething(f);

}

void Base64( char ** result, IdxTy * len, const StrTy & d)
{
Ss ss= string_input(d);
StrTy f= ss.str()+" tempfunc | "+m_base64;
DoSomething( result, len,  f );
}
void Html( char ** result, IdxTy * len, const StrTy & d)
{
Ss ss= string_input(d);
StrTy f= ss.str()+" tempfunc | "+m_html;
DoSomething( result, len,  f );
}
void Qp( char ** result, IdxTy * len, const StrTy & d)
{
Ss ss= string_input(d);
StrTy f= ss.str()+" tempfunc | "+m_qp;
DoSomething( result, len,  f );
}

void GetOutput( char ** result, IdxTy * len, const StrTy & d, const StrTy & cmd)
{
//Ss ss= string_input(d);
//StrTy f= ss.str()+" tempfunc | "+cmd;
StrTy f= cmd;
//MM_ERR(" XXXX "<<MMPR(f))
DoSomething( result, len,  f );
}

void GetOutputWithInput( char ** result, IdxTy * len, const StrTy & d, const StrTy & cmd)
{
Ss ss= string_input(d);
StrTy f= ss.str()+" tempfunc | "+cmd;
//MM_ERR(" XXXX "<<MMPR(f))
//StrTy f= cmd;
DoSomething( result, len,  f );
}

void GetOutputWithInput( char ** result, IdxTy * len, const Blob & d, const StrTy & cmd)
{
Ss ss= bin_input(d);
StrTy f= ss.str()+" tempfunc | "+cmd;
//MM_ERR(" XXXX "<<MMPR(f))
//StrTy f= cmd;
DoSomething( result, len,  f );
}






StrTy DoSomething( const StrTy & something )
{
StrTy s; 
//MM_ERR(MMPR(something))
const char * cmd= something.c_str();
const IdxTy BUFSZ=OUTBUFSZ; // 1<<12;
    char* buffer = new char [BUFSZ];
// TODO no idea on threading here, hang on full buffer ? 
    std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error("popen() failed!");
    try {
        while (fgets(buffer, BUFSZ, pipe) != NULL) {
            result += buffer;
        }
    } catch (...) {
        pclose(pipe);
	delete buffer;
        throw;
    }
    pclose(pipe);
	s=result;
	delete buffer;
return s; 
}

IdxTy  DoSomething(char ** result, IdxTy * len,  const StrTy & something )
{
StrTy s; 
if( Bit(m_debug,0)) MM_ERR(MMPR(something))
const char * cmd= something.c_str();
const IdxTy BUFSZ=OUTBUFSZ; // 1<<12;
IdxTy cnt=0;
//int fd[2]; // mjm
//pipe(fd);
//dup2(2,fd[1]);
//MM_ERR(" doing "<<MMPR(something))
std::vector<char> buffer(BUFSZ);
std::vector<char> errbuffer(BUFSZ);
//char * buffer= new char[BUFSZ];
//    char buffer[BUFSZ];
// TODO no idea on threading here, hang on full buffer ? 
   // std::string result = "";
    FILE* pipe = popen(cmd, "r");
    if (!pipe) throw std::runtime_error(" mjm_pawnoff popen() failed!");
    try {
	while ( !feof(pipe)&&!ferror(pipe)) { 
	IdxTy sz=fread(&buffer[cnt],1,buffer.size()-cnt,pipe);
	cnt+=sz;
	IdxTy seff=cnt+BUFSZ;
	if (seff>buffer.size()) buffer.resize(seff);
}  

//        while (fgets(buffer, BUFSZ, pipe) != NULL) {
//            result += buffer;
 //       }
    } catch (...) {
        pclose(pipe);
        throw;
    }
//close(fd[0]); close(fd[1]);
    pclose(pipe);
char * b = new char[cnt+STRING_GUARD];
memcpy(b,&buffer[0],cnt);
*result=b;
*len=cnt;
//	s=result;
return 0; 
}

// MEMBERS
bool m_verify_temps_gone;

}; // mjm_pawnoff


#endif // MJM_pawnoff_H__

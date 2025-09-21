
#ifndef MJM_datascope_H__
#define MJM_datascope_H__
 
#include "mjm_globals.h"
// 
#include "mjm_thread_util.h"
#include "mjm_canned_methods.h"
 
#ifndef NO_DSCOPE_GRAPHICS
// complains if included later doh 
// although this should only be in glut_scope NOT
// here in datascope. 
// save is all deprecated crap needs to be fixed 
//#include "mjm_glut_saver.h" 
#include "mjm_glut_scope_ii.h"
#include "mjm_glut_rags.h" 
//#include "mjm_svg_render.h" 

// the pool makes a problem... 
// for the rags objects
#include "mjm_glut_helpers.h"
#include "mjm_object_pool.h"
#endif

#include "mjm_data_model_error_log.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include "mjm_strings.h"
#include "mjm_string_kvp.h"


#include "mjm_dscope_pktfifo.h"
#include "mjm_dscope_rawfifo.h"
#include "mjm_dscope_dgram.h"

// g++  -std=gnu++11 -DTEST_MJM_GLUT_SAVER -I. -I../../mjm/hlib -I../../mjm/num  -gdwarf-3 -O0  -x c++ mjm_glut_saver.h  -o mjm_glut_saver.out -lpthread -lreadline  -lGL -lGLU -lglut -lpng -lavcodec -lswscale -lavutil 


#include "mjm_cli_ui.h"

#include "mjm_tokenized_collections.h"


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
made from wizard script on Thu Mar 16 07:24:28 EDT 2023
datascope

// https://tex.stackexchange.com/questions/254610/how-can-i-use-bibtex-to-cite-a-software
needs to come from "about" too needs version in name etc.  
@software{mjm_cpp_datascope,
  author = {Michael J Marchywka},
  title = {datascope},
abstract=(),
institution={},
license={Knowledge sir should be free to all },
publisher={Mike Marchywka},
email={marchywka@hotmail.com},
authorid={orcid.org/0000-0001-9237-455X},
  filename = {datascope},
  url = {},
  version = {0.0.0},
  date-started = {Thu Mar 16 07:24:28 EDT 2023}
}

*/


////////////////////////////////////////////////////////////////

class datascope_params : public mjm_logic_base
{
typedef mjm_logic_base Super;
typedef Super::Tr Tr;
typedef  Tr::D D;
typedef Tr::IdxTy IdxTy;
typedef Tr::StrTy StrTy;
typedef Tr::Ss Ss;

public:
datascope_params( const StrTy & nm) : Super(nm) {}
datascope_params() : Super() {}
// should be geneated code, do not edit  to make every entry for one name  
bool log_commands() const { return m_map.get_bool("log_commands",!true); }
bool exit_on_err() const { return m_map.get_bool("exit_on_err",!true); }
//StrTy protrait_eol() const { return m_map.get_string("protrait_eol","\r"); }
IdxTy omaxmax() const { return m_map.get_uint("omaxmax",5); } // // 100;
StrTy ragged_params() const { return m_map.get_string("ragged_params","."); }
IdxTy maxcnt() const { return m_map.get_uint("maxcnt",100); } // // 100;
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
//ss<<"protrait_eol="<<protrait_eol().c_str()[0]<<sep;
ss<<"omaxmax"<<omaxmax()<<sep;
ss<<"ragged_params"<<ragged_params()<<sep;
ss<<"maxcnt"<<maxcnt()<<sep;
return ss.str();
}


}; // datascope_params


namespace datascope_traits
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
typedef IdxTy FlagTy;
typedef  data_model_error_log Dmel; 
// NB NOTE wow reducing the size to match need cut memory usage
// and bumped speed a lot with compact DS. 
typedef uint64_t KeyCode;
}; // 

typedef Tr::IdxTy IdxTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::FlagTy FlagTy;
// ns types should come from trits of something 
typedef std::vector<StrTy> Words;


}; // datascope_traits
///////////////////////////////////////////////////////////////




class mjm_datascope 
{
typedef datascope_traits::Tr  Tr;
//typedef linc_graph_traits::util  Util;
typedef mjm_datascope Myt;
typedef mjm_cli_ui<Myt> CliTy;
typedef Tr::IdxTy IdxTy;
typedef Tr::FlagTy FlagTy;
typedef Tr::D D;
typedef Tr::Ss Ss;
typedef Tr::IsTy IsTy;
typedef Tr::OsTy OsTy;
typedef Tr::Ofs Ofs;
typedef Tr::MyBlock  MyBlock;

// TYPEDEFS
typedef std::map<StrTy, StrTy> LocalVar;
typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef void (Myt:: * CompleteFunc) ( CliTy::list_type & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

typedef mjm_canned_methods Canned;
typedef datascope_params ParamGlob;
typedef mjm_logic_base VariableStore;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line;
typedef std::map<StrTy, Ragged> RaggedMap;
typedef mjm_tokenized_ragged_table TokRagged;
typedef std::map<StrTy, TokRagged> TokRaggedMap;

typedef std::vector<StrTy> Words;
typedef data_model_error_log Dmel;

//typedef string_tokenizer Tokenizer;
////////////////////////////////////////////////////
//typedef mjm_tax_tree TaxTree;
//typedef std::map<StrTy, TaxTree> TaxTrees;
typedef std::map<IdxTy,IdxTy> Mii;
typedef std::map<StrTy, Mii> MiiMap;
//typedef std::map<StrTy,TaxTree> TaxTrees;

typedef std::vector<IdxTy> LocTy;

typedef mjm_string_base_params<Tr> BaseParams;
// traits? 
enum {BAD=~0};

#ifndef NO_DSCOPE_GRAPHICS
// changed the indexing from Line to normal StrTy
// using first line "# id stuff ... name"
typedef mjm_glut_rags<Tr> RagScene;
//typedef mjm_glut_rags_map<Tr,Line> SceneMap;
typedef mjm_glut_rags_map<Tr,StrTy> SceneMap;
typedef mjm_glut_scope_ii<Myt,Tr> GlutScope;
typedef mjm_glut_saver<Tr> Saver;
//typedef mjm_svg_render<Tr> SvgRender;
typedef typename Saver::save_params_type Sp;
// added for recycling 
// while this is the model it is defeined in view stuff
typedef mjm_glut_helpers<Tr> GlutUtil;
typedef typename  GlutUtil::junk_bin_t ModelInfo;
// apparently this is just a pool shared among
// index_buffers that each claim their own 
// slots . just a memory coherence thing ... 
typedef mjm_object_pool<Tr,ModelInfo> ModelPool;


#endif

typedef mjm_dscope_dgram<Tr> DgramID;
typedef mjm_dscope_pktfifo<Tr> PktFifoID;
typedef mjm_dscope_rawfifo<Tr> RawFifoID;
typedef DgramID::data_type Data;


typedef mjm_thread_util<Tr> ThreadTy;
typedef pthread_t ThreadId;
typedef ThreadTy::ThParam<Myt> ThParam;

public:
// FIXME doh put this somwhere lol 
int myatoi(const StrTy & s ) const { return Canned::myatoi(s.c_str()); } 
int myatoi(const char * c) const { return Canned::myatoi(c); }


// API
public :
mjm_datascope():m_dmel(new Dmel()) {Init();}
mjm_datascope(int argc,char **_args) : m_dmel(new Dmel())
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
if (i==istart) { MM_ERR(" did nothing with "<<args[i]) ++i;  } 

}
}
~mjm_datascope()
{
//clear_handlers();
delete m_dmel;
}
////////////////////////////////////////////////////////
// command block

// this should be in the parameters map, nothing special here... 
 void configi(IdxTy & dest, const StrTy & cmd, int  & i, int argc, char ** args)
{
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if ( s==cmd)
{
++i;
//const StrTy nm=StrTy(args[i]);
dest=myatoi(args[i]);
MM_ERR(" setting "<<s<<" to "<<dest)
++i; // consume param and cmd
}
}
 void cmdlcmd( int  & i, int argc, char ** args)
{
const bool confirm=true;
if (argc<=i) return; 
const StrTy s=StrTy(args[i]);
if (s=="-source") { ++i; command_modef(args[i]); ++i; }
if (s=="-start") { ++i; start_date(StrTy(args[i])); ++i; }
if (s=="-end") { ++i; end_date(StrTy(args[i])); ++i; }
if (s=="-cmd") { ++i; command_mode(StrTy(args[i])); ++i; }
if (s=="-quit") { MM_ERR(" QUIT "<<MMPR4(i,argc,args[i],s))  ++i; clean_up(); }
if (s=="-about") { ++i; about(); }
if (confirm) {}
} // cmdlcmd
void arg_cmd(int & i,  char ** args, const IdxTy n, const char *  base, const bool confirm)
{
StrTy cmd=StrTy(base);
for (IdxTy j=0; j<n ; ++j)  { ++i; cmd=cmd+StrTy(" ")+StrTy(args[i]); }
if(confirm) {MM_ERR(" executing  "<<cmd) } 
command_mode(cmd);
++i; 
} 
void start_date(const StrTy & d) { m_flp.set("start_date",d); }
void end_date(const StrTy & d) { m_flp.set("end_date",d); }

void command_modef(const char * fn)
{ std::ifstream fin(fn); CommandInterpretter li(&fin); command_mode(li); }
void command_mode() { CommandInterpretter li(&std::cin); command_mode(li); }
void command_mode(const StrTy & cmd) 
{ CommandInterpretter li; li.set(cmd,1); command_mode(li); }

CmdMap m_cmd_map;
CompMap m_comp_map;
 void cli_cmd( CliTy::list_type & choices,  const char * frag)
{
//MM_ERR("cli_cmd"<<MMPR(frag))
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}

}
 void cli_param( CliTy::list_type & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
const StrTy cmd=CliTy::word(StrTy(_cmd),0);
auto ii=m_comp_map.find(cmd);
if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag); 
}


/*
void cmd_stream_edit_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_stream_edit_fasta( cip ,  lv ); 
}

void cmd_read_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_read_fasta(cip ,  lv,  m_fasta_map  ); }

void cmd_write_fasta(Cip & cip , LocalVar & lv ) { Canned::cmd_write_fasta(cip ,  lv,  m_fasta_map  );

}
*/


// if (m_p_app) if (m_p_app->glut_key(key,x,y,shift,ctrl,alt)) return;
IdxTy glut_key(const int key, const int x, const int y, const bool shift, const bool ctrl, const bool alt) 
{


return 0;
} // glut_key
IdxTy glut_key(const int key, const int x, const int y) 
{


return 0;
} // glut_key


bool flag_bit(const IdxTy flag, const IdxTy bit, const bool pol=true)
{
const bool x=((flag&(1<<bit))!=0);
return pol?x:!x;
}

#ifndef NO_DSCOPE_GRAPHICS
void clear(const IdxTy flags=0)
{
MM_ERR(" app clearing ")
// this does nothing alone
//m_scenes.clear();
m_glut.clear(flags);
MM_ERR(" scenese clearing ")

m_scenes.clear();
} // clear
void cmd_clear(Cip & cip , LocalVar & lv )
{

const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
clear();
} // cmd_clear 
#define CC(x) myatoi(cip.wif(x)) 
void cmd_scope(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
GlutScope & x=m_glut;
if (cmd=="dump") { auto wtf=(x.dump()) ; MM_ERR(wtf) }
if (cmd=="saver") {x.set_saver(&m_saver);  }
if (cmd=="sig") {x.signal(s,flags);  }
if (cmd=="nsaver") {x.set_saver(0);  }
if (cmd=="bounds") {x.bound_box(CC(2),CC(3),CC(4),CC(5),0);  }
if (cmd=="spr") {Sp sp(s); x.set_saver_params(sp);  }
if (cmd=="launch") { auto wtf=(x.launch(s)); MM_ERR(wtf);  }
//if (cmd=="stream") { MM_ERR(x.stream(s)) }
// no real reson to call this directly.. 
if (cmd=="start") { auto wtf=(x.start(s)); MM_ERR(wtf)  }
if (cmd=="wait") { StrTy x; std::cin>> x ; }
MM_ERR(" return from cmd_scope")
} // cmd_scope


void cmd_saver(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
//GlutScope & x=m_glut;
Saver & x= m_saver;
//if (cmd=="main") { MM_ERR(x.main(0,0)) }
if (cmd=="main") {auto wtf=(x.main(0,0)); MM_ERR(wtf); }
//if (cmd=="dump") { MM_ERR(x.dump()) }
// not important here but for text search LOL WTF
if (cmd=="dump") { auto wtf=x.dump() ; MM_ERR(wtf);}
if (cmd=="set") { x.set(s,cip.wif(3)); MM_ERR("set");}
//if (cmd=="launch") { MM_ERR(x.launch(s,flags)) }
//if (cmd=="stream") { MM_ERR(x.stream(s)) }
//if (cmd=="start") { MM_ERR(x.start(s)) }
//if (cmd=="wait") { StrTy x; std::cin>> x ; }
MM_ERR(" return from cmd_dgram")
} // cmd_saver


void cmd_display(Cip & cip , LocalVar & lv )
{
const StrTy src=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,src,s,flags))
// this assures we can leave scope and let user delete 
ThParam* tp= new ThParam(this, &Myt::Displayt, s,flags ) ;
tp->invoke();
MM_ERR(" return from "<<MMPR(__FUNCTION__) )
} // cmd_display

void cmd_display_fifo(Cip & cip , LocalVar & lv )
{
const StrTy src=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,src,s,flags))
// this assures we can leave scope and let user delete 
ThParam* tp= new ThParam(this, &Myt::Displaytfifo, s,flags ) ;
tp->invoke();
MM_ERR(" return from "<<MMPR(__FUNCTION__) )
} // cmd_display_fifo



void cmd_modify(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.cmd();
const StrTy sin=cip.p1;
//const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(2));
MM_ERR(MMPR4(__FUNCTION__,sin,cmd,flags))
m_glut.modify(sin,flags);

} // cmd_modify

void cmd_show(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.cmd();
const StrTy sin=cip.p1;
//const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(2));
MM_ERR(MMPR4(__FUNCTION__,sin,cmd,flags))
StrTy s=m_glut.show(sin,flags);
MM_MSG(MMPR(s))
} // cmd_show



#endif

void cmd_dgram(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
//GlutScope & x=m_glut;
DgramID & x= m_dgrams;
//if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="dump") { const auto wtf=x.dump() ; MM_ERR(wtf) }
//if (cmd=="launch") { MM_ERR(x.launch(s,flags)) }
if (cmd=="launch") { const auto wtf=(x.launch(s,flags)); MM_ERR(wtf)  }
//if (cmd=="stream") { MM_ERR(x.stream(s)) }
//if (cmd=="start") { MM_ERR(x.start(s)) }
//if (cmd=="wait") { StrTy x; std::cin>> x ; }
MM_ERR(" return from cmd_dgram")
} // cmd_dgram



void cmd_pktfifo(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
//GlutScope & x=m_glut;
PktFifoID & x= m_pktfifo;
if (cmd=="dump") { MM_ERR(x.dump()) }
//if (cmd=="launch") { MM_ERR(x.launch(s,flags)) }
if (cmd=="launch") { const auto wtf=(x.launch(s,flags)); MM_ERR(wtf);  }
//if (cmd=="stream") { MM_ERR(x.stream(s)) }
//if (cmd=="start") { MM_ERR(x.start(s)) }
//if (cmd=="wait") { StrTy x; std::cin>> x ; }
MM_ERR(" return from cmd_dgram")
} // cmd_pktfifo

void cmd_rawfifo(Cip & cip , LocalVar & lv )
{
const StrTy cmd=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,cmd,s,flags))
//GlutScope & x=m_glut;
RawFifoID & x= m_rawfifo;
// DOH turning off the error messages prevents this from getting
// evaluted ROFL 
//if (cmd=="dump") { MM_ERR(x.dump()) }
//if (cmd=="launch") { MM_ERR(x.launch(s,flags)) }
if (cmd=="dump") { auto wtf=x.dump() ; MM_ERR(wtf) }
if (cmd=="launch") { auto wtf=x.launch(s,flags); MM_ERR(wtf)  }

//if (cmd=="stream") { MM_ERR(x.stream(s)) }
//if (cmd=="start") { MM_ERR(x.start(s)) }
//if (cmd=="wait") { StrTy x; std::cin>> x ; }
MM_ERR(" return from cmd_rawfifo")
} // cmd_rawfifo
// send stdin to a  ragged and send somewhere  
// this cn fork is the commands are not coming from stdin,
// not sure if source works from cmdline...
void cmd_send(Cip & cip , LocalVar & lv )
{
const StrTy src=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,src,s,flags))
// this assures we can leave scope and let user delete 
ThParam* tp= new ThParam(this, &Myt::Sendt, s,flags ) ;
tp->invoke();
MM_ERR(" return from "<<MMPR(__FUNCTION__) )
} // cmd_send
void cmd_sendfifo(Cip & cip , LocalVar & lv )
{
const StrTy src=cip.p1;
const StrTy s=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
MM_ERR(MMPR4(__FUNCTION__,src,s,flags))
// this assures we can leave scope and let user delete 
ThParam* tp= new ThParam(this, &Myt::Sendtfifo, s,flags ) ;
tp->invoke();
MM_ERR(" return from "<<MMPR(__FUNCTION__) )
} // cmd_sendfifo

/*
void cmd_zymo_rags(Cip & cip , LocalVar & lv )
{

const StrTy ragin=cip.p1;
const IdxTy flags=myatoi(cip.p2); // wif(3);
const StrTy cmd1=cip.wif(3);
const StrTy cmd2=cip.wif(4);
const Ragged & qio=m_ragged_map[ragin];
const Ragged & params=m_ragged_map[cmd2];
const IdxTy nval=myatoi(cmd1);
const IdxTy maxcnt=m_flp.maxcnt();
const bool output_latex=(((1<<0)&flags)!=0);
const bool  output_rank_table=(((1<<1)&flags)!=0);
const bool  output_ranks=(((1<<2)&flags)!=0)&&!output_rank_table;
const bool output_summary=(((1<<3)&flags)!=0);

} // zymo_rags
//void cmd_add_to_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_add_to_fasta(cip ,  lv,  m_fasta_map  ); }

//void cmd_zymo_merge_fasta(Cip & cip , LocalVar & lv )
//{ Canned::cmd_zymo_merge_fasta(cip ,  lv,  m_fasta_map, m_ragged_map  ); }

void cmd_write_svg_ragged(Cip & cip , LocalVar & lv ) 
{
const StrTy cmd=cip.cmd();
const StrTy fn=cip.p1;
const StrTy name=cip.p2;
const IdxTy flags=myatoi(cip.wif(3));
const StrTy prag=(cip.wif(4));
Ragged & r=m_ragged_map[name];
Ragged & pr=m_ragged_map[prag];
MM_ERR(MMPR4(cmd,fn,name,flags)<<MMPR3(prag,pr.size(),r.size()))

}

*/

void cmd_transpose_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_transpose_ragged(cip ,  lv, m_ragged_map  ) ; } // transpose
void cmd_read_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_read_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_write_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_write_ragged( cip ,  lv, m_ragged_map  ) ; }
void cmd_dump_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_dump_ragged( std::cout, cip ,  lv, m_ragged_map  ) ; }



void cmd_add_ragged(Cip & cip , LocalVar & lv ) 
{ Canned::cmd_add_ragged( cip ,  lv, m_ragged_map  ) ; }
//void cmd_transpose_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_transpose_ragged(cip ,  lv, m_tokragged_map  ) ; } // transpose
//void cmd_read_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_read_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_write_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_write_ragged( cip ,  lv, m_tokragged_map  ) ; }
//void cmd_add_tragged(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_add_ragged( cip ,  lv, m_tokragged_map  ) ; }



//void cmd_tt(Cip & cip , LocalVar & lv ) 
//{ Canned::cmd_tt( cip ,  lv, m_tax_trees ) ; }


void cmd_source(Cip & cip , LocalVar & lv ) 
{ const char * fn=cip.p1.c_str();  command_modef(fn); }

void cmd_quit(Cip & cip , LocalVar & lv )
{ clean_up(); return;  }
void cmd_cm(Cip & cip , LocalVar & lv )
{ dump_cm(); return;  }
void cmd_banner(Cip & cip , LocalVar & lv )
{ config_banner(); return;  }

void cmd_set_param(Cip & cip , LocalVar & lv )
{  //kkkkkkkkkkk 
//const StrTy cmd=cip.cmd();
//const StrTy fn=cip.p1;
//const StrTy name=cip.p2;
//const IdxTy flags=myatoi(cip.wif(3));
//const StrTy prag=(cip.wif(4));
 if (cip.li().cmd_ok(3))  m_flp.set(cip.li().cmd_set());
        return;  }
void cmd_get_param(Cip & cip , LocalVar & lv )
{
 if (cip.li().cmd_ok(2))
        std::cout<<lv["local_label"]<<" "<<cip.li().word(1)<<" "
                        <<m_flp.get(cip.li().word(1))<<CRLF;
         }





void cmd_help(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_cmd_map)
{
MM_ERR(MMPR((*ii).first))

} 

}


void _quiet() { mjm_global_flags::mm_err_enable=!true; }
void _nquiet() { mjm_global_flags::mm_err_enable=true; }

void cmd_err_msg(Cip & cip , LocalVar & lv ) 
{
const IdxTy flags=myatoi(cip.p1); // wif(3);
if (Bit(flags,0)) _nquiet(); else _quiet(); 

}


void cmd_set(Cip & cip , LocalVar & lv ) 
{
const StrTy& v=(cip.p1); // wif(3);
const IdxTy flags=myatoi(cip.p2); // wif(3);
//if (Bit(flags,0)) _nquiet(); else _quiet(); 
if (v=="debug" ) 
	{ m_debug=flags; MM_ERR(" setting "<<MMPR3(v,flags,m_debug)); }
else if (v=="done" ) 
	{ m_done=flags; MM_ERR(" setting "<<MMPR3(v,flags,m_done)); }
else if (v=="stop-send" ) 
	{ m_stop_send=flags; MM_ERR(" setting "<<MMPR3(v,flags,m_stop_send)); }

MM_ERR(" done with cmd_set")
}


void cmd_list(Cip & cip , LocalVar & lv ) 
{
MM_LOOP(ii,m_ragged_map) { MM_MSG("m_ragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_LOOP(ii,m_tokragged_map) { MM_MSG("m_tokragged_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_fasta_map) { MM_MSG("m_fasta_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_dig_map) { MM_MSG("m_dig_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_pheno_map) { MM_MSG("m_pheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_hbpheno_map) { MM_MSG("m_hbpheno_map "<<MMPR2((*ii).first,(*ii).second.size())) }
//MM_LOOP(ii,m_tax_trees) { MM_MSG("m_tax_trees "<<MMPR2((*ii).first,(*ii).second.size())) }
MM_MSG(" configuration "<<m_flp.to_string())
dump_cm();


}


static const IdxTy & bad()  { static IdxTy i=~0U; return i; } 

////////////////////////////////////////////
void SetupCmdMap()
{
m_cmd_map[StrTy("?")]=&Myt::cmd_help;
m_cmd_map[StrTy("help")]=&Myt::cmd_help;
m_cmd_map[StrTy("list")]=&Myt::cmd_list;
m_cmd_map[StrTy("err-msg")]=&Myt::cmd_err_msg;
m_cmd_map[StrTy("source")]=&Myt::cmd_source;

m_cmd_map[StrTy("quit")]=&Myt::cmd_quit;
m_cmd_map[StrTy("cm")]=&Myt::cmd_cm;
m_cmd_map[StrTy("banner")]=&Myt::cmd_banner;
m_cmd_map[StrTy("set-param")]=&Myt::cmd_set_param;
m_cmd_map[StrTy("get-param")]=&Myt::cmd_get_param;

m_cmd_map[StrTy("read-ragged")]=&Myt::cmd_read_ragged;
m_cmd_map[StrTy("dump-ragged")]=&Myt::cmd_dump_ragged;
m_cmd_map[StrTy("write-ragged")]=&Myt::cmd_write_ragged;
m_cmd_map[StrTy("transpose-ragged")]=&Myt::cmd_transpose_ragged;
m_cmd_map[StrTy("add-ragged")]=&Myt::cmd_add_ragged;
m_cmd_map[StrTy("dgram")]=&Myt::cmd_dgram;
m_cmd_map[StrTy("pktfifo")]=&Myt::cmd_pktfifo;
m_cmd_map[StrTy("rawfifo")]=&Myt::cmd_rawfifo;
m_cmd_map[StrTy("send")]=&Myt::cmd_send;
m_cmd_map[StrTy("sendfifo")]=&Myt::cmd_sendfifo;
m_cmd_map[StrTy("set")]=&Myt::cmd_set;

#ifndef NO_DSCOPE_GRAPHICS
m_cmd_map[StrTy("modify")]=&Myt::cmd_modify;
m_cmd_map[StrTy("show")]=&Myt::cmd_show;
m_cmd_map[StrTy("scope")]=&Myt::cmd_scope;
m_cmd_map[StrTy("clear")]=&Myt::cmd_clear;
m_cmd_map[StrTy("saver")]=&Myt::cmd_saver;
m_cmd_map[StrTy("display")]=&Myt::cmd_display;
m_cmd_map[StrTy("displayfifo")]=&Myt::cmd_display_fifo;
#endif

//m_cmd_map[StrTy("read-tragged")]=&Myt::cmd_read_tragged;
//m_cmd_map[StrTy("write-tragged")]=&Myt::cmd_write_tragged;
//m_cmd_map[StrTy("transpose-tragged")]=&Myt::cmd_transpose_tragged;
//m_cmd_map[StrTy("add-tragged")]=&Myt::cmd_add_tragged;

m_cmd_map[StrTy("string-ragged")]=&Myt::cmd_read_ragged;
//m_cmd_map[StrTy("write-svg-ragged")]=&Myt::cmd_write_svg_ragged;
//m_cmd_map[StrTy("read-dig")]=&Myt::cmd_read_dig;
//m_cmd_map[StrTy("read-fasta")]=&Myt::cmd_read_fasta;
//m_cmd_map[StrTy("stream-edit-fasta")]=&Myt::cmd_stream_edit_fasta;
//m_cmd_map[StrTy("write-fasta")]=&Myt::cmd_write_fasta;
//m_cmd_map[StrTy("add-to-fasta")]=&Myt::cmd_add_to_fasta;
//m_cmd_map[StrTy("zymo-merge-fasta")]=&Myt::cmd_zymo_merge_fasta;
//m_cmd_map[StrTy("zymo-rags")]=&Myt::cmd_zymo_rags;
//m_cmd_map[StrTy("group-stats")]=&Myt::cmd_group_stats;

//m_cmd_map[StrTy("linc-graph")]=&Myt::cmd_linc_graph;

//m_cmd_map[StrTy("query-aln")]=&Myt::cmd_query_aln;
//m_cmd_map[StrTy("tt")]=&Myt::cmd_tt;

}
 
void command_mode(CommandInterpretter & li)
{
SetupCmdMap();
//TaxTree & tt = m_tax_tree;
StrTy local_label="tat";
//typedef void (Tsrc::* TargCmd)( ListTy & choices,  const char * frag);
//typedef void (Tsrc::* TargParam)( ListTy & choices, const char *  cmd, const char * frag);
m_cli.set_target(*this);
//void set_command_handler(TargCmd * p ) { m_targ_cmd=p; }
//void set_param_handler(TargParam * p ) { m_targ_param=p; }
m_cli.set_command_handler(&Myt::cli_cmd);
m_cli.set_param_handler(&Myt::cli_param);
//std::vector<StrTy> some;
//some.push_back(StrTy("load-tax-nodes"));
//m_cli.set_list(some);
m_cli.activate();
LocalVar mloc;
mloc["local_label"]=local_label;
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
//MM_ERR(" processing "<<li.dump())
if (m_flp.log_commands()) 
		if (m_dmel!=0) (*m_dmel).event("normal command",li.dump(),"NOTDATA");
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd.c_str()[0]=='#' ) continue; 

if (m_cmd_map.find(cmd)!=m_cmd_map.end())
{
 CommandInterpretterParam  cip(li); 
(this->*m_cmd_map[cmd])(cip,mloc);
if ( m_done){ final_exit();  return; } 
continue;

}
const StrTy p1=(sz>1)?li.word(1):StrTy("");
const StrTy p2=(sz>2)?li.word(2):StrTy("");
if (cmd=="about") { about();  continue; } 
//if (cmd=="solve") { solve(); continue; } 
if (cmd=="print") { MM_MSG(li.line())  continue; } 
if (cmd=="err") { MM_ERR(li.line())  continue; } 
if (cmd=="status") { MM_STATUS(li.line())  continue; } 
if (cmd=="set-param") { if (li.cmd_ok(3))  m_flp.set(li.cmd_set());  continue; } 
if (cmd=="get-param") { if (li.cmd_ok(2))  std::cout<<local_label<<" "<<li.word(1)<<" "<<m_flp.get(li.word(1))<<CRLF;  continue; } 
if (cmd=="banner") { config_banner();  continue; } 
if (cmd=="cm") { dump_cm();  continue; } 
if (cmd=="test") { test(li);  continue; } 
if (cmd=="quit") { clean_up(); return; } 
if (cmd.length()==0) { continue;}
if (cmd.c_str()[0]=='#') { continue; }
MM_ERR(" unrecignized command "<<li.line()<<" "<<li.dump())
if (m_dmel!=0) (*m_dmel).event("bad command ",li.line(),"NOTDATA");
if (m_flp.exit_on_err())
{
MM_ERR(" quiting "<<MMPR(m_flp.exit_on_err()))
clean_up();
return; 

}
}


} //command_mode
// when exiting the interpretter
void clean_up()
{
m_done=true;

} 
void final_exit()
{
MM_ERR(" final exit ") 
#ifndef NO_DSCOPE_GRAPHICS
if (m_glut.alive()) m_glut.exit_loop();
IdxTy nwait=0;
while (m_glut.alive())
{
if (!m_glut.alive()) break;
if (nwait) MM_ERR(" wait for glt die "<<MMPR(nwait))
sleep(1); 
++nwait;
} 
#endif

} // final_exit

void about()
{
Ss ss;
ss<<" mjm_datascope built  "<<__DATE__<<" "<<__TIME__<<CRLF;
ss<<" Mike Marchywka marchywka@hotmail.com started  Thu Mar 16 07:24:28 EDT 2023 "<<CRLF;
ss<<" A free standing application to display data from multiple sources  "<<CRLF;
ss<<" with a nominally consistent but limited format.  "<<CRLF;
ss<<" Intended largely for devlopement or monituring applications  "<<CRLF;
ss<<"  continuation of a work started in 2007 or so ( probably in  glutp.zip  in my   "<<CRLF;
ss<<"  archives ) for simple real time plotting     "<<CRLF;
ss<<"  beyond ASCII art and motivated now by freefem and       "<<CRLF;
ss<<"  interest in simple diffusive systems.       "<<CRLF;
ss<<" See also similar applications, "<<CRLF;
ss<<" https://www.thregr.org/wavexx/software/trend/ Iñaki Ucar iucar@fedoraproject.org"<<CRLF ; 
ss<<" glut skeleton from  Stanford University, CS248, Fall 2000"<<CRLF;
ss<<" various links in source code,"<<CRLF;

ss<<" http://reality.sgi.com/mjk_asd/spec3/spec3.html"<<CRLF;
ss<<" https://community.khronos.org/t/what-are-the-codes-for-arrow-keys-to-use-in-glut-keyboard-callback-function/26457"<<CRLF;
ss<<" https://github.com/cirosantilli/cpp-cheat/blob/19044698f91fefa9cb75328c44f7a487d336b541/ffmpeg/encode.c"<<CRLF;
ss<<" https://github.com/cirosantilli/cpp-cheat/blob/19044698f91fefa9cb75328c44f7a487d336b541/png/open_manipulate_write.c"<<CRLF;
ss<<" https://lmb.informatik.uni-freiburg.de/people/reisert/opengl/doc/glOrtho.html"<<CRLF;
ss<<" https://openshot.org/files/libopenshot/FFmpegUtilities_8h_source.html"<<CRLF;
ss<<" https://stackoverflow.com/questions/14378/using-the-mouse-scrollwheel-in-glut"<<CRLF;
ss<<" https://stackoverflow.com/questions/3191978/how-to-use-glut-opengl-to-render-to-a-file"<<CRLF;
ss<<" https://stackoverflow.com/questions/538661/how-do-i-draw-text-with-glut-opengl-in-c"<<CRLF;
ss<<" https://tex.stackexchange.com/questions/254610/how-can-i-use-bibtex-to-cite-a-software"<<CRLF;
ss<<" https://www3.ntu.edu.sg/home/ehchua/programming/opengl/CG_Examples.html"<<CRLF;
ss<<" https://www.cs.cmu.edu/afs/cs/academic/class/15492-f07/www/pthreads.html#SYNCHRONIZATION"<<CRLF;
ss<<" https://www.ece.lsu.edu/koppel/gpup/code/gpup/demo-4-simple-ogl.cc.html"<<CRLF;
ss<<" https://www.geeksforgeeks.org/udp-server-client-implementation-c/"<<CRLF;
ss<<" https://www.khronos.org/opengl/wiki/Viewing_and_Transformations"<<CRLF;
ss<<" https://www.thregr.org/wavexx/software/trend/"<<CRLF;
ss<<" http://www.gnu.org/licenses/>."<<CRLF;


std::ostream & os=std::cout;
os<<ss.str();

}


// tests, document what this crap should do lol
// for hierarchaila commands 
void test( CommandInterpretter & li )
{
// this turns out to cluttter up other code..  
li.push(); // shift commands and remove invoking one, 
// use CommandInterpretterParam 
//const StrTy & cmd=li.word(0);

//if (cmd=="int") { test_int(li);} //   continue; } 
//else if (cmd=="diff") { test_diff(li);} //   continue; } 
//else if (cmd=="fl") { test_fl(li);} //   continue; } 
//void test_gr_integrator( CommandInterpretter & li )
if (false) {}else { MM_ERR(" unrecignized TEST command "<<li.dump()) } 

li.pop();
} // test

void dump_cm()
{
 m_cm.dump("solve_step"  , std::cout);
 MM_MSG("solve_times"<<CRLF<<m_cm.time_table("solver_times"))
}

void config_banner()
{
MM_INC_MSG(m_cm,"test test" )
MM_MSG(" configuration "<<m_flp.to_string())
//MM_MSG(" logic "<<m_tree.to_string())
//MM_MSG(" membeers "<<MMPR(m_ord_col)<<MMPR(m_t_col))
}
bool done() const  { return m_done; } 

void  dump_dmel(OsTy & os )  const
{

if (m_dmel==0) { os<<" m_dmel Data Model Error Log is NULL"<<CRLF; return; }
os<<(*m_dmel).string()<<CRLF;
}

//////////////////////////////////////////////////////////////////////////
// end of skeleton, now for the meat 
/////////////////////////////////////////////////////////////////////////
// need to move and use better impl faster or more general

IdxTy Sendt(const StrTy& s, const IdxTy flags)
{
IdxTy rc=0;
BaseParams sk(s);
IdxTy hbeat=0;
StrTy src="-";
StrTy chunks="";
StrTy chars="";
sk.get(hbeat,"hbeat");
sk.get(src,"src");
sk.get(chars,"chars");
sk.get(chunks,"chunks");
// right now TODO FIXME nexted strings not suppored
// or could pass chunks 
if (chunks!="") return SignalChunks(s,flags);
if (hbeat!=0) return Hbeat(s,flags);
std::istream * pifs=0;
try { 
if (chars.length()){ pifs= new std::stringstream(chars); }
else if (src=="-") pifs=&std::cin;
else pifs= new std::ifstream(src);
rc=Send(*pifs,s,flags);
} catch (...)
{
if (pifs!=&std::cin) delete pifs;
MM_ERR(" throwing  "<<MMPR(__FUNCTION__))
throw; 
}

if (pifs!=&std::cin) delete pifs;
MM_ERR(" returning "<<MMPR(__FUNCTION__))
return rc; 
} // Sendt

IdxTy Sendtfifo(const StrTy& s, const IdxTy flags)
{
IdxTy rc=0;
BaseParams sk(s);
IdxTy hbeat=0;
StrTy src="-";
StrTy chunks="";
StrTy chars="";
sk.get(hbeat,"hbeat");
sk.get(src,"src");
sk.get(chars,"chars");
sk.get(chunks,"chunks");
// right now TODO FIXME nexted strings not suppored
// or could pass chunks 
if (chunks!="") return SignalChunksfifo(s,flags);
if (hbeat!=0) return Hbeat(s,flags);
std::istream * pifs=0;
try { 
if (chars.length()){ pifs= new std::stringstream(chars); }
else if (src=="-") pifs=&std::cin;
else pifs= new std::ifstream(src);
rc=Send(*pifs,s,flags);
} catch (...)
{
if (pifs!=&std::cin) delete pifs;
MM_ERR(" throwing  "<<MMPR(__FUNCTION__))
throw; 
}

if (pifs!=&std::cin) delete pifs;
MM_ERR(" returning "<<MMPR(__FUNCTION__))
return rc; 
} // Sendtfifo


#ifndef NO_DSCOPE_GRAPHICS

IdxTy Displaytfifo(const StrTy& s, const IdxTy flags)
{
IdxTy rc=0;
bool debug_fifo=(m_debug!=0);
MM_ERR(MMPR(__FUNCTION__))
BaseParams sk(s);
IdxTy hbeat=0;
StrTy src="-";
StrTy name="name";
StrTy chars="";
sk.get(hbeat,"hbeat");
sk.get(src,"src");
// this is the name of a single channel not a source
sk.get(name,"name");
sk.get(chars,"chars");
sk.get(debug_fifo,"debug_fifo");
//RagScene rs;
while (true)
{
//MM_ERR(MMPR(__FUNCTION__))
Ragged r;
// some one has to sort these out... 
// TODO FIXME needs to poll of have other menas to check pktfifo
// or other sources doh 
m_rawfifo.get(r,flags);
AddNewRagged(r,name,flags);
#if 0
//MM_ERR(MMPR(r.dump_ssv()))
if (false) MM_ERR(" dgram received "<<MMPR(r.size()))
if (r.size() == 0 ) continue;
// Line  k=r[0];
// this finds one with same r[0] line or creates it. 
RagScene * rs=  m_scenes[r];
if (rs==0)
{
MM_ERR(" no p ")
return 0;
} 
if (debug_fifo)  MM_ERR(MMPR(rs->dump()))
// add this data rag to the scene. 
 rs->add(r,0);

// this needs to create a glut_rag 
// 
//rs.add(r,0);
// the use of name here is channel specific not source
// specific as of properly locked the channels are source
// muxed. Name is also in the ragged which reflects the source.
// Multiple lesteners/server ports may exist too
// but generally display does not care... 
// infact, only one port of each type is allowed in the app right
// now 
// this is dumb it then allows only a source name for this scene
m_glut.add(rs,name,flags);
#endif

}
MM_ERR(MMPR(__FUNCTION__))

return rc; 
} // Displaytfifo



// "name" is the name of the channel ( fifo, socket et )
// but each channel is multiplexed and needs to get id
// for adding to dscope 
IdxTy AddNewRagged(const Ragged & r,const StrTy & name,const IdxTy flags)
{

//MM_ERR(MMPR(r.dump_ssv()))
if (false) MM_ERR(name<< "  received "<<MMPR(r.size()))
if (r.size() == 0 ) return 0; //  continue;
//const Line & k=r[0];
// this needs to get a string key from r[0][len-1] 
// or last entry on line 0 
RagScene * rs=m_scenes[r];
rs->pool(&m_pool);
const StrTy srcname=rs->srcname();
if (false) MM_ERR(MMPR(rs->dump()))
// fuck this needs a model info here ... 
rs->add(r,0);
// this needs to create a glut_rag 
// 
//rs.add(r,0);

m_glut.add(rs,srcname,flags);

return 0; 
} // AddNewRagged


IdxTy Displayt(const StrTy& s, const IdxTy flags)
{
IdxTy rc=0;
MM_ERR(MMPR(__FUNCTION__))
BaseParams sk(s);
IdxTy hbeat=0;
StrTy src="-";
StrTy name="name";
StrTy chars="";
sk.get(hbeat,"hbeat");
sk.get(src,"src");
sk.get(name,"name");
sk.get(chars,"chars");
//RagScene rs;
while (true)
{
//MM_ERR(MMPR(__FUNCTION__))
Ragged r;
// some one has to sort these out... 
// TODO FIXME needs to poll of have other menas to check pktfifo
// or other sources doh 
m_dgrams.get(r,flags);
//src name 
AddNewRagged(r,name,flags);

#if 0 
//MM_ERR(MMPR(r.dump_ssv()))
if (false) MM_ERR(" dgram received "<<MMPR(r.size()))
if (r.size() == 0 ) continue;
//const Line & k=r[0];
RagScene * rs=m_scenes[r];
MM_ERR(MMPR(rs->dump()))

rs->add(r,0);
// this needs to create a glut_rag 
// 
//rs.add(r,0);

m_glut.add(rs,name,flags);

#endif

}


MM_ERR(MMPR(__FUNCTION__))

return rc; 
} // Displayt 
#endif

IdxTy SignalChunksfifo(const StrTy s, const IdxTy flags)
{ return SignalChunks(s,flags|(1<<8)); } 
IdxTy SignalChunks(const StrTy s, const IdxTy flags)
{
const bool rawf=Bit(flags,8);
/// TODO FIXME and AGAIN reparse s into sk... doh 
IdxTy rc=0;
IdxTy sleep=1000000;
BaseParams sk(s);
MM_ERR(MMPR2(__FUNCTION__,rawf))
//DgramIO m_dgrams;
IdxTy line=0;
D time=0;
D dt=1;
StrTy idn="";
sk.get(idn,"idn");
while (!m_stop_send)
{
Ragged r; 
{ Line l; l.push_back("#"); l.push_back("id"); l.push_back("test pattern chunks "); l.push_back(idn);  r.add(l); } 
{ Line l; l.push_back("#"); l.push_back("oscilloscope"); l.push_back("chunks"); r.add(l); } 

{ Line l; l.push_back("#"); l.push_back("trace"); Ss ss; ss<<line;  l.push_back(ss.str()); r.add(l); } 
for(IdxTy i=0; i<200; ++i)
{
Line l;
D y=sin(time*.1);
{ Ss ss; ss<<time; l.push_back(ss.str()); }  
{ Ss ss; ss<<y; l.push_back(ss.str()); }  
r.add(l);
time=time+dt;
} // i 
MM_ERR(MMPR4(line,__FUNCTION__,rawf,r.size()))
if (rawf) m_rawfifo.send(r);
else m_dgrams.send(r);
usleep(sleep);
++line;
} // while
MM_ERR(MMPR(__FUNCTION__))
return rc; 
} //SignalChunks 



IdxTy Hbeat(const StrTy s, const IdxTy flags)
{
/// TODO FIXME and AGAIN reparse s into sk... doh 
IdxTy rc=0;
IdxTy sleep=1000000;
BaseParams sk(s);
MM_ERR(MMPR(__FUNCTION__))
//DgramIO m_dgrams;
IdxTy line=0;
while (!m_stop_send)
{
Ragged r; 
Line l;
Ss ss;
IdxTy sht=0;
ss<<MMPR(line);
l.push_back(ss.str());
l.push_back(" vulgar terms here  ");
{ Ss rr; rr<<sht; ++sht;  l.push_back(rr.str()); } 
r.add(l);
{ Ss rr; rr<<sht; ++sht;  l.push_back(rr.str()); } 
r.add(l);
{ Ss rr; rr<<sht; ++sht;  l.push_back(rr.str()); } 
r.add(l);
{ Ss rr; rr<<sht; ++sht;  l.push_back(rr.str()); } 
r.add(l);
MM_ERR(" sending ragged "<<MMPR3(line, r.size(),r.dump_ssv()))
m_dgrams.send(r);
usleep(sleep);
++line;
} // while
MM_ERR(MMPR(__FUNCTION__))
return rc; 
} // Hbeat

IdxTy Send(std::istream & is, const StrTy s, const IdxTy flags)
{
BaseParams sk(s);
IdxTy bufsz=1<<16;
IdxTy bufrd=bufsz-2;
sk.get(bufsz,"bufsz");
sk.get(bufrd,"bufrd");
MM_ERR(MMPR4(__FUNCTION__,bufsz,bufrd,flags))
IdxTy line=0;
IdxTy bads=0;
// get this from the dgram allocator? 
Data * buf= new Data[bufsz];
memset(buf,0,bufsz);
try{
while (!is.eof()&&!m_stop_send)
{
// sending one line at a time, dumb need timers erc
is.getline((char*)(buf),bufrd);
MM_ERR(" have input line "<<MMPR2(is.bad(),buf))
if (is.bad())
{
MM_ERR(MMPR2(bads,is.bad()))
is.clear();
++bads;
continue;
} // is.bad
Ragged r;
Ss ss((const char *)buf);
r.load(ss);
m_dgrams.send(r);
MM_ERR("sent from datascopr send")

++line;
} // wihle 
} catch (...)
{
if (bufsz>11) buf[10]=0; else 
if (bufsz>0) buf[bufsz-1]=0;
MM_ERR(" sending threw "<<MMPR4(__FUNCTION__,line,bufrd,s)<<MMPR(bads))
delete [] buf;
throw ;
}

MM_ERR(" sending ends "<<MMPR3(__FUNCTION__,line,is.eof()))
delete [] buf;
return 0;
} // Send

//////////////////////////////////////////////////
void parse_dmel( const StrTy & d, const StrTy & n, const IdxTy start,  const CommandInterpretter & li)
{
	const IdxTy sz=li.size();
	const Words & w=li.words();
	if (sz<=start) return;
	const StrTy & name=w[start];
// do nothing with this for now 
		if (m_dmel!=0) 
		{ (*m_dmel).event("userdmelentry" ,name, d,w); } 

}




/////////////////////////////////////////////////////////////////////
private:
void Init()
{
m_debug=0;
#ifndef NO_DSCOPE_GRAPHICS
GlutScope::p(& m_glut);
m_pool.set_size(100);
m_glut.app(this);
#endif
m_done=false;
m_stop_send=false;
}
bool Bit(const FlagTy flags, const IdxTy b) { return (flags&(1<<b))!=0; }

/*
void DMel(const StrTy & e)
{


}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<" "<<code<<CRLF;
    ss.str(StrTy(""));
}

*/

void DMel(const StrTy & e)
{
MM_ERR(e)
if (m_dmel!=0) {m_dmel->event("wtf",e); }
}
//void DMel(Ss & ss)
void DMel(OsTy & _ss)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    std::cerr<<ss.str()<<CRLF;
	if (m_dmel!=0) {m_dmel->event("wtf",ss.str()); }
    ss.str(StrTy(""));
}
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; }
	if (m_dmel!=0) {m_dmel->event(code,ss.str()); }
    ss.str(StrTy(""));
}
void DumpDMel() // (OsTy & os)
{
if (m_dmel!=0) { MM_ERR(MMPR((m_dmel->string(1)))) } 
else { MM_ERR(" m_dmel is null ") }

}




// MEMBERS members
IdxTy m_debug;
#ifndef NO_DSCOPE_GRAPHICS
GlutScope m_glut;
// TODO the pool belongs in the SceneMap maybe?
SceneMap m_scenes;
ModelPool m_pool;
Saver m_saver;
#else
public:
#endif

DgramID m_dgrams;
PktFifoID m_pktfifo;
RawFifoID m_rawfifo;
bool m_done,m_stop_send;
ParamGlob m_flp;
VariableStore m_variables;
Dmel * m_dmel;
//Ss m_ss;
//FastaMap m_fasta_map;
RaggedMap m_ragged_map;
TokRaggedMap m_tokragged_map;
//DigMap m_dig_map;
//PhenoMap m_pheno_map;
//HbMap m_hbpheno_map;
//TestStringMap m_queries;
//CharMat m_char_mat;
//MiiMap m_luts;
//TaxTree m_tax_tree; // now have sequences ID's to taxon 
//TaxTrees m_tax_trees;
CounterMap m_cm;
CliTy m_cli;
}; //mjm_datascope 

/////////////////////////////////////////////////////////

#ifdef  TEST_datascope__
int main(int argc,char **args)
{
typedef mjm_datascope Myt;
///typedef double D;
//typedef unsigned int IdxTy;
Myt x(argc,args);
if (!x.done()) x.command_mode();

return 0;
}

#endif


#endif

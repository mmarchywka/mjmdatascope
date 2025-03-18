#ifndef MJM_CLI_UI_H__
#define MJM_CLI_UI_H__

#include "mjm_globals.h"
#include "mjm_data_model_error_log.h"
//#include "mjm_generic_iterators.h"
//#include "mjm_instruments.h"
//#include "mjm_logic_base.h"
//#include "mjm_calendar.h"
#include "mjm_strings.h"

#include <stdio.h>
#ifndef CLI_EXCLUDE_NCURSES
#include <ncurses.h>
#endif

#include <readline/readline.h>
#include <readline/history.h>

#include <algorithm>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>

// if (s.length()<4) { DMel(m_ss<<MM_STR_LOC<<MMPR2(level,s),false); ++ii;
#define MM_DMEL(code,x) DMel(code, m_ss<<MM_STR_LOC<<" "<<x); 
#define MM_DMELF(fm,code,x) f.DMel(code, f.m_ss<<MM_STR_LOC<<" "<<x); 


/*

http://web.mit.edu/gnu/doc/html/rlman_2.html

// https://stackoverflow.com/questions/13168833/c-command-line-interface

*/


namespace mjm_cli_ui_stuff 
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
//typedef mjm_sparse_matrix<D> MySparse;
}; // 

}; // mjm_cli_ui_stuff
#ifndef CLI_EXCLUDE_NCURSES
template <  class Traits=mjm_cli_ui_stuff::Tr > class mjm_ncurses_links
{
typedef mjm_ncurses_links<Traits> Myt;
//typedef mjm_cli_ui_stuff::Tr Tr; 
typedef Traits Tr; 
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

public:

class layout_result { public: layout_result(const void * p, const IdxTy _x, const IdxTy _y)
: item(p),x(_x),y(_y) {} const void * item; IdxTy x,y; };

// http://tldp.org/HOWTO/NCURSES-Programming-HOWTO/printw.html
void box(const StrTy & s, const IdxTy x0, const IdxTy y0, const IdxTy w, const IdxTy h)
{
const IdxTy a=w*h;
const IdxTy sz=s.length();
const IdxTy xf=x0+w;
const IdxTy yf=y0+h;
for (IdxTy i=0; i<w; ++i)
for (IdxTy j=0; j<h; ++j)
{
const IdxTy k=j*w+i;
if ( k>=sz) break; 
const IdxTy x= xf-i-1;
const IdxTy y= yf-j-1;
const IdxTy p=sz-k-1;
 mvprintw(y,x,"%c",s.c_str()[p]);
}

} // box
void start() { initscr(); raw(); keypad(stdscr,TRUE);noecho();  }
void finish() { endwin(); }
// TODO returning char was not flagged with -Wall 
int  show() { refresh(); return getch(); } 
void clear() { ::clear(); } 
typedef std::vector<layout_result> ScreenLocations;
void branch(const StrTy & base, const std::vector<StrTy> & kids,const IdxTy firstchild,  ScreenLocations * sl=0)
{
 int row,col;			
 getmaxyx(stdscr,row,col);		/* get the number of rows and columns */
const IdxTy baselen=base.length();
IdxTy longest=0;
MM_LOOP(ii,kids)
{
const IdxTy len=(*ii).length();

if (len>longest) longest=len;  
}
IdxTy ckids=baselen+1;
IdxTy extra=0,eb=0;
IdxTy coldiv=col>>1;
if ( longest<coldiv) extra=coldiv-longest;
if ( baselen<coldiv) eb=coldiv-baselen;
IdxTy y=row>>1;
IdxTy x=0;
IdxTy ww=(col>>1)-eb;
IdxTy hh=row>>2;
box(base,x,y,ww,hh);
if (sl!=0) { sl->push_back(layout_result(&base,x+ww-1,y+hh-1)); } 
IdxTy i=0;
MM_LOOP(ii,kids)
{
if ( i<firstchild){ ++i;  continue; } 
x=ckids;
y=i-firstchild;
if (y>=IdxTy(row)) break;
ww=col-ckids-extra;
hh=1;
box((*ii),x,y,ww,hh);
if (sl!=0) { sl->push_back(layout_result(&(*ii),x+ww-1,y+hh-1)); } 

++i;
}

}

 
int main()
{
 char mesg[]="Just a string";		/* message to be appeared on the screen */
 int row,col;				/* to store the number of rows and *
					 * the number of colums of the screen */
 initscr();				/* start the curses mode */
 getmaxyx(stdscr,row,col);		/* get the number of rows and columns */
 mvprintw(row/2,(col-strlen(mesg))/2,"%s",mesg);
                                	/* print the message at the center of the screen */
 mvprintw(row-2,0,"This screen has %d rows and %d columns\n",row,col);
 printw("Try resizing your window(if possible) and then run this program again");
 refresh();
 getch();
 endwin();

 return 0;
}
private:

}; // mjm_ncurses_links

#endif

template < class Tsrc, class Traits=mjm_cli_ui_stuff::Tr >
class mjm_cli_ui 
{
typedef mjm_cli_ui<Tsrc,Traits>  Myt;
protected:
typedef Traits Tr;
//typedef Tobj To;
typedef typename Tr::IdxTy IdxTy;
typedef typename Tr::D D;
typedef typename Tr::Ss Ss;
typedef typename Tr::IsTy IsTy;
typedef typename Tr::OsTy OsTy;
typedef typename Tr::MyBlock  MyBlock;
typedef typename Tr::Dmel  Dmel;

//typedef string_tokenizer St;

//typedef typename To::const_iterator CItor;
//typedef typename To::iterator Itor;



//http://web.mit.edu/gnu/doc/html/rlman_2.html
typedef void VFunction ();

//static char * command_generator (char * text, int state)//  char *text; int state;
typedef char *( CPFunction) (const char *, int );
typedef char **( CPPFunction) (const char*,int,int);

typedef std::vector<StrTy> ListTy;

typedef void (Tsrc::* TargCmd)( ListTy & choices,  const char * frag);
typedef void (Tsrc::* TargParam)( ListTy & choices, const char *  cmd, const char * frag);

StrTy m_app; // add name for .inputrc 
bool m_set_app,m_active;

Tsrc * m_targ;
TargCmd  m_targ_cmd;
TargParam  m_targ_param;
CPPFunction * m_old_complete;
ListTy m_commands;

static Myt *&  current()
{
static Myt * p=0;
return p; 
}

public:
mjm_cli_ui(): m_set_app(false),m_active(false),m_targ(0),m_targ_cmd(0),m_targ_param(0),m_dmel(0)  {}
mjm_cli_ui(Tsrc & src): m_set_app(false),m_active(false),m_targ(&src),m_targ_cmd(0),m_targ_param(0),m_dmel(0)  {}
virtual ~mjm_cli_ui() {}
typedef ListTy list_type;
void set_target(Tsrc &  p) { m_targ=&p; } 
void set_command_handler(TargCmd  p ) { m_targ_cmd=p; } 
void set_param_handler(TargParam  p ) { m_targ_param=p; } 
void activate()
{
if (m_active) return;
current()=(this);
initialize_readline();
m_active=true;
}
void deactivate()
{
if (!m_active) return;
rl_attempted_completion_function=m_old_complete;
current()=0; // this needs to pop FIXME TODO 
m_active=true;

}
void set_list(list_type & list) { m_commands=list; } 

//http://web.mit.edu/gnu/doc/html/rlman_2.html
void initialize_readline ()
{
  /* Allow conditional parsing of the ~/.inputrc file. */
  rl_readline_name = m_set_app?m_app.c_str():0; // "FileMan";
  /* Tell the completer that we want a crack first. */
  m_old_complete=rl_attempted_completion_function;
  rl_attempted_completion_function = (CPPFunction *)Myt::fileman_completion;
}

static char **
fileman_completion (char * text, int start, int end) // char *text; int start, end;
{
  char **matches;

//MM_ERR(" completion  "<<MMPR4(text,rl_line_buffer,start,end))
  matches = (char **)NULL;

  /* If this word is at the start of the line, then it is a command
     to complete.  Otherwise it is the name of a file in the current
     directory. */
//  if (start == 0)
    matches = rl_completion_matches (text, Myt::command_generator);

  return (matches);
}


/* Generator function for command completion.  STATE lets us know whether
   to start from scratch; without any state (i.e. STATE == 0), then we
   start at the top of the list. */
static char *
command_generator (const char * text, int state)//  char *text; int state;
{
//MM_ERR(" cmd gen "<<MMPR2(text,state))

  /* If this is a new word to complete, initialize now.  This includes
     saving the length of TEXT for efficiency, and initializing the index
     variable to 0. */
  //    list_index = 0;
  if (state!=0) { return ((char *)NULL); } 
const IdxTy len=strlen(text);
  /* Return the next name which partially matches from the command list. */
// TODO FIXME binary search iterator for stings arrays 
  //while (name = commands[list_index].name)
ListTy * commands=&(current()->m_commands); // &m_commands;
  //if(commands==0) return ((char *)NULL);
 //if (state<10)  if(commands==0) return ("adfasdfasfasdfa");
  if(commands==0) return ((char *)NULL);
	if ((current()->m_targ)!=0) 
{

if (state==0)
{
TargCmd  fudd=current()->m_targ_cmd;
Tsrc * shot=current()->m_targ;
 (shot->*fudd)((*commands), text); 
} else
{
TargParam  fudd=current()->m_targ_param;
Tsrc * shot=current()->m_targ;
 (shot->*fudd)((*commands), rl_line_buffer, text); 


} 

}
//  while (0!=(name = (*commands)[list_index].c_str()))
	MM_LOOP(ii,(*commands))
    {
		const char * name=(*ii).c_str();
//      list_index++;
// the url code calls xmalloc, I guess new is ok ... 
      if (strncmp (name, text, len) == 0)
        //return (dupstr(name));
		{
			const IdxTy nlen=strlen(name);
			char * x = new char[nlen+1];
			strcpy(x,name);
			x[nlen]=0;
			//MM_ERR(" returning "<<MMPR(x))
			return x; 
		} 
    }

  /* If no names matched, then return NULL. */
  return ((char *)NULL);
}




static StrTy trim(const StrTy & s)
{
Ss ss;
ss<<s;
StrTy x;
while ( true) { StrTy y; ss>>y; if ( y=="") break; if (x!="") x=x+StrTy(" "); x=x+y; }
return x;
}
static StrTy word(const StrTy & s, const IdxTy n)
{
Ss ss;
ss<<s;
StrTy x;
IdxTy i=0; 
while ( true) { StrTy y; ss>>y; if ( y=="") break; if (i==n) return y;  ++i;  }
return x;
}





private:

int myatoi(const StrTy & s ) const { return myatoi(s.c_str()); }
int myatoi(const char * c) const { return ::strtol(c,0,0); }

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
void DMel(const StrTy & code, OsTy & _ss, const bool print=true)
{
    Ss& ss=dynamic_cast<Ss& >(_ss);
    if ( print ) { std::cerr<<ss.str()<<" "<<code<<CRLF; } 
    ss.str(StrTy(""));
}


Dmel * m_dmel;
// 2023-12-25 complaining about copy ctor on assignment? lol
//Ss m_ss;

};  // mjm_collection

#undef MM_DMEL

#endif


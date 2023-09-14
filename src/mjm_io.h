
#ifndef MJM_IO_H__
#define MJM_IO_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>


#include "mjm_globals.h"

class mjm_io 
{
private:
typedef mjm_io Myt;
typedef std::runtime_error ErTy;
typedef std::runtime_error CmdEr;

public:
typedef std::ostream OsTy;
typedef std::istream IsTy;
typedef std::string StrTy;
typedef unsigned int FlagTy;
typedef unsigned int IdxTy;
OsTy * m_ss[4];



~mjm_io()
{
// don't want to really risk killing
// on bad files, 
//std::cout<<MM_MARK<<" dtor called at "<<" " <<" lines/sec="<<CRLF; 
// counting on dtor to flush files and close
/*
if( m_ss[0]!= &std::cout) delete m_ss[0];
if( m_ss[1]!= &std::cerr) delete m_ss[1];
if( m_ss[2]!= &std::cout) delete m_ss[2];
if( m_ss[3]!= &std::cerr) delete m_ss[3];
*/

}
void setOs(OsTy * os,const IdxTy i ) { m_ss[i]=os; } 

OsTy & os(const IdxTy i) { return *(m_ss[i]); } 

static void stream_diagnose(const StrTy& lbl, std::ios & s, OsTy & os)
{
os<<MM_MARK<<" "<<lbl<<" good="<<s.good()<<" fail="<<s.fail()<< 
" eof="<<s.eof()<<" rdstate="<<std::hex<<s.rdstate()<<std::dec<<CRLF;


}

static IsTy * get_in(const char * f, OsTy & os )
{

	IsTy * is = new std::ifstream(f,std::ios::binary);
	if ( !(*is).good())
	{
		os<<MM_MARK<<" no luck open "<<f<<CRLF;
		stream_diagnose(StrTy(f),*is,os);
		return NULL;
	}
return is;
}

static IsTy * get_in(const char * f )
{

	IsTy * is = new std::ifstream(f,std::ios::binary);
	if ( !(*is).good())
	{
		return NULL;
	}
return is;
}








static OsTy * get_out(const char * f, OsTy & os )
{

	OsTy * is = new std::ofstream(f,std::ios::binary);
	if ( !(*is).good())
	{
		os<<MM_MARK<<" no luck open "<<f<<CRLF;
		stream_diagnose(StrTy(f),*is,os);
		return NULL;
	}
return is;
}


/*
IdxTy vatou(const char * g, const IdxTy m )
{

	bool ok=true;
	if ( ::strlen(g)==0) ok=false;
	if (ok)
	{	const IdxTy lj=::strlen(g);
		for ( IdxTy j=0; j<lj; ++j)
			if ( !::isdigit(g[j])) { ok=false ; break; }
	}
	if ( !ok)
	{ 
	esx()<<MM_MARK<<" bad os number "<<g<<" for output file "<<" "<< CRLF;
	throw CmdEr(ext_app(" file spec not good ",g));
	}
// this relies on unsigned being >=0
	IdxTy idx= ::atoi(g);
if ( m==0) return idx;
	if (idx>m) 
	{
		esx()<<MM_MARK<<" bad os number "<<g<<" parsed to "<<idx
		<<" limit is "<<m<<CRLF;
		throw CmdEr(ext_app(" bad parse of input number  ",g));
	}
return idx;

}

*/

/*
void os_cmd(int argc, char ** argv, int & i, MyClass & myc)
{
	const char * g=argv[++i];
	const char * f=argv[++i];
	OsTy * is = get_out(f,myc.os(1));
	osx()<<MM_MARK<<" now to open "<<g<<" "<<f<< CRLF; 	
	if (is==NULL)
	{
		myc.os(1)<<MM_MARK<<" can't open os file "<<f<<CRLF;
		throw CmdEr(ext_app(" faile to open", f)); 
	}
		
	IdxTy idx=vatou(g,3);
	myc.setOs(is,idx);
// (*is)<<MM_MARK<<"test write"<<CRLF;	
	if ( PEDANTIC_AND_CONDESCENDING) 
 		stream_diagnose("cmd line file opn",*is,osx());
	osx()<<MM_MARK<<" seemed to set output file ok "
			<<idx<<" for "<<g<<" to "<<f<<CRLF;

}
*/

};  // mjm_io


#endif


#ifndef MJM_DMEL_H__
#define MJM_DMEL_H__

#include "mjm_globals.h"
#include "mjm_generic_iterators.h"
#include "mjm_block_matrix.h"
#include "mjm_instruments.h"
#include "mjm_logic_base.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>

/*
*/

// record exceptions and warnings for data or code problems diagnosis 
class data_model_error_log 
{
typedef std::string StrTy;
typedef std::stringstream Ss;
typedef double D;
typedef unsigned int IdxTy;
typedef std::ostream Os;
//typedef std::map<StrTy, D > DoseMap;
typedef std::vector<StrTy > Words;
class entry
{
public:
StrTy dump(const IdxTy flags=0) const
{
Ss ss;
switch (flags)
{
case 0:{ ss<<MMPR3(word,date,problem); break; }
case 1:{ ss<<MMPR2(word,problem); break; }
default:
ss<<MM_STR_LOC<<" bad case "<<MMPR(flags);
} // switch

for (IdxTy i=0; i<adjectives.size(); ++i) ss<<" "<<adjectives[i];
return ss.str();
}
StrTy word;
Words adjectives;
StrTy date;
StrTy problem;
 
}; // entry 
typedef entry Event;
typedef std::vector<Event > EventList;
typedef std::map<StrTy, IdxTy> ProblemCount;
public:
void event(const StrTy & p , const StrTy & word)
{ Event x; x.problem=p; x.word=word; event(x); }
void event(const StrTy & p , const StrTy & word, const StrTy & d)
{ Event x; x.problem=p; x.word=word; x.date=d; event(x); }
void event(const StrTy & p , const StrTy & word, const StrTy & d,const Words & words)
{ Event x; x.problem=p; x.word=word; x.date=d; x.adjectives=words; event(x); }

void event(const Event & x)
{
	m_events.push_back(x);
	++m_problems[x.problem];
}

StrTy string(const IdxTy flags=0) const
{
Ss ss;
const IdxTy sz=m_events.size();
for (IdxTy i=0; i<sz; ++i)
{
ss<<i<<" "<<m_events[i].dump(flags)<<CRLF;
}

return ss.str();

}

EventList m_events;
ProblemCount m_problems;
}; // data_model_error_log

#endif



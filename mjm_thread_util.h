#ifndef MJM_THREAD_UTIL_H__
#define MJM_THREAD_UTIL_H__

#include <pthread.h>
#include <vector>
// https://stackoverflow.com/questions/21091000/how-to-get-thread-id-of-a-pthread-in-linux-c-program

#include <unistd.h>
#include <sys/types.h>


template <class Tr>
class mjm_thread_util 
{
 typedef mjm_thread_util Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
// typedef typename Tr::MyBlock  MyBlock;

typedef pthread_mutex_t Mutex;

public:

class mutex
{
Mutex  m_mutex;
Mutex * pmutex() { return &m_mutex; }

public:
void enter_serial() { pthread_mutex_lock( pmutex() ); }
void exit_serial() { pthread_mutex_unlock( pmutex() ); }

mutex ()
{
m_mutex = PTHREAD_MUTEX_INITIALIZER;
// m_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
}
};  // mutex


class reentrant_mutex
{
Mutex  m_mutex;
Mutex * pmutex() { return &m_mutex; }

public:
void enter_serial() { pthread_mutex_lock( pmutex() ); }
void exit_serial() { pthread_mutex_unlock( pmutex() ); }

reentrant_mutex ()
{
 //m_mutex = PTHREAD_MUTEX_INITIALIZER;
 m_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
}
};  // reentrant_mutex



class mutex_vector
{

typedef std::vector<mutex> MuVector;
public:
mutex_vector() : m_v() {}
mutex_vector(const IdxTy n) : m_v(n) {}
IdxTy size() const { return m_v.size(); } 
void enter_serial(const IdxTy n) { m_v[n].enter_serial(); } 
void exit_serial(const IdxTy n) { m_v[n].exit_serial(); } 
typename MuVector::iterator begin() { return m_v.begin(); } 
typename MuVector::iterator end() { return m_v.end(); } 
private:

MuVector m_v;

}; // mutex_vector

typedef mutex_vector MutexVector; 

class reentrant_mutex_vector
{

typedef std::vector<mutex> MuVector;
public:
reentrant_mutex_vector() : m_v() {}
reentrant_mutex_vector(const IdxTy n) : m_v(n) {}
void enter_serial(const IdxTy n) { m_v[n].enter_serial(); } 
void exit_serial(const IdxTy n) { m_v[n].exit_serial(); } 

private:

MuVector m_v;

}; // reentrant_mutex_vector
typedef reentrant_mutex_vector ReentrantMutexVector; 

class wait_notify
{

public:
wait_notify() { Init(); } 
void wait()
{
      pthread_mutex_lock( &condition_mutex );
      pthread_cond_wait( &condition_cond, &condition_mutex );
      pthread_mutex_unlock( &condition_mutex );
} // wait
void notify()
{
       pthread_mutex_lock( &condition_mutex );
          pthread_cond_signal( &condition_cond );
       pthread_mutex_unlock( &condition_mutex );
} // notify 
private:
void Init()
{
 condition_mutex = PTHREAD_MUTEX_INITIALIZER;
 condition_cond  = PTHREAD_COND_INITIALIZER;
} 

pthread_mutex_t condition_mutex; //  = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  condition_cond ; //  = PTHREAD_COND_INITIALIZER;
}; // wait_notify
/*
pthread_mutex_t count_mutex     = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t condition_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t  condition_cond  = PTHREAD_COND_INITIALIZER;

void *functionCount1();
void *functionCount2();
int  count = 0;
#define COUNT_DONE  10
#define COUNT_HALT1  3
#define COUNT_HALT2  6

main()
{  
   pthread_t thread1, thread2;

   pthread_create( &thread1, NULL, &functionCount1, NULL);
   pthread_create( &thread2, NULL, &functionCount2, NULL);
   pthread_join( thread1, NULL);
   pthread_join( thread2, NULL);

   exit(0);
}

void *functionCount1()
{
   for(;;)
   {
      pthread_mutex_lock( &condition_mutex );
      while( count >= COUNT_HALT1 && count <= COUNT_HALT2 )
      {
         pthread_cond_wait( &condition_cond, &condition_mutex );
      }
      pthread_mutex_unlock( &condition_mutex );




*/





// this would be rare as it is static to the base class
/*template <IdxTy MU_SZ>
class static_locks
{
public:
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
}; // static_locks
*/

class locks {
protected:
locks(const IdxTy n)  { Init(n); } 
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }

void Init(const IdxTy n )
{
m_mutex_vector = MutexVector(n);
}
private:
 mutable MutexVector m_mutex_vector;

}; // locks

class reentrant_locks {
typedef ReentrantMutexVector Mv;
protected:
reentrant_locks(const IdxTy n)  { Init(n); } 
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }

void Init(const IdxTy n )
{
m_mutex_vector = Mv(n);
}
private:
 mutable Mv m_mutex_vector;

}; // locks


template<class Ty> class ThParam
{
// this is passed as a new'ed void* and needs to be
// deleted by recipient.. 
public:
typedef IdxTy (Ty::*  pFunc)(const StrTy & s,const IdxTy flags );
ThParam() :p(0),flags(0),func(0) {}
ThParam(Ty * pt,pFunc pf,const StrTy & _s,const IdxTy flags ) 
:p(pt),s(_s),flags(0),func(pf) {}
pthread_t invoke()
{
pthread_t thread=launch(& ThParam::go,this);
return thread; 
}
//static void *  _Launch( void * x) { ThParam * tp =(ThParam*)x;
//Myt * p = tp->p; StrTy s=tp->s; IdxTy flags=tp->flags;
//auto pFunc=tp->func;
static void*  go(void * x)
{
ThParam * pp=(ThParam*)x; 
pp->start();
return 0;
} // go 
void start()
{
//Myt * p = tp->p; StrTy s=tp->s; IdxTy flags=tp->flags;
Ty * pt = p; StrTy _s=s; IdxTy _flags=flags;
delete this; 
//auto pFunc=tp->func;
//delete tp;
//p->Start(s,flags); 
try {
(pt->*func)(_s,_flags);
} catch (...) { MM_ERR(" throw "<<MMPR(__FUNCTION__)) }
pthread_exit(NULL);
//delete this;
} // start 

Ty * p;
StrTy s;
IdxTy flags;
//typedef IdxTy (Myt:: * pFunc)(const StrTy & s,const IdxTy flags );
pFunc func;

}; // ThParam



static void launch
(const IdxTy nthreads, void *(*thread_function)(void* ) ,void **pv)
{
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   { pthread_create( &thread_id[i], NULL, thread_function, pv[i] ); }
   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}

static pthread_t launch
(void *(*thread_function)(void* ) ,void *pv)
{
   pthread_t thread_id;
   pthread_create( &thread_id, NULL, thread_function, pv ); 
return thread_id;
}



static void fire_and_forget
( const IdxTy nthreads, void *(*thread_function)(void* ), void *msgp)
{
//MM_ERR(" firing")
   //auto pfunc=thread_function; //  Myt::thread_entry_func;
   //IdxTy nthreads=1;
   pthread_t thread_id[nthreads];
   IdxTy i, j;
   for(i=0; i < nthreads; i++)
   //{ pthread_create( &thread_id[i], NULL, pfunc, msgp+i ); }
   { pthread_create( &thread_id[i], NULL, thread_function, msgp ); }
//   for(j=0; j < nthreads; j++) { pthread_join( thread_id[j], NULL); }
}






static IdxTy  thread_id() 
{
pid_t tid = gettid();
return tid;
}











private:


}; // mjm_thread_util

#endif // MJM_THREAD_UTIL_H__ 

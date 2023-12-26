#ifndef MJM_GLUT_SCOPE_II_H__
#define MJM_GLUT_SCOPE_II_H__

#include "mjm_globals.h"
#include "mjm_thread_util.h"

// the include order here is a problem
// TODO find out wth saver does to make render harder to compile
#include "mjm_glut_saver.h"
#include "mjm_svg_render.h"

#include "mjm_glut_helpers.h"
#include "mjm_glut_rags.h"
#include "mjm_collections.h"
#include "mjm_string_kvp.h"
#include "mjm_so_loader.h"
//#include "mjm_glut_elements.h"

// for the rags objects
#include "mjm_object_pool.h"

#include <freeglut.h> 

#include <map> 
#include <vector> 
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <string>
#include <fstream>
#include <signal.h>
#include <stdlib.h>
#include <stdint.h>


// Wed Mar  1 12:24:36 EST 2023
// generated by /home/documents/cpp/scripts/cpputil -classhdr mjm_glut_scope_ii   
// g++  -Wall -std=gnu++11  -DTEST_MJM_GLUT_SCOPE_II -I. -I../../mjm/hlib -I../../mjm/num   -I/usr/include/GL -gdwarf-3 -O0  -x c++ mjm_glut_scope_ii.h  -o mjm_glut_scope_ii.out -lpthread -lreadline -lGL -lglut -lGLU


mjm_global_credits::credit __credit__mjm_glut_scope_ii("mjm_glut_scope_ii"
, "  ");

// a lot of code snippets from 
// glut_example.c
// Stanford University, CS248, Fall 2000
//
// Demonstrates basic use of GLUT toolkit for CS248 video game assignment.
// More GLUT details at http://reality.sgi.com/mjk_asd/spec3/spec3.html
// Here you'll find examples of initialization, basic viewing transformations,
// mouse and keyboard callbacks, menus, some rendering primitives, lighting,
// double buffering, Z buffering, and texturing.
//
// Matt Ginzton -- magi@cs.stanford.edu
#include <glut.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _WIN32
#   include <windows.h>
#else
#   include <sys/time.h>
#endif
//#include "texture.h"
#define VIEWING_DISTANCE_MIN  3.0
#define TEXTURE_ID_CUBE 1




template <class Tr>
class mjm_glut_scope_ii 
{
 typedef mjm_glut_scope_ii Myt;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;

typedef pthread_t ThreadId;
typedef mjm_thread_util<Tr> ThreadTy;
typedef typename ThreadTy::template ThParam<Myt> ThParam;




typedef GLfloat Gf;

typedef mjm_glut_saver<Tr> Saver;
typedef mjm_glut_helpers<Tr> GlutUtil;
typedef typename  GlutUtil::draw_info_t DrawInfo;
typedef typename  GlutUtil::key_info_t KeyInfo;
typedef typename  GlutUtil::ptr_info_t PtrInfo;

// not used here just for pooling 
typedef typename  GlutUtil::junk_bin_t ModelInfo;
typedef mjm_object_pool<Tr,ModelInfo> ModelPool;


//typedef typename  GlutUtil::draw_info_t scope_draw_param_type;
//typedef typename  GlutUtil::key_info_t scope_key_param_type;
//typedef typename  GlutUtil::ptr_info_t scope_ptr_param_type;



typedef mjm_string_kvp<Tr> StrKvp;
typedef mjm_string_base_params<Tr> BaseParams;
typedef mjm_so_loader<Tr> Loader;
typedef mjm_ragged_table Ragged;
typedef Ragged::Line Line ;
// this should be in traits... 
enum {BAD=~0};


//typedef mjm_glut_elements<Tr> Elements;
//typedef  typename Elements::Element Element ;

//typedef std::map<StrTy,Elements> Scenes;


//typedef std::map<StrTy, View> Views;

enum {
  MENU_LIGHTING = 1,
  MENU_POLYMODE,
  MENU_TEXTURING,
  MENU_EXIT
};
typedef int BOOL;
#define TRUE 1
#define FALSE 0
typedef mjm_glut_rags<Tr> RagScene;
typedef std::map<StrTy,RagScene * > RagSceneMap;

struct _mjm_gl_status
{



bool b1down() const { return m_b1_down;}
bool b2down() const { return m_b2_down;}

bool b1down(const bool x) {m_b1_down=x;  return m_b1_down;}
bool b2down(const bool x) {m_b2_down=x;  return m_b2_down;}

bool b1down(const bool ud,const int x, const int y ) {m_b1_down=ud;  
if (ud) { m_b1_down_x=x; m_b1_down_y=y; }
else { m_b1_up_x=x; m_b1_up_y=y; }
last(x,y);
return m_b1_down;}
bool b2down(const bool ud,const int x, const int y ) {m_b2_down=ud;  
if (ud) { m_b2_down_x=x; m_b2_down_y=y; }
else { m_b2_up_x=x; m_b2_up_y=y; }
last(x,y);
return m_b2_down;}
void last(const int x, const int y ) {m_last_x=x; m_last_y=y; }  




void setup()
{
g_bLightingEnabled = TRUE;
g_bFillPolygons = TRUE;
g_bTexture = FALSE;
//g_bButton1Down = FALSE;
m_b1_down=false;
m_b2_down=false;
g_fTeapotAngle = 0.0;
g_fTeapotAngle2 = 0.0;
//g_fViewDistance = 3 * VIEWING_DISTANCE_MIN;
g_nearPlane = 1;
g_farPlane = 1000;
g_Width = 600;                          // Initial window width
g_Height = 600;                         // Initial window height
//g_yClick = 0;
m_b1_down_x=0; m_b1_down_y=0;
m_b1_up_x=0; m_b1_up_y=0;
m_b2_down_x=0; m_b2_down_y=0;
m_b2_up_x=0; m_b2_up_y=0;
m_last_x=0;
m_last_y=0;

g_lightPos =new float[5] { 10, 10, -100, 1 };  // Position of light

} // Setup 

BOOL g_bLightingEnabled ;
BOOL g_bFillPolygons ;
BOOL g_bTexture ;
//BOOL g_bButton1Down ;
bool m_b1_down;
bool m_b2_down;
GLfloat g_fTeapotAngle;
GLfloat g_fTeapotAngle2;
//GLfloat g_fViewDistance;
GLfloat g_nearPlane ;
GLfloat g_farPlane ;
int g_Width ;                          // Initial window width
int g_Height ;                         // Initial window height
//int g_yClick ;
int m_b1_down_x,m_b1_down_y;
int m_b1_up_x,m_b1_up_y;
int m_b2_down_x,m_b2_down_y;
int m_b2_up_x,m_b2_up_y;
int m_last_x,m_last_y;


//float g_lightPos[4] ;  // Position of light
const float* g_lightPos ;  // Position of light
#ifdef _WIN32
DWORD last_idle_time;
#else
struct timeval last_idle_time;
#endif
}; // _mjm_gl_status
typedef _mjm_gl_status MyGLStatus;



//typedef typename DrawInfo::hot_type::layout_type lt;
typedef typename DrawInfo::hot_type GuiLayout;
typedef mjm_svg_render<Tr> SvgRender;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
// API 
typedef typename Saver::save_params_type Sp;
mjm_glut_scope_ii() {Init();}
~mjm_glut_scope_ii() {Free(); }
// TODO FIXME thes is a hazard here but not sure it is worth
// creating a posisble deadlock on panic exit until
// threading worked out.  
void exit_loop() {MM_ERR(MMPR(alive())) 
if (alive()){  glutLeaveMainLoop(); 
// it does nothing until this happens doh 
// this STILL does not fing work f 
// wtf? 
try { 
  glutIdleFunc (Myt::_AnimateScene);
glutPostRedisplay(); } catch (...) {}
 }} 
bool alive() const { return m_alive; } 
// this is a source of confiusion.. 
static Myt * p( Myt * pin=0)  { static Myt * x=(pin)?pin: new Myt(); return x; } 
// add a scene - this creates ptr ownership issues, smartptr?
// this comes from datascope and it has not parse
RagScene *  add( RagScene * p, const StrTy name, const IdxTy flags) 
{ return Add(p,name,flags);  } 
IdxTy actives() const { return m_actives.size();}
IdxTy scenes() const { return m_scenes.size();}

// set active scene
void activate( const StrTy name, const IdxTy flags){Activate(name,flags); } 
void deactivate( const StrTy name, const IdxTy flags){Deactivate(name,flags); } 
void set_saver(Saver * p ) { m_psaver=p; }
void set_saver_params(Sp & sp) { m_save_params=sp; } 
IdxTy launch(const StrTy & s) { return Launch(s);  } 
IdxTy start(const StrTy & s) { return Start(s);  } 
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:

RagScene *  Add( RagScene * p, const StrTy name, const IdxTy flags) 
{ 
EnterSerial(0);
// TODO FIXME need a mutex here.. 
//MM_ERR(" in Add  "<<MMPR2(this,p))
auto ii=m_scenes.find(name);
RagScene* pold=0;
if (ii!=m_scenes.end()) pold=(*ii).second;
m_scenes[name]=p; 
//if (actives()==0) 
if (m_deactive.find(name)==m_deactive.end()) {	activate(name,0); 
if (pold==0) BuildPopupMenu(); } 
glutPostRedisplay();
ExitSerial(0);
return pold; 
} // Add  

void  Activate( const StrTy name, const IdxTy flags) 
{ 
// TODO FIXME need a mutex here.. 
auto ii=m_scenes.find(name);
if (ii!=m_scenes.end())
	if (m_actives.find(name)==m_actives.end())  m_actives[name]=(*ii).second;
{auto ii=m_deactive.find(name);
if (ii!=m_deactive.end()) m_deactive.erase(ii);
}
glutPostRedisplay();
} // Activate  


void  Deactivate( const StrTy name, const IdxTy flags) 
{ 
// TODO FIXME need a mutex here.. 
{auto ii=m_actives.find(name);
if (ii!=m_actives.end()) {  m_actives.erase(ii); m_deactive[name]=(*ii).second; } 

}
glutPostRedisplay();
} // Deactivate  


void  Select( const StrTy name, const IdxTy flags) 
{ 
// TODO FIXME need a mutex here.. 
auto ii=m_scenes.find(name);
if (ii!=m_scenes.end())
	if (m_select.find(name)==m_select.end())  m_select[name]=(*ii).second;
{auto ii=m_deselect.find(name);
if (ii!=m_deselect.end()) m_deselect.erase(ii);
}
glutPostRedisplay();
} // Activate  


void  Deselect( const StrTy name, const IdxTy flags) 
{ 
// TODO FIXME need a mutex here.. 
{auto ii=m_select.find(name);
if (ii!=m_select.end()) {  m_select.erase(ii); m_deselect[name]=(*ii).second; } 

}
glutPostRedisplay();
} // Deactivate  




bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);

static void _Display(void) { p()->Display(); } 
static void _Reshape(GLint width, GLint height) { p()->Reshape(width,height); } 
static void _MouseButton(int button, int state, int x, int y)
{ p()->MouseButton(button,state,x,y); } 
static void _MouseMotion(int x, int y) { p()->MouseMotion(x,y); } 
static void _AnimateScene(void) { p()->AnimateScene(); } 
static void _SelectFromMenu(int idCommand) { p()->SelectFromMenu(idCommand); } 
static void _SelectFromSubMenu(int idCommand) { p()->SelectFromSubMenu(idCommand); } 
static void _SpecialInput(int key, int x, int y) { p()->SpecialInput(key,x,y); } 
static void _Keyboard(unsigned char key, int x, int y) { p()->Keyboard(key,x,y); } 

void DrawHotZones(DrawInfo & di,const IdxTy flags)
{


} // DrawHotzones



void DrawTexts(DrawInfo& di)
{
typedef typename GlutUtil::color_t color_t; // ct(sf);
typedef typename GlutUtil::view_info_t ViewInfo; // ct(sf);
// pass things like w and h etc 
// but the screen view port is normalized to +/- 1 or should be 
// otherwise it fails without one active and in any case
// picks up that view. 
ViewInfo  v=ViewInfo(); // m_view;
GlutUtil::start_view(v);


MyGLStatus & gls= m_gl_status;
// TODO need an atomic flip here? 
m_gui=di.m_hz;
GuiLayout& gui=m_gui;
color_t white(1,1,1);

//gui.add_text(0,0,24,"MJMDatascope",white,0,BAD);
gui.add_text(0,0,24,m_product.c_str(),white,0,BAD);
//typename DrawInfo::hot_type::layout_type lt;
typename GuiLayout::layout_type lt;
lt.screen_size(gls.g_Width, gls.g_Height); 
gui.layout(lt);
// needs a 'view"?
// no need to end its the end now... 
GlutUtil::begin_screen_coords();

for(auto ii=gui.laid_begin(); ii!=gui.laid_end(); ++ii)
{
const auto & sf=((*ii));
//typename GlutUtil::color_t 
color_t ct(sf);
Gf r=sf.r;
Gf g=sf.g;
Gf b=sf.b;
glColor3f(r,g,b);
if (!false) 
	{ GlutUtil::render_string(sf.x, sf.y, 0 , sf.text.c_str(),ct); 
//MM_ERR(MMPR3(sf.x,sf.y,sf.text))

} // ,ctv[tc]);
if (false) { GlutUtil::draw_string(sf.x, sf.y, sf.z, sf.text); } // ,ctv[tc]);
if (false)
{
  glRasterPos2f(sf.x, sf.y);
//  glRasterPos3f(x, y,0.0);
//if (font==0)
  glutBitmapString(GLUT_BITMAP_TIMES_ROMAN_24, (const unsigned char *)(sf.text.c_str()));
//  else glutBitmapString(font, (const unsigned char *) string);
}


} // ii 

GlutUtil::end_screen_coords();

} // DrawTexts
void DealWithLocals(const DrawInfo & di)
{
// this is purely for testing or maybe some other thing late
// the working one is in rags.h 
if (m_svg.data())
{
// upside down wtf.. 
//glPixelZoom(1,-1);
// red and blue seem to be interchanged 
GLenum ty= GL_UNSIGNED_INT_8_8_8_8_REV;
//GLenum ty= GL_UNSIGNED_INT_8_8_8_8;
//GLenum ty= GL_UNSIGNED_BYTE;
//glDrawPixels(m_svg.w(),m_svg.h(),GL_RGBA, ty,m_svg.data());
glDrawPixels(m_svg.w(),m_svg.h(),GL_BGRA, ty,m_svg.data());
//glPixelZoom(1,1);

} // data 

} // DealWithLocals

void Display(void)
{
GlutUtil::start_view_zed();
DrawInfo di;
di.actives(m_actives.size());
//MM_ERR(" in display  "<<MMPR(this))
// TODO FIXME needs mutex 
IdxTy i=0;
MM_LOOP(ii,m_actives){
di.active(i);
 (*ii).second->draw(&di,0);
++i;
}
// finally add the  hot zone texts 
DrawTexts(di);
// svg Svg
DealWithLocals(di);
glutSwapBuffers();
//glFlush();
//MM_ERR(" flushing ")
m_displayed=true;
// Draw any hot items... 


if (m_psaver!=0)
{
IdxTy rc= m_psaver->display(m_save_params); 
// if 1, cancel active save... 
if (rc==1) m_save_params.cancel_capture();
}
//  RenderObjects();
} // Display

void Reshape(GLint width, GLint height)
{
// TODO FIXME this needs to use the same code as initially invoked...
MyGLStatus & gls= m_gl_status;
   gls.g_Width = width;
   gls.g_Height = height;
   glViewport(0, 0, gls.g_Width, gls.g_Height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   //gluPerspective(65.0, (float)g_Width / g_Height, g_nearPlane, g_farPlane);
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
  glutPostRedisplay();
}
IdxTy ResetViews(const IdxTy flags)
{
IdxTy  dispu=false;
const bool reset_all=!Bit(flags,0);
EnterSerial(0);
if (reset_all) 
	{ MM_LOOP(ii,m_scenes) {++dispu;  ((*ii).second->view()).reset(); }  } 
else 
	{ MM_LOOP(ii,m_actives) {++dispu;  ((*ii).second->view()).reset(); }  }
//if (m_active!=0 ) { auto & v=m_active->view(); v.reset(); }
 ExitSerial(0);
if (dispu)   glutPostRedisplay();
return dispu;
} // ResetViews 


// https://stackoverflow.com/questions/14378/using-the-mouse-scrollwheel-in-glut
void MouseButton(int button, int state, int x, int y)
{
//MM_ERR(MMPR4(button,state,x,y))
MyGLStatus & gls= m_gl_status;
// pressing wheel gives this 
//bool dispu=false;
// this is apparently the GLTU middle butteon and now
// is a menu trigger 
if (button==1) if (state==GLUT_DOWN)
{
ResetViews(0);
//EnterSerial(0);
//MM_LOOP(ii,m_actives) {dispu=true;  ((*ii).second->view()).reset(); } 
//if (m_active!=0 ) { auto & v=m_active->view(); v.reset(); }
// ExitSerial(0);

//if (dispu)   glutPostRedisplay();

return; 
}// reset

//if ( m_active==0) return; // TODO FIXME threading, globals etc 
//auto & v=m_active->view();
  // Respond to mouse button presses.
  // If button1 pressed, mark this state so we know in motion function.
  if (button == GLUT_LEFT_BUTTON)
    {
		gls.b1down((state==GLUT_DOWN),x,y);
		if (state==GLUT_UP)
		{
		    D zsfx=1; // 100.0/ gls.g_Width;	
		    D zsfy=1; // 100.0/ gls.g_Height;	
			const auto *  ev= m_gui.mouse_input(x*zsfx,y*zsfy,0);
			if (ev) 
			{
			(*((RagScene *)ev->po)).event(ev->code); 
/* typedef void * Ptr; Ptr po; IdxTy code; // event details. 
D x,y; IdxTy flags; */
			} // ev
		} // UP	
//      gls.g_bButton1Down = (state == GLUT_DOWN) ? TRUE : FALSE;
//      gls.g_yClick = y - 3 * gls.g_fViewDistance;
    }
  if (button == GLUT_RIGHT_BUTTON) gls.b2down((state==GLUT_DOWN),x,y);

EnterSerial(0);
//if (m_active==0 ) { ExitSerial(0); return; } 
bool dispu=false;
MM_LOOP(ii,m_actives)
{
dispu=true;
auto & v=(*ii).second->view();
// Wheel reports as button 3(scroll up) and button 4(scroll down)
   if ((button == 3) || (button == 4)) // It's a wheel event
   { 
       if (state == GLUT_DOWN)  // Disregard redundant GLUT_UP events
{
      if ( button ==3) v.distance*= 1.1; else  v.distance*= 1.0 / 1.1;
   // MM_ERR(" wheel fac 1.1 or .9  "<<MMPR2(v.distance, button))	

} // down 
       // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
     //  if (state == GLUT_UP) return; // Disregard redundant GLUT_UP events
     //  printf("Scroll %s At %d %d\n", (button == 3) ? "Up" : "Down", x, y);
  //  }//else{  // normal button event
     //  printf("Button %s At %d %d\n", (state == GLUT_DOWN) ? "Down" : "Up", x, y);
   }

}

ExitSerial(0);
if (dispu) {   glutPostRedisplay(); }
}
void MouseMotion(int x, int y)
{
MyGLStatus & gls= m_gl_status;
EnterSerial(0);
//if ( m_active==0)
//{
//ExitSerial(0);
// return; // TODO FIXME threading, globals etc 
//}
MM_LOOP(ii,m_actives)
{
auto & v=(*ii).second->view();
  // If button1 pressed, zoom in/out if mouse is moved up/down.
  //if (gls.g_bButton1Down)
  if (gls.b1down())
    {
//      gls.g_fViewDistance = (y - gls.g_yClick) / 3.0;
// ambiguous wrt rationals.. 
		v.distance=v.distance *::pow(1.1,double(y-gls.m_last_y));
//    MM_ERR(" mouse move zoom  "<<MMPR3(v.distance, y,gls.m_last_y))	


    }
  if (gls.b2down())
{
// TODO FIXME this needs to use GL coords doh 
// static IdxTy nproject(Gc & x, Gc & y, Gc & z, const Gc xin, const Gc yin)
D xw,yw,zw;
D xwl,ywl,zwl;
GlutUtil::unproject( xw,yw,zw,x,y);
//MM_ERR(MMPR4(x,y,xw,yw)<<MMPR(zw))
GlutUtil::unproject( xwl,ywl,zwl,gls.m_last_x,gls.m_last_y);
//MM_ERR(MMPR4(gls.m_last_x,gls.m_last_y,xwl,ywl)<<MMPR(zwl))
// this is really dum. It needs to invert the projection matrix
// and move tangential to the screen...
//auto & v=m_active->view();
//D fu=-1.0; // 1.0/2.0; // v.distance;
D fu=-1.0/ v.distance;
v.m_c[0]+=(xw-xwl)*fu;
v.m_c[1]-=(yw-ywl)*fu;
v.m_c[2]+=(zw-zwl)*fu;
//if (ctrl) {MM_ERR("y_c up "<<MMPR2(v.m_c[1],dy))  v.m_c[1]+=dy;} else
//if (alt) { MM_ERR(" rotate theta  ")   v.rot(-3.1415/20*f,0); }
//else { MM_ERR("roate phi ")   v.rot(0,-3.1415/20*f);}
} // gls.b2down

} // ii 

gls.last(x,y);
ExitSerial(0);
//  if (gls.g_bButton1Down)
if (gls.b1down()||gls.b2down())  glutPostRedisplay();
}
void AnimateScene(void)
{
/*
  float dt;
#ifdef _WIN32
  DWORD time_now;
  time_now = GetTickCount();
  dt = (float) (time_now - last_idle_time) / 1000.0;
#else
  // Figure out time elapsed since last call to idle function
  struct timeval time_now;
  gettimeofday(&time_now, NULL);
  dt = (float)(time_now.tv_sec  - last_idle_time.tv_sec) +
  1.0e-6*(time_now.tv_usec - last_idle_time.tv_usec);
#endif
  // Animate the teapot by updating its angles
  g_fTeapotAngle += dt * 30.0;
  g_fTeapotAngle2 += dt * 100.0;
  // Save time_now for next time
  last_idle_time = time_now;
 */
 // Force redraw
//MM_ERR(" animate redisplay ")  
usleep(200000);
if (!false)  glutPostRedisplay();
}
enum { POPUP_RESET,POPUP_TOGGLE_VIS };
//static 
// needs to be rebuilt when new scenes added. 
// TODO needs a mutex for m_scenees... 
// active, deactive to avoid auto-active, 
// select for manipulation 

int BuildPopupMenu (void)
{
glutDetachMenu(GLUT_MIDDLE_BUTTON);
//  int menu;
  m_vismenu = glutCreateMenu (Myt::_SelectFromSubMenu);
int i=0;
MM_LOOP(ii,m_scenes) { glutAddMenuEntry((*ii).first.c_str(),i); ++i; }
  m_selmenu = glutCreateMenu (Myt::_SelectFromSubMenu);
MM_LOOP(ii,m_scenes) { glutAddMenuEntry((*ii).first.c_str(),i); ++i; }

  m_menu = glutCreateMenu (Myt::_SelectFromMenu);
glutAddMenuEntry("Reset Views \tl",POPUP_RESET);
//if (m_scenes.size()) glutAddSubMenu("Visibility \tl",POPUP_TOGGLE_VIS);
if (m_scenes.size()) glutAddSubMenu("Visibility \tl",m_vismenu);
if (m_scenes.size()) glutAddSubMenu("Select \tl",m_selmenu);

 // glutAddMenuEntry ("Toggle polygon fill\tp", MENU_POLYMODE);
 // glutAddMenuEntry ("Toggle texturing\tt", MENU_TEXTURING);
 // glutAddMenuEntry ("Exit demo\tEsc", MENU_EXIT);
if (!false)   glutAttachMenu (GLUT_MIDDLE_BUTTON);
  return m_menu;
}


void SelectFromSubMenu(int idCommand)
{
const IdxTy nscenes=m_scenes.size(); // TODO hazard? 
MM_ERR(" sbimenu  "<<MMPR2(idCommand,nscenes))
bool sel=(idCommand>=nscenes);

auto ii=m_scenes.begin();
if (sel) idCommand-=nscenes;
while (idCommand) { ++ii; --idCommand;  } 
if (sel)
{
if (m_select.find((*ii).first)!=m_select.end())
Deselect((*ii).first,0);
else Select((*ii).first,0);

} // sel 
else
{
if (m_actives.find((*ii).first)!=m_actives.end())
Deactivate((*ii).first,0);
else Activate((*ii).first,0);
} // ! sel 


} // SelectFromSubMenu
void SelectFromMenu(int idCommand)
{
MM_ERR(" could have a menu for this ")
  switch (idCommand)
    {
		case POPUP_RESET : { ResetView(0); break; } 
 /*   case MENU_LIGHTING:
      g_bLightingEnabled = !g_bLightingEnabled;
      if (g_bLightingEnabled)
         glEnable(GL_LIGHTING);
      else
         glDisable(GL_LIGHTING);
      break;
    case MENU_POLYMODE:
      g_bFillPolygons = !g_bFillPolygons;
      glPolygonMode (GL_FRONT_AND_BACK, g_bFillPolygons ? GL_FILL : GL_LINE);
      break;      
    case MENU_TEXTURING:
      g_bTexture = !g_bTexture;
      if (g_bTexture)
         glEnable(GL_TEXTURE_2D);
      else
         glDisable(GL_TEXTURE_2D);
      break;    
    case MENU_EXIT:
      exit (0);
      break;
 */
    }
 // Almost any menu selection requires a redraw
  glutPostRedisplay();
}
// https://community.khronos.org/t/what-are-the-codes-for-arrow-keys-to-use-in-glut-keyboard-callback-function/26457

#define FOREACHA(x)  MM_LOOP(ii,m_actives){ auto & v=(*ii).second->view();  x } 
void SpecialInput(int key, int x, int y)
{
MyGLStatus & gls= m_gl_status;
//if ( m_active==0) return; // TODO FIXME threading, globals etc 
int sw=glutGetModifiers();
bool shift=sw&GLUT_ACTIVE_SHIFT;
bool ctrl=sw&GLUT_ACTIVE_CTRL;
bool alt=sw&GLUT_ACTIVE_ALT;
MM_ERR(MMPR4(key,shift,ctrl,alt))
Gf f=1.0;
Gf dx=.1;
Gf dy=.1;
if (shift) f=.1;
switch(key)
{
case GLUT_KEY_UP:
//do something here
{ MM_ERR(" down arrow ")
EnterSerial(0);
FOREACHA(
//auto & v=m_active->view();
if (ctrl) {MM_ERR("y_c up "<<MMPR2(v.m_c[1],dy))  v.m_c[1]+=dy;} else
if (alt) { MM_ERR(" rotate theta  ")   v.rot(-3.1415/20*f,0); } 
else { MM_ERR("roate phi ")   v.rot(0,-3.1415/20*f);} 
MM_ERR(MMPR(v.dump())) 
)
ExitSerial(0);
break;}
case GLUT_KEY_DOWN: { 
//do something here
MM_ERR(" down arrow ")
EnterSerial(0);
FOREACHA(
//auto & v=m_active->view();
//if (ctrl) {  v.m_c[1]-=dy;} else
//if (alt) v.rot(3.1415/20*f,0);
//else v.rot(0,3.1415/20*f);
if (ctrl) {MM_ERR("y_c down "<<MMPR2(v.m_c[1],dy))  v.m_c[1]-=dy;} else
if (alt) { MM_ERR(" rotate theta  ")   v.rot(3.1415/20*f,0); } 
else { MM_ERR("roate phi ")   v.rot(0,3.1415/20*f);} 
MM_ERR(MMPR(v.dump())) 
)
ExitSerial(0);
break;}
case GLUT_KEY_LEFT:
{MM_ERR(" left arrow ")
EnterSerial(0);
FOREACHA(
//auto & v=m_active->view();
if (ctrl) {  v.m_c[0]+=dx;} else
v.vidx(1);
MM_ERR(MMPR(v.dump()) )
)
ExitSerial(0);
//do something here
break;}
case GLUT_KEY_RIGHT:
{MM_ERR(" right arrow ")
EnterSerial(0);
FOREACHA(
//auto & v=m_active->view();
if (ctrl) {  v.m_c[0]-=dx;} else
v.vidx(-1);
MM_ERR(MMPR(v.dump()) )
)
ExitSerial(0);
//do something here
break;}
// don't do this, it passes shift as 112 which is "p" lol 
//default: Keyboard(key,x,y); 
}
  glutPostRedisplay();
} // SpecialInput

void Keyboard(unsigned char key, int x, int y)
{
MyGLStatus & gls= m_gl_status;
  switch (key)
  {
 /* case 27:             // ESCAPE key
	  exit (0); // not trapped? k
	  break;
  case 'l':
	  _SelectFromMenu(MENU_LIGHTING);
	  break;
  case 'p':
	  _SelectFromMenu(MENU_POLYMODE);
	  break;
  case 't':
	  _SelectFromMenu(MENU_TEXTURING);
	  break;
*/
case 'x': 
{
EnterSerial(0);
if (m_psaver) m_psaver->stop_capture(); 
ExitSerial(0);
}
case 'c': 
{
EnterSerial(0);
if (m_psaver) m_psaver->start_capture(); 
ExitSerial(0);
}

case 'v': 
{
EnterSerial(0);
if (m_psaver) m_psaver->pause_capture(); 
ExitSerial(0);
}



  case 'z':
EnterSerial(0);
MM_ERR(" z")
//if (m_active) m_views[""].m_scale[2]*=1.1;
//if (m_active) m_active->view().m_scale[2]*=1.1;
FOREACHA( v.m_scale[2]*=1.1; ) 

ExitSerial(0);
break;

  case 'Z':
EnterSerial(0);
MM_ERR(" z")
//m_views[""].m_scale[2]*=.9;
//if (m_active) m_active->view().m_scale[2]*=0.9;
FOREACHA( v.m_scale[2]*=0.9; ) 
ExitSerial(0);
break;

case ' ' : { ResetView(0); break; } 


  }
  glutPostRedisplay();
} // Keyboardj

IdxTy  ResetView(const IdxTy flags)
{
EnterSerial(0);
FOREACHA( v.reset(); 
//if (m_active==0) return 0; 
//auto & v=m_active->view();
//v.reset();
//if (ctrl) {  v.m_c[1]-=dy;} else
//if (alt) v.rot(3.1415/20*f,0);
//else v.rot(0,3.1415/20*f);
MM_ERR(MMPR(v.dump())) 
)
ExitSerial(0);
return 1; 
} // ResetView 



// don't pass this, need to put it somewhere safe, 
IdxTy Launch(const StrTy s, const IdxTy flags=0 )
{ 
ThParam* tp= new ThParam(p(), &Myt::Start, s,flags ) ;
tp->invoke();
// wait upto a second or two for display to draw 
IdxTy cnt=0; 
while (!m_displayed ) 
{ usleep(10000); ++cnt; if (cnt>200) { MM_ERR(" no display ") break; }  } 
return 0; 
}

IdxTy Start(const StrTy & s,const IdxTy flags=0 )
{
//int main(int argc, char** argv)
  // GLUT Window Initialization:
int argc=0;
char ** argv=0;
MyGLStatus & gls= m_gl_status;
  glutInit (&argc, argv);
  glutInitWindowSize (gls.g_Width, gls.g_Height);
glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
//  glutInitDisplayMode ( GLUT_RGB |  GLUT_DEPTH);
  //glutCreateWindow ("MJM datascope");
  glutCreateWindow (m_product.c_str());
   //void* pTextureImage=0;
   glEnable(GL_DEPTH_TEST);
   glDepthFunc(GL_LESS);
   glShadeModel(GL_SMOOTH);
   glEnable(GL_LIGHTING);
   glEnable(GL_LIGHT0);

// Register callbacks:
  glutDisplayFunc (Myt::_Display);
  glutReshapeFunc (Myt::_Reshape);
	glutSpecialFunc(Myt::_SpecialInput);
  glutKeyboardFunc (Myt::_Keyboard);
  glutMouseFunc (Myt::_MouseButton);
  glutMotionFunc (Myt::_MouseMotion);
// needed to redraw????
// 100 pct cpu  
// fing postredisoplay does not work from Add wtf? 
  glutIdleFunc (Myt::_AnimateScene);
 // glutIdleFunc (0);
  // Create our popup menu
// TODO create hot-area menu 
// can do both now... 
  BuildPopupMenu ();
//if (false)   glutAttachMenu (GLUT_RIGHT_BUTTON);
//if (!false)   glutAttachMenu (GLUT_MIDDLE_BUTTON);
  // Get the initial time, for use by animation
#ifdef _WIN32
  last_idle_time = GetTickCount();
#else
  gettimeofday (&m_gl_status.last_idle_time, NULL);
#endif
  // Turn the flow of control over to GLUT
m_alive=true;
// never sets matrix?? wtf TODO FIXME
MM_ERR(" glut loop START") std::cerr.flush();
try { 
  glutMainLoop ();
} catch (...) { m_alive=false;  MM_ERR(" glutMainLoop exit throws  " ); std::cerr.flush();  throw ;  } 
MM_ERR(" glut loop done ") std::cerr.flush();
m_alive=false; 
  return 0;
} // Start
 

void Free()
{
// TODO FIXME menu etc?


//delete[] g_lightPos;
} //Free

void Init()
{
m_product="MJMDatascope";
m_mutex_vector=MutexVector(4);
m_alive=false;
m_displayed=false;
m_thread=0;
m_gl_status.setup();
//m_active=0;
m_psaver=0;
m_menu=BAD;
//m_submenu=BAD;
m_vismenu=BAD;
m_selmenu=BAD;
const StrTy svgstr="";
// merely for testing stuff can remove m_svg etc or use for bakcground 
if (false) { m_svg.render(svgstr,0); }
} // Init


// MEMBERS - didn't know inits allowed here? wtf
StrTy m_product;
volatile bool m_alive=false; 
volatile bool m_displayed=false; 
ThreadId m_thread;
int m_menu,m_vismenu,m_selmenu;
//volatile 
//Scenes m_scenes;
//StrTy m_scene;
//Views m_views;
//RagScene * m_active;
//ModelPool m_pool;
// These are just pointers... 
RagSceneMap m_scenes,m_actives,m_deactive,m_select,m_deselect;
MyGLStatus m_gl_status;
Saver * m_psaver;
Sp m_save_params;
GuiLayout m_gui;
// tacked in
SvgRender m_svg;
}; // mjm_glut_scope_ii

//////////////////////////////////////////////

template <class Tr>
class mjm_glut_scope_ii_map : public std::map<typename Tr::StrTy, mjm_glut_scope_ii< Tr > >  
{
 typedef mjm_glut_scope_ii_map Myt;
typedef typename std::map<typename Tr::StrTy, mjm_glut_scope_ii< Tr> >   Super;
 typedef typename Tr::IdxTy IdxTy;
 typedef typename Tr::D D;
 typedef typename Tr::StrTy StrTy;
 typedef typename Tr::Ss Ss;
 typedef typename Tr::IsTy IsTy;
 typedef typename Tr::OsTy OsTy;
 typedef typename Tr::Ofs Ofs;
//typedef typename Tr::FlagTy; 
// typedef typename Tr::MyBlock  MyBlock;
public:
mjm_glut_scope_ii_map() {}
StrTy dump(const IdxTy flags=0) { return Dump(flags); }
private:
bool Bit(const IdxTy f, const IdxTy b) const  { return  ((f>>b)&1)!=0; }
// should loop over map now 
//StrTy Dump(const IdxTy flags=0) {Ss ss;  return ss.str(); }
typedef typename mjm_thread_util<Tr>::mutex_vector MutexVector;

enum { MAP_MU=0 , MU_SZ};
mutable MutexVector m_mutex_vector;
void EnterSerial(const IdxTy i)const  {  m_mutex_vector.enter_serial(i ); }
void ExitSerial(const IdxTy i)const  {  m_mutex_vector.exit_serial(i ); }
//m_mutex_vector = MutexVector(MU_SZ);


//StrTy dump(const IdxTy flags=0) { return Dump(flags); }

private:

void Init()
{
//m_mutex_vector = MutexVector(10);
}

StrTy Dump(const IdxTy flags=0)
{
Ss ss;
MM_LOOP(ii,(*this))
{
ss<<(*ii).first<<CRLF;
ss<<(*ii).second.dump()<<CRLF;


}
return ss.str();
// return Dump(flags); 

}




private:

}; // mjm_glut_scope_ii_map




////////////////////////////////////////////
#ifdef  TEST_MJM_GLUT_SCOPE_II
class Tr {
public:
// typedef mjm_string_picker Myt;
 typedef unsigned int IdxTy;
 typedef double  D;
 typedef std::string StrTy;
 typedef std::stringstream Ss;
 typedef std::istream  IsTy;
 typedef std::ostream  OsTy;
 typedef std::ofstream  Ofs;
// typedef typename Tr::MyBlock  MyBlock;
}; // 


#include "mjm_instruments.h"
#include "mjm_cli_ui.h"
typedef Tr::StrTy StrTy;
typedef Tr::IdxTy IdxTy;

template <class Tt> class tester_ {
typedef tester_<Tt> Myt;
typedef mjm_cli_ui<Myt> Cli;
//typedef tester Myt;
//typedef mjm_cli_ui<Myt> Cli;
typedef std::map<StrTy, StrTy> LocalVar;

typedef CommandInterpretterParam Cip ;
typedef void (Myt:: * CmdFunc)(Cip &, LocalVar &  ) ;
typedef std::map<StrTy, CmdFunc> CmdMap;
typedef std::vector<StrTy> Choices;
//typedef void (Myt:: * CompleteFunc) ( Cli::list_type & choices,  const char * cmd, const char * frag);
typedef void (Myt:: * CompleteFunc) ( Choices & choices,  const char * cmd, const char * frag);
typedef std::map<StrTy, CompleteFunc> CompMap;

public:
 //void cli_cmd( Cli::list_type & choices,  const char * frag)
 void cli_cmd( Choices & choices,  const char * frag)
{
const IdxTy nfrag=strlen(frag);
MM_LOOP(ii,m_cmd_map)
{
const StrTy & v=(*ii).first;
if (strncmp(v.c_str(),frag,nfrag)==0)  choices.push_back(v);
}
}

 //void cli_param( Cli::list_type & choices,  const char * _cmd, const char * frag)
 void cli_param( Choices & choices,  const char * _cmd, const char * frag)
{
MM_ERR("cli_param"<<MMPR2(_cmd,frag))
//const StrTy cmd=CliTy::word(StrTy(_cmd),0);
//auto ii=m_comp_map.find(cmd);
//if ( ii!=m_comp_map.end()) ((this)->*(*ii).second)(choices,cmd.c_str(),frag);
}

CmdMap m_cmd_map;


 }; // tester_
typedef tester_< mjm_glut_scope_ii <Tr>  > tester;

typedef mjm_cli_ui<tester> Cli;


typedef Tr::Ss Ss;
void about()
{
Ss ss;
ss<<" MJM_GLUT_SCOPE_II "<<__DATE__<<" "<<__TIME__<<CRLF;
MM_ERR(ss.str())

}

int main(int argc,char **args)
{
about();
typedef mjm_glut_scope_ii<Tr>  Myt;
//Myt x(argc,args);
Myt&  x= *Myt::p();

//if (!x.done()) x.command_mode();
Cli cli;
tester tester;
CommandInterpretter li(&std::cin);
li.push(args,argc);
cli.set_target(tester);
cli.set_command_handler(&tester::cli_cmd);
cli.set_param_handler(&tester::cli_param);
cli.activate();
li.set_split(1,' ');
while (li.nextok())
{
const IdxTy sz=li.size();
if (sz<1) continue;
const StrTy cmd=li.word(0);
if (cmd=="") continue;
if (cmd=="about"){ about();  continue; } 
CommandInterpretterParam  cip(li);

if (cmd=="quit") break;
if (cmd=="dump") { MM_ERR(x.dump()) }
if (cmd=="launch") { MM_ERR(x.launch(cip.p1)) }
//if (cmd=="stream") { MM_ERR(x.stream(cip.p1)) }
//if (cmd=="start") { MM_ERR(x.start(cip.p1)) }
if (cmd=="wait") { StrTy x; std::cin>> x ; }
//else if (cmd=="load") { x.load(li.words(),1); }
//else if (cmd=="clear") { x.clear(); }

} // nextok

//if (!x.done()) x.command_mode();
return 0;
}

#endif // main

#endif // MJM_GLUT_SCOPE_II_H__ 
#undef TRUE 
#undef FALSE 


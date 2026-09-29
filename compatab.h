/*compatab.h : borland's type of graphic for WATCOM, GNU compilers */
/*Astrakharchik Gregory '98 gregory@imp.inm.ras.ru */

#ifndef __COMPATAB_H_
#define __COMPATAB_H_

/* undefine to disable graphics in linux */
//#define LINUX_GRAPHICS 

//operating systems: _WIN32, _WIN64, __unix__, __APPLE__
//compilers: __INTEL_COMPILER, __GNUC__, _MSC_VER

#define __DUMB__
//#define __MVC__ // Windows graphics

#ifdef __WATCOMC__
# undef __DUMB__
# define _found
#endif

#ifdef __linux__
# define INLINE static inline
# define _found
#
# ifdef LINUX_GRAPHICS
#   undef __DUMB__
# endif
# define _found
#endif

#ifdef __TURBOC__
# define _found
#endif

/*#ifdef __GNUC__
# define _found
# define INLINE inline
#endif*/

#ifdef __SUN__
# define _found
#endif

/*#ifndef _found
//#  error "don't know such a compiler"
# define _found
// Borland C Builder assumed
//#define BCBuilder 
#  define INLINE __inline
#endif*/

#ifdef __linux__
#  define SLASH "/"
#else
#  define SLASH "\\"
#endif

/******************************** BORLAND *********************************/
#ifdef __TURBOC__
#include <graphics.h>
#endif

/********************************* WATCOM **********************************/
#ifdef __WATCOMC__

#include <graph.h>
#include <math.h>
#include <time.h>
#include <conio.h>

#define DETECT 0
#define grOk 0
#define M_PI 3.14159265358979

#define graphresult() 0

#define grapherrormsg(a) "OK"

extern videoconfig _videocfg;
extern short _fill;

inline void setfillstyle(int a, int b) {
_setcolor(b);
if(a)
  _fill = _GFILLINTERIOR;
else
  _fill = _GBORDER;
}

inline void initgraph(int *a, int *b, char *c) {
        /*_setvideomode(_MAXRESMODE);*/
        /*_setvideomode(_SVRES16COLOR);*/
        _setvideomode(_XRES256COLOR);
        _setcharsize(11,9);
        _setcharspacing(0);
        _getvideoconfig(&_videocfg);
}

inline void line(int x1, int y1, int x2, int y2) {
_moveto(x1, y1);
_lineto(x2,y2);
}

inline void putpixel(int x, int y, int colour) {
int colour2;
colour2 = _getcolor();
_setcolor(colour);
_setpixel(x,y);
_setcolor(colour2);
}

#define getmaxx() _videocfg.numxpixels
#define getmaxy() _videocfg.numypixels
#define getmaxcolor() _videocfg.numcolors
#define closegraph() _setvideomode(_DEFAULTMODE)
#define cleardevice() _clearscreen(0)
//#define ellipse(x, y, alpha, beta, xradius, yradius) _pie(_GBORDER, x-xradius, y-yradius, x+xradius, y+yradius, x+cos(alpha*M_PI/180), y+sin(alpha*M_PI/180), x+cos(beta*M_PI/180), y+sin(beta*M_PI/180))
//#define ellipse(x, y, alpha, beta, xradius, yradius) _pie(_GBORDER, x-xradius, y-yradius, x+xradius, y+yradius, x+1, y+1, x+2, y+2)
#define ellipse(x, y, a, b, xradius, yradius) _ellipse(_GBORDER, x-xradius, y+yradius, x+xradius, y-yradius)
#define fillellipse(x, y, xradius, yradius) _ellipse(_GFILLINTERIOR, x-xradius, y-yradius, x+xradius, y+yradius)
#define circle(x,y,radius) _ellipse(_GBORDER, x-radius, y-radius, x+radius, y+radius)
#define bar(a,b,c,d) _rectangle (_GFILLINTERIOR,a,b,c,d)
#define outtextxy _grtext
#define moveto _moveto
#define lineto _lineto
#define getpixel _getpixel
#define floodfill _floodfill
#define rectangle(a,b,c,d) _rectangle(_fill,a,b,c,d)
#define setcolor _setcolor

#define randomize() srand((unsigned)time(NULL))
#define random(a) (a <= 0) ? 0 : (rand() % a)

#endif

/********************************* LINUX ***********************************/
#ifdef __linux__
#define sound(a)
#define nosound()
#define randomize() srand((unsigned)time(NULL))
//#define random(a) (a <= 0) ? 0 : (rand() % a)

//void delay(unsigned long milisec);

#include<time.h>  
#include<signal.h>  

/*void delay(unsigned long milisec) {  
  struct timespec req={0},rem={0};  
  time_t sec=(int)(milisec/1000);  
  milisec=milisec-(sec*1000);  
  req.tv_sec=sec;  
  req.tv_nsec=milisec*1000000L;  
  nanosleep(&req,&rem);  
} */ 

//#define __STRICT_ANSI__

/********************************* LINUX SVGA lib ****************************/
#ifdef LINUX_GRAPHICS

#include <vga.h>

#define initgraph(a,b,c) vga_setmode(4)
#define getmaxx vga_getxdim
#define getmaxy vga_getydim
#define getmaxcolor vga_getcolors
#define closegraph() vga_setmode(vga_getdefaultmode())
#define cleardevice vga_clear
#define setcolor  vga_setcolor
#define setfillstyle(a,b)
#define bar(a,b,c,d)
#define outtextxy
#define circle(x,y,radius) {line(x-radius, y, x, y-radius);line(x,y-radius, x+radius,y);line(x+radius,y,x,y+radius);line(x,y+radius,x-radius,y);}
#define moveto
#define lineto
#define line vga_drawline
#define getpixel vga_getpixel
#define getpalette vga_getpalette
#define setpalette vga_getpalette
#define getgraphmode vga_getcurrentmode

#define getch vga_getch
#define kbhit vga_getkey

static __inline__ void putpixel(int x, int y, int colour) {
  vga_setcolor(colour);
  vga_drawpixel(x, y);
  return;
}

#define graphresult() 0
#define grapherrormsg

#define _setcharsize(a,b)
#define _setlinestyle(a)
#define _grtext(a,b,c)
#endif
#endif

/******************************* SUN GNU c ************************************/
#ifdef __SUN__

#include <curses.h>

#define randomize() srand((unsigned)time(NULL))
#define random(a) (a <= 0) ? 0 : (rand() % a)
#define itoa(value, str, radix) sprintf((str),"%d",(value))

#endif

/*#ifdef __GNUC__

#include <libbcc.h>

#define randomize() srand((unsigned)time(NULL))
#define random(a) (a <= 0) ? 0 : (rand() % a)
#define itoa(value, str, radix) sprintf((str),"%d",(value))
#define getch() getkey()

#endif*/

/************************ Borland C Builder ********************************/
#ifdef BCBuilder

#include <vcl\Classes.hpp>
#include <vcl\Controls.hpp>
#include <vcl\StdCtrls.hpp>
#include <vcl\Forms.hpp>

class TForm1 : public TForm
{
__published:	// IDE-managed Components
        void __fastcall FormPaint(TObject *Sender);
private:	// User declarations
public:		// User declarations
        __fastcall TForm1(TComponent* Owner);
};

extern TForm1 *Form1;

#define DETECT 0
#define grOk 0
#define closegraph()
#define initgraph(a,b,c)
#define graphresult() 0
#define grapherrormsg
#define delay(a)
#define outtext(a)
#define outtextxy Form1->Canvas->TextOut
#define getmaxx() Form1->Width
#define getmaxy() Form1->Height
#define lineto Form1->Canvas->LineTo
#define moveto Form1->Canvas->MoveTo

inline int __color(int i) {

int clr[16]={
0, //clBlack (TColor)(0)
//711680, //clBlue (TColor)
//5, clRed (TColor)
128,//Maroon (TColor)
32768, //clGreen (TColor)
32896,//clOlive (TColor)
8388608,//#define clNavy (TColor)
8388736,//#define clPurple (TColor)
8421376,//#define clTeal (TColor)
8421504,//#define clGray (TColor)
12632256,//#define clSilver (TColor)
65280,//#define clLime (TColor)
65535,//#define clYellow (TColor)
16711935,//#define clFuchsia (TColor)
16776960,//#define clAqua (TColor)
12632256,//#define clLtGray (TColor)
8421504,//#define clDkGray (TColor)
16777215//#define clWhite (TColor)
};

//(536870911),//#define clNone (TColor)

/*#define clDefault (TColor)(536870912)
#define cmBlackness (Byte)(66)
#define cmDstInvert (int)(5570569)
#define cmMergeCopy (int)(12583114)
#define cmMergePaint (int)(12255782)
#define cmNotSrcCopy (int)(3342344)
#define cmNotSrcErase (int)(1114278)
#define cmPatCopy (int)(15728673)
#define cmPatInvert (int)(5898313)
#define cmPatPaint (int)(16452105)
#define cmSrcAnd (int)(8913094)
#define cmSrcCopy (int)(13369376)
#define cmSrcErase (int)(4457256)
#define cmSrcInvert (int)(6684742)
#define cmSrcPaint (int)(15597702)
#define cmWhiteness (int)(16711778)
*/
return clr[i];
}

inline void setcolor(int c) {
Form1->Canvas->Pen->Color = (TColor)(__color(c));
}

inline circle(int x,int y,int r) {
Form1->Canvas->Brush->Style = bsClear;
Form1->Canvas->Brush->Color = clNone;
Form1->Canvas->Ellipse(x-r,y-r,x+r,y+r);
/*bsSolid, bsClear, bsHorizontal, bsVertical, bsFDiagonal, bsBDiagonal,
bsCross, bsDiagCross*/
}

#include <math.h>

inline ellipse(int x,int y,int alpha,int beta,int rx,int ry) {
Form1->Canvas->Brush->Style = bsClear;
Form1->Canvas->Brush->Color = clNone;
Form1->Canvas->Arc(x-rx,y-ry,x+rx,y+ry,x+rx*cos(alpha/180*3.14159265468979),y+ry*sin(alpha/180*3.14159265468979),x+rx*cos(beta/180*3.14159265468979),y+ry*sin(beta/180*3.14159265468979));
}

inline void line(int x1, int y1, int x2, int y2) {
Form1->Canvas->MoveTo(x1, y1);
Form1->Canvas->LineTo(x2, y2);
}

inline void rectangle(int x1, int y1, int x2, int y2) {
Form1->Canvas->Rectangle(x1,y1,x2,y2);
}

inline int getpixel(int x, int y) {
 int c;
 int i;
 c=(int)Form1->Canvas->Pixels[x][y];
 for (i=0;i<16;i++) if (__color(i)==(int)c) return i;
 return 1;
}

inline int floodfill(int x, int y,int c) {
Form1->Canvas->FloodFill(x,y,__color(c),fsBorder);
}

#define cleardevice Form1->Canvas->Refresh
#define bar Form1->Canvas->Rectangle

inline void setfillstyle(int a, int b){
if (a)
   Form1->Canvas->Brush->Style = bsSolid;
else
   Form1->Canvas->Brush->Style = bsClear;

Form1->Canvas->Brush->Color = (TColor) __color(b);
}

inline void drawpoly(int a, int *b) {
   int i;

   moveto(b[0], b[1]);

   for(i=1; i<a; i++)
      lineto(b[2*i], b[2*i+1]);

   lineto(b[0], b[1]);
}

#define fillpoly drawpoly
#define putpixel(x,y,c) Form1->Canvas->Pixels[x][y]=__color(c)

#define main Main
#endif

/*Definitions for the dumb display */
#ifdef __DUMB__

#define initgraph(a,b,c)
#define graphresult() 0
#define grapherrormsg(a) "This is simulation of the dumb display"
#define getmaxx() 1
#define getmaxy() 1
#define getmaxcolor() 1
#define closegraph() 
#define cleardevice() 
#define ellipse(x, y, alpha, beta, xradius, yradius)
#define line(a,b,c,d) 
#define fillellipse(x, y, xradius, yradius)
#define circle(x,y,radius)
#define bar(a,b,c,d)
#define outtextxy
#define moveto(a,b)
#define lineto(a,b)
#define getpixel(a,x) 1
#define floodfill(a)
#define rectangle(a,b,c,d)
#define setcolor(c)
#define putpixel(a,b,c)
#define kbhit() 0
#define getch() 0

#define _setlinestyle(a)
#define _grtext(a,b,c)
#define _setcharsize(a,b)

#endif

/* General constants */

#define ON 1
#define OFF 0

#define _E         2.71828182845904523536
#define _LOG2E     1.44269504088896340736
#define _LOG10E    0.434294481903251827651
#define _LN2       0.693147180559945309417
#define _LN10      2.30258509299404568402

#define _PI        3.14159265358979323846264338327950288
#define _PI_2      1.57079632679489661923
#define _PI_4      0.785398163397448309616
#define _2PI       6.2831853071795864769252867665590058
#define _1_PI      0.318309886183790671538
#define _2_PI      0.636619772367581343076
#define _SQRTPI    1.772453850905515840000
#define _1_SQRTPI  0.564189583547756286948
#define _2_SQRTPI  1.12837916709551257390 // 2 / sqrt(pi)

#define _SQRT2     1.41421356237309504880
#define _SQRT_2    0.707106781186547524401
#define gammaE     0.577215664901532860606512090082402431042 // Euler-Mascheroni constant

#define pi         _PI
#define PI         _PI
#define PI2        9.8696044010893586188
#define PI3        31.006276680299820175
#define PI4        97.409091034002437236
#define PI5        306.01968478528145326
#define PI6        961.38919357530443703

#define SOLID  0xFFFF
#define DASHED 0xF0F0

#define BLACK          0
#define BLUE           1
#define GREEN          2
#define CYAN           3
#define RED            4
#define MAGENTA        5
#define BROWN          6
#define WHITE          7
#define GRAY           8
#define LIGHTBLUE      9
#define LIGHTGREEN     10
#define LIGHTCYAN      11
#define LIGHTRED       12
#define LIGHTMAGENTA   13
#define YELLOW         14
#define BRIGHTWHITE    15


//Definitions for the Microsoft Visual Studio
#ifdef __MVC__

#define _CRT_SECURE_NO_WARNINGS // no fopen() warning, define in compiler options

#include "Windows.h"
#include <stdio.h>

//#define MAXX 1920
//#define MAXY 1150
extern int MAXX;
extern int MAXY;

extern HWND MyHWND;
extern HDC MyHDC;
extern const SMALL_RECT MyRect;

LOGPALETTE MyPLTE;
HANDLE MyConsoleHndl;
//LPPAINTSTRUCT MylpPaint;
PAINTSTRUCT MylpPaint;

__inline void initgraph(int *a, int *b, char *c) {

  MyHWND = FindWindow("tty", "main");
  MyHDC = GetWindowDC(MyHWND);
  //MyHDC= BeginPaint(MyHWND, &MylpPaint);
  //MylpPaint.rcPaint.left);

  SetBkColor(MyHDC, 0);
  SetTextColor(MyHDC, RGB(255,255,255));

  MAXX = GetDeviceCaps(MyHDC, HORZRES);
  MAXY = GetDeviceCaps(MyHDC, VERTRES);

  printf("Initializing screen:\n  width  %i pixels, %i mm\n", GetDeviceCaps(MyHDC, HORZRES), GetDeviceCaps(MyHDC, HORZSIZE));
  printf("  height %i pixels, %i mm\ndone\n", GetDeviceCaps(MyHDC, VERTRES), GetDeviceCaps(MyHDC, VERTSIZE));
}

HPEN MyPen;

__inline COLORREF EGA2RGB(int c) {
  switch(c%16) {
    case 0 : return RGB(000,000,000);
    case 1 : return RGB(000,000,168);
    case 2 : return RGB(000,168,000);
    case 3 : return RGB(000,168,168);
    case 4 : return RGB(168,000,000);
    case 5 : return RGB(168,000,168);
    case 6 : return RGB(168,168,000);
    case 7 : return RGB(208,208,208);
    case 8 : return RGB(168,168,168);
    case 9 : return RGB(000,000,255);
    case 10: return RGB(000,255,000);
    case 11: return RGB(000,255,255);
    case 12: return RGB(255,000,000);
    case 13: return RGB(255,000,255);
    case 14: return RGB(255,255,000);
    case 15: return RGB(255,255,255);
    default: return RGB(255,255,255);
  }		
}
__inline void setcolor(int c) {
  MyPen = CreatePen(PS_SOLID, 0, EGA2RGB(c));
  SelectObject(MyHDC, MyPen);
}

#define putpixel(a,b,c) SetPixel(MyHDC, a, b, EGA2RGB(c))
#define line(a,b,c,d) {MoveToEx(MyHDC, a, b, NULL); LineTo(MyHDC, c, d);}
#define ellipse(x, y, alpha, beta, xradius, yradius) Ellipse(MyHDC, x-xradius, y+yradius, x+xradius, y-yradius)
#define circle(x,y,radius) Ellipse(MyHDC, x-radius, y+radius, x+radius, y-radius);
#define bar(a,b,c,d) PatBlt (MyHDC, a, b, c, d, BLACKNESS)
#define outtextxy(x,y,string) TextOut(MyHDC, x, y, string, strlen(string))
#define cleardevice() SetBkMode(MyHDC, OPAQUE)

#define getmaxcolor() 16
#define getmaxx() MAXX
#define getmaxy() MAXY

#define graphresult() 0
#define grapherrormsg(a) "This is simulation of the dumb display"
#define closegraph() 
#define fillellipse(x, y, xradius, yradius)
#define floodfill(a)
#define getpixel(a,x) 1
#define kbhit() 0
#define getch() 0
#define moveto(a,b)
#define lineto(a,b)
#define rectangle(a,b,c,d)

#define _setlinestyle(a)
#define _grtext(a,b,c)
#define _setcharsize(a,b)

//#define _CRT_SECURE_NO_DEPRECATE
//Or add _CRT_SECURE_NO_DEPRECATE to Project Properties -> Configuration Properties -> C/C++ -> Preprocessor -> Preprocessor Definitions.
#endif
/////////// visual c

// math functions
#define tg(x) tan(x)
#define arctg(x) atan(x)
#define th(x) tanh(x)
#define sh(x) sinh(x)

#ifdef __INTEL_COMPILER
#ifdef __linux__
#  define MATHINCLUDE "mathimf.h"
#else
#  define MATHINCLUDE "math.h"
#endif
#else // not INTEL
#ifdef __GNUC__
#  define MATHINCLUDE "math.h"
#else
#  define MATHINCLUDE "math.h"
#  define erfc(x) erfc_my(x) 
#endif
#endif

///////// OPENMP

#ifndef _OPENMP
#  define omp_get_thread_num() 0
#endif

///// Windows
// define isnan() in old versions of visual studio
#if _MSC_VER == 1310 //(Visual Studio 2003)
# define isnan_not_defined
#endif

#if _MSC_VER == 1800 // (Visual Studio 2013)
# define isnan_not_defined
#endif

#if _MSC_VER == 1700 // (Visual Studio 2012)
# define isnan_not_defined
#endif

#if _MSC_VER == 1600 // (Visual Studio 2010)
# define isnan_not_defined
#endif

#if _MSC_VER == 1500 // (Visual Studio 2008)
# define isnan_not_defined
#endif

#if _MSC_VER == 1400 // (Visual Studio 2005)
# define isnan_not_defined
#endif

#ifdef _WIN32
#include <float.h>
#ifdef isnan_not_defined
#  define isnan _isnan
#endif
#ifdef _MSC_VER // visual studio
//#define _CRT_SECURE_NO_WARNINGS
#endif
#endif


#endif

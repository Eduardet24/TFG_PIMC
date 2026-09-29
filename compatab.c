/*compatab.cpp*/

/*include "compatab.h" in <main>.cpp

add
1) <main>.cpp
2) compatab.cpp
to the project*/

#include "compatab.h"

#ifdef BCBuilder

#include <vcl\vcl.h>
#pragma hdrstop
#pragma resource "*.dfm"

TForm1 *Form1;
__fastcall TForm1::TForm1(TComponent* Owner)
        : TForm(Owner){}


void __fastcall TForm1::FormPaint(TObject *Sender)
{
Form1->Width = 639;
Form1->Height = 479;
/*Form1->Width = 800;
Form1->Height = 600;*/

int Main(void);
Main();
}

#endif

#ifdef __WATCOMC__

videoconfig _videocfg;
char* _error_message = "graphics uses emulator for Watcom";
short _fill = _GBORDER;

#endif

#ifdef __MVC__

HWND MyHWND;
HDC MyHDC;

/*const int _clr[16] = {
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
};*/

//const SMALL_RECT MyRect = {0,0,MAXX,MAXY};
int MAXX=800;
int MAXY=600;

const SMALL_RECT MyRect = {0,0,800,600};

#endif

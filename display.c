/*display.c*/

#include <stdio.h>
#include "main.h"
#include "display.h"
#include "utils.h"
#include "optimiz.h"
#include "compatab.h"
#include MATHINCLUDE

/******************************* Close Graph ********************************/
void CloseGraph(void) {
  closegraph();
}

/********************************** Picture *********************************/
void Picture(void) {
  int i;
  static int first_time = ON;
  static int ncolors;
  static int maxx;
  static int maxy;
  static int maxx2;
  static int maxy2;
  static int t=0;
  extern DOUBLE E;

  if(kbhit()) {
    if(getch() == 27) {
      CloseGraph();
      Exit(0, "Escape pressed\n");
    }
  }

  if(first_time == ON) {
    i = 0;
    initgraph(&i, &i, "");
    if(graphresult()) {
      Error("Graphics error: %s\n", grapherrormsg(graphresult()));
    }
    ncolors = getmaxcolor();
    maxx2 = getmaxx()/2;
    maxy2 = getmaxy()/2;
    maxx = getmaxx();
    maxy = getmaxy();
  }

  if(optimization == 2) { // energy optimization
    if(first_time) {
      setcolor(250);
      bar(0, 0, maxx, maxy);
      setcolor(15);
    }
    //PicturePlot(0, "Energy (green), mean energy (yellow)", 1, t, E, 0, maxx2, 0, maxy, first_time, LIGHTGREEN);
    PicturePlot(0, "Energy (green), mean energy (yellow)", 1, t, E, 0, maxx2, 0, maxy2, first_time, LIGHTGREEN);
    PicturePlot(4, "sigma (green), mean energy (yellow)", 4, t, opt_sigma, 0, maxx2, maxy2, maxy, first_time, LIGHTGREEN);
#ifdef CRYSTAL
#ifdef CRYSTAL_WIDTH_ARRAY
    PicturePlot(1, "positions", 5, t, Crystal.Rz[0], maxx2, maxx, 0, maxy2, first_time, LIGHTBLUE);
#endif
#endif
    /*PicturePlot(1, "positions", 5, t, Crystal.Rz[1], maxx2, maxx, 0, maxy2, first_time, LIGHTGREEN);
    PicturePlot(1, "positions", 5, t, Crystal.Rz[2], maxx2, maxx, 0, maxy2, first_time, LIGHTCYAN);
    PicturePlot(1, "positions", 5, t, Crystal.Rz[3], maxx2, maxx, 0, maxy2, first_time, LIGHTRED);
    PicturePlot(1, "positions", 5, t, Crystal.Rz[4], maxx2, maxx, 0, maxy2, first_time, LIGHTMAGENTA);
    //PictureDistribution(1, 4, maxx2, maxx, 0, maxy, first_time, GREEN);
    PicturePlot(2, "positions", 6, t, fabs(Crystal.z[0]), maxx2, maxx, maxy2, maxy, first_time, LIGHTBLUE);
    PicturePlot(2, "positions", 6, t, fabs(Crystal.z[1]), maxx2, maxx, maxy2, maxy, first_time, LIGHTGREEN);
    PicturePlot(2, "positions", 6, t, fabs(Crystal.z[2]), maxx2, maxx, maxy2, maxy, first_time, LIGHTCYAN);
    PicturePlot(2, "positions", 6, t, fabs(Crystal.z[3]), maxx2, maxx, maxy2, maxy, first_time, LIGHTRED);
    PicturePlot(2, "positions", 6, t, fabs(Crystal.z[4]), maxx2, maxx, maxy2, maxy, first_time, LIGHTMAGENTA);*/
  }
  else if(boundary == THREE_BOUNDARY_CONDITIONS) {
    PicturePlot(0, "Energy (green), mean energy (yellow)", 1, t, E, 0, maxx, 40, maxy2, first_time, LIGHTGREEN);
    PicturePlot(0, "", 1, t, Evar/(DOUBLE)Nmeasured, 0, maxx, 40, maxy2, OFF, YELLOW);
    if(MC == VARIATIONAL) {
      PicturePlot(0, "", 1, t, EFF,   0, maxx, 40, maxy2, OFF, LIGHTCYAN);
      if(t*Nmeasure%Niter == 0) {// i.e. draw each block
        //PictureCoordinates('x', 'y', 0, (int)(maxx*0.33-2), maxy2+2, maxy-1);
        //PictureCoordinates('x', 'z', (int)(maxx*0.33), (int)(0.66*maxx)-2, maxy2+2, maxy-1);
        //PictureCoordinates('y', 'z', (int)(maxx*0.66), maxx-1, maxy2+2, maxy-1);
      }
      if(t*Nmeasure%Niter == 0) {
        PictureDistribution(1, 2, 0, maxx2, maxy2, maxy, first_time, GREEN);
        PictureDistribution(2, 1, maxx2, maxx, maxy2, maxy, first_time, GREEN);
      }
    }
    else if(MC == DIFFUSION) { // i.e. DMC
      PicturePlot(1, "Population of walkers", 2, t, Nwalkersw, 0, maxx, maxy2, maxy, OFF, LIGHTBLUE);
      PicturePlot(1, "Population of walkers", 2, t, (DOUBLE) Nwalkers, 0, maxx, maxy2, maxy, first_time, YELLOW);
    }
  }
  else {
    if(MC == VARIATIONAL) {
      if(t*Nmeasure%Niter == 0) {
        //PictureCoordinates('z', 'y' ,0, maxx*0.32, maxy-0.33*maxx, maxy-1);
        //PictureCoordinates('z', 'x', maxx*0.33, maxx*0.64, maxy-0.33*maxx, maxy-1);
        //PictureCoordinates('y', 'x', (int)(maxx*0.67), (int)(maxx-1), (int)(maxy-0.33*maxx), (int)(maxy-1));
      }
      /* plot energy */
      PicturePlot(0, "Energy", 1, t, E, 0, maxx, 40, maxy2, first_time, LIGHTGREEN);
      PicturePlot(0, "", 1, t, EFF, 0, maxx, 40, maxy2, OFF, LIGHTRED);
      //PicturePlot(0, "", 1, t, Evar/(DOUBLE)Nmeasured, 0, maxx, 40, maxy2, OFF, YELLOW);
      /*plot <r2> and <z2> */
      if(measure_R2) {
        //PicturePlot(1, "Sqrt(<r2>)", 4, t, RD.r2recent, 0, (int)(maxx*0.33), maxy2, maxy, first_time, LIGHTGREEN);
        //PicturePlot(2, "Sqrt(<z2>)", 3, t, RDz.r2recent, (int)(maxx*0.33), (int)(maxx*0.66), maxy2, maxy, first_time, LIGHTGREEN);
        PicturePlot(2, "Sqrt(<z2>)", 3, t, RDz.r2recent, 0, (int)(maxx*0.33), maxy2, maxy, first_time, LIGHTGREEN);
      }
      //if(t*Nmeasure%Niter == 0) {
        //PictureDistribution(1, 2, maxx*0.66, maxx, maxy2, maxy, first_time, GREEN);
        //PictureDistribution(1, 2, (int)(maxx*0.55), maxx, maxy2, maxy, first_time, GREEN);
        //PictureDistribution(2, 4, 0, (int)(maxx*0.55), maxy2, maxy, first_time, GREEN);
      //}
      if(t*Nmeasure%Niter == 0) {
        //PictureCoordinates('x', 'y', (int)(maxx*0.33), (int)(maxx-1), (int)(maxy2), (int)(maxy-1));
        PictureCoordinates('x', 'y', (int)(maxx*0.66), (int)(maxx-1), (int)(maxy2), (int)(maxy-1));
      }
    }
    else if(MC == DIFFUSION){ // i.e. DMC
      PicturePlot(0, "Energy", 1, t, E, 0, maxx, 40, maxy2, first_time, LIGHTGREEN);
      PicturePlot(0, "", 1, t, Evar/(DOUBLE)Nmeasured, 0, maxx, 40, maxy2, OFF, YELLOW);
      PicturePlot(1, "Population of walkers", 2, t, Nwalkersw, 0, (int)(maxx*0.66), maxy2, maxy, OFF, LIGHTBLUE);
      PicturePlot(1, "Population of walkers", 2, t, (DOUBLE) Nwalkers, 0, (int)(maxx*0.66), maxy2, maxy, first_time, YELLOW);
      if(t*Nmeasure%Niter == 0) {
        //PictureCoordinates('z', 'y' ,0, maxx*0.32, maxy-0.33*maxx, maxy-1);
        //PictureCoordinates('z', 'x', maxx*0.33, maxx*0.64, maxy-0.33*maxx, maxy-1);
        PictureCoordinates('y', 'x', (int)(maxx*0.67), (int)(maxx-1), (int)(maxy-0.33*maxx), (int)(maxy-1));
      }
    }
    else if (MC == CLASSICAL) {
      // plot energy
      PicturePlot(0, "Energy", 1, t, E, 0, (int)(maxx*0.66), 0, maxy2, first_time, LIGHTGREEN);
      if(t*Nmeasure%Niter == 0) {
        PictureCoordinates('x', 'y', (int)(maxx*0.), (int)(maxx-1), (int)(maxy2), (int)(maxy-1));
      }
      if(measure_R2) {
        PicturePlot(2, "<z2>", 3, t, RDz.r2recent, (int) (maxx*0.66), maxx, 0,  maxy2, first_time, LIGHTGREEN);
      }
    }
  }

  if(first_time) first_time = OFF;
  t++;
}

/********************************** Picture *********************************/
void PicturePlot(int plotn, char *title, int plottype, int t, DOUBLE f, int minxi, int maxxi, int minyi, int maxyi, int initialize, int colour) {
  static int height=11, width=9;
  static DOUBLE xscale[10], yscale[10];
  int i;
  char text[100];
  int exponent = 0;
  DOUBLE F, Fmax;
  int Fint;
  int minx, maxx, miny, maxy;
  int xc, yc;

  minx = minxi+5*width;
  maxy = maxyi-2*height;
  maxx = maxxi-width;
  miny = minyi+2*height;

  if(initialize) {
    setcolor(250);
    bar(minxi, minyi, maxxi, maxyi);
    setcolor(15);
    line(minx, maxy, maxx, maxy);
    line(minx, miny, minx, maxy);
    line(minx, miny, minx-3, miny+15);
    line(minx, miny, minx+3, miny+15);
    line(maxx, maxy, maxx-15, maxy-3);
    line(maxx, maxy, maxx-15, maxy+3);

    F = f;
    if(F == 0) F = 1.;

    while(fabs(F)<1. && F != 0.) {
      exponent--;
      F *= 10.;
    }
    while(fabs(F)>10.) {
      exponent++;
      F *= 0.1;
    }

    Fint = (int)(2*F);
    Fmax = Fint*pow(10, exponent);

    yscale[plotn] = (maxy-miny)/Fmax;
    xscale[plotn] = (DOUBLE)((maxx-minx)*Nmeasure)/(DOUBLE) (blck*Niter);

    _setcharsize(height, width);
    setcolor(50);
    _grtext(minx-20, miny-15, "x10");
    sprintf(text, "%i", exponent);
    _setcharsize(height/2, width/2);
    _grtext(minx-20+width*3, miny-15, text);
    _setcharsize(height, width);

    setcolor(WHITE);
    _setcharsize(1.2*height, 2*width);
    outtextxy(minx+40, miny-15, title);
    _setcharsize(height, width);

    for(i=0; i<10; i++) {
      setcolor(WHITE);
      F = Fmax/10.* (DOUBLE) i;
      line(minx-5, maxy-(int)(F*yscale[plotn]), minx+5, maxy-(int)(F*yscale[plotn]));
      _setlinestyle(DASHED);
      setcolor(150);
      if(i) line(minx+5, maxy-(int)(F*yscale[plotn]), maxx, maxy-(int)(F*yscale[plotn]));
      _setlinestyle(SOLID);
      setcolor(LIGHTCYAN);
      sprintf(text, "%.2" LF "", (DOUBLE) Fint * 0.1*(DOUBLE) i);
      outtextxy(minx-5*width, maxy-F*yscale[plotn]-height/2, text);
    }
    for(i=1; i<=blck; i++) {
      line(minx+(int)(xscale[plotn]*i*Niter/Nmeasure), maxy-4, minx+(int)(xscale[plotn]*i*Niter/Nmeasure), maxy+4);
      _setlinestyle(DASHED);
      setcolor(150);
      line(minx+(int)(xscale[plotn]*i*Niter/Nmeasure), maxy-4, minx+(int)(xscale[plotn]*i*Niter/Nmeasure), miny);
      _setlinestyle(SOLID);
      setcolor(LIGHTCYAN);
      sprintf(text, "%i", i);
      outtextxy(minx+(int)(xscale[plotn]*i*Niter/Nmeasure-width/2), maxy+(int)(4+height/2), text);
    }

    if(plottype == 1) { // Energy
#ifdef BC_3DPBC   // infinite system
      // plot Mean-field energy
      yc = (int) (maxy-(4.*PI*n*N)*yscale[plotn]);
      setcolor(LIGHTCYAN);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
#else // a trap and a tube
      // plot Mean-fielf energy
#ifdef TRIAL_3D // Energy > 1
      yc = (int) (maxy-(1+0.3*pow(3*(N-1)*lambda*a, 2./3.))*yscale[plotn]);
#else           // subtract 1 from the energy
      yc = (int) (maxy-0.3*pow(3*(N-1)*lambda*a, 2./3.)*yscale[plotn]);
#endif
      setcolor(LIGHTCYAN);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
      // plot Tonks energy
#ifdef TRIAL_3D // Energy > 1
      yc = (int) (maxy-(1+N*lambda*0.5)*yscale[plotn]);
#else           // subtract 1 from the energy
      yc = (int) (maxy-N*lambda*0.5*yscale[plotn]);
#endif
      setcolor(LIGHTMAGENTA);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
#endif
    }
    else if(plottype == 2) { // Size of the population
      // plot the control borders of the population size
      setcolor(CYAN);
      yc = (int) (maxy-Npop_max*yscale[plotn]);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
      yc = (int) (maxy-Npop_min*yscale[plotn]);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
    }
    else if(plottype == 3) { // Mean <z2>
#ifdef BC_ABSENT
      // Mean-field <z2>
      yc = (int) (maxy-pow(3*(N-1)*lambda*a, 1./3.)/Sqrt(5)/lambda*yscale[plotn]);
      // Tonks <z2>
      yc = (int) (maxy-Sqrt(0.5*N/lambda)*yscale[plotn]);
      setcolor(LIGHTCYAN);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
      setcolor(LIGHTMAGENTA);
      if(yc>miny && yc<maxy) line(minx, yc, maxx, yc);
#endif
    }
  }

  xc = (int) (minx+t*xscale[plotn]);
  yc = (int) (maxy-f*yscale[plotn]);
  if(xc>minx && xc<maxx && yc>miny && yc<maxy) 
    putpixel(xc, yc, colour);
  else {
    if(xc<minx) xc = minx;
    if(xc>maxx) xc = maxx;
    if(yc<miny) yc = miny;
    if(yc>maxy) yc = maxy;
    putpixel(xc, yc, YELLOW);
  }
}

/********************************** Picture *********************************/
void PictureCoordinates(char absciss, char ordinate, int minx, int maxx, int miny, int maxy) {
  int maxx2, maxy2;
  int w,i;
#ifdef TRIAL_2D
  DOUBLE maximum = 0;
#endif
#ifndef TRIAL_1D
  DOUBLE yscale;
#endif
  DOUBLE xscale;
  int radius = 5;
  DOUBLE x,y;
  int xc, yc;
  char ch[2] = " ";

  setcolor(LIGHTCYAN);
  line(minx, miny, maxx, miny);
  line(minx, maxy, maxx, maxy);
  line(maxx, miny, maxx, maxy);
  line(minx, miny, minx, maxy);
  setcolor(0);
  bar(minx+1, miny+1, maxx-1, maxy-1);
  minx += 2*radius;
  miny += 2*radius;
  maxx -= 2*radius;
  maxy -= 2*radius;

  maxx2 = (minx+maxx)/2;
  maxy2 = (miny+maxy)/2;

#ifdef TRIAL_2D  // zig-zag
  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      if(fabs(W[w].y[i]) > maximum) maximum = fabs(W[w].y[i]);
    }
  }
  maximum *= 4;
#ifdef  BC_1DPBC_X
  xscale = (maxx-minx)/Lx;
  yscale = (maxy-miny)/maximum;

  setcolor(15);
  line(minx, maxy2, maxx, maxy2);
  line(maxx2, miny, maxx2, maxy);

  for(i=-10; i<10; i++) {
    //line(minx+i*xscale*0.1, maxy2-5, minx+i*xscale*0.1, maxy2+5);
    //line(maxx2-5, maxy2-i*xscale*0.1, maxx2+5, maxy2-i*xscale*0.1);
    if(0.1*i<maximum/2) line(maxx2-5, maxy2-(int)(0.1*i*yscale), maxx2+5,  maxy2-(int)(0.1*i*yscale));
  }
#else //BC_2DPBC
  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      if((absciss == 'x' || ordinate == 'x') && fabs(W[w].x[i]) > maximum) maximum = fabs(W[w].x[i]);
      if((absciss == 'y' || ordinate == 'y') && fabs(W[w].y[i]) > maximum) maximum = fabs(W[w].y[i]);
      if((absciss == 'z' || ordinate == 'z') && fabs(W[w].z[i]) > maximum) maximum = fabs(W[w].z[i]);
    }
  }
  //maximum *= 2; //trap
  xscale = (maxx-minx)/maximum;
  yscale = (maxy-miny)/maximum;

  setcolor(15);
  line(minx, maxy2, maxx, maxy2);
  line(maxx2, miny, maxx2, maxy);

  ellipse(maxx2, maxy2, 0, 0, xscale, yscale);
  for(i=-10; i<10; i++) {
    line(maxx2-i*xscale*0.1, maxy2-5, maxx2-i*xscale*0.1, maxy2+5);
    line(maxx2-5, maxy2-i*xscale*0.1, maxx2+5, maxy2-i*xscale*0.1);
  }
#endif
#endif


#ifdef TRIAL_3D
  for(i=0; i<10; i++) {
    line(minx+(int)(i*xscale*0.1), maxy, minx+(int)(i*xscale*0.1), maxy-5);
    line(minx, maxy-(int)(i*xscale*0.1), minx+5, maxy-(int)(i*xscale*0.1));
  }
  xscale = (maxx-minx)/L;
  yscale = (maxy-miny)/L;
#endif

#ifdef TRIAL_1D
  xscale = (maxx-minx)/L;
#endif

  setcolor(LIGHTCYAN);

  ch[0] = absciss;
  outtextxy(maxx, maxy, ch);
  ch[0] = ordinate;
  outtextxy(minx-radius, miny, ch);

  for(w=0; w<Nwalkers; w++) {
    for(i=0; i<N; i++) {
      if(absciss == 'x') 
        x = W[w].x[i];
      else if(absciss == 'y') 
        x = W[w].y[i];
      else
        x = W[w].z[i];
      if(ordinate == 'x') 
        y = W[w].x[i];
      else if(ordinate == 'y')
        y = W[w].y[i];
      else
        y = W[w].z[i];

#ifdef TRIAL_3D
      xc = minx+(int)(x*xscale);
      yc = maxy-(int)(y*yscale);
#endif

#ifdef TRIAL_2D  // zig-zag
#ifdef BC_1DPBC_X
      xc = minx+(int)(x*xscale);
      yc = maxy2-(int)(y*yscale);
#else
      xc = minx+(int)(x*xscale); //hom
      yc = maxy-(int)(y*yscale);
#endif
#endif

#ifdef TRIAL_1D
      xc = minx+(int)(x*xscale);
      yc = 0;
#endif

      setcolor(1+i);
      if(xc>minx && xc<maxx && yc>miny && yc<maxy) 
        circle(xc, yc, radius);
    }
  }
}

/********************************** Picture Distribution ********************/
//plot type:
//  0) general
//  1) radial distribution function
//  2) pair distribution function
//  3) one body density matrix
//  4) z distribution
void PictureDistribution(int plotn, int plottype, int minxi, int maxxi, int minyi, int maxyi, int initialize, int colour) {
  static int height=11, width=9;
  static DOUBLE xscale[10], yscale[10];
  int i, size;
  char text[100];
  DOUBLE F;
  static int minx, maxx, miny, maxy;
  int xc, yc, xold, yold;
  DOUBLE dr, r, f, rold, fold, normalization, rmin;

  minx = minxi+5*width;
  maxy = maxyi-2*height;
  maxx = maxxi-width;
  miny = minyi+2*height;

  setcolor(250);
  //bar(minxi, minyi, maxxi, maxyi);
  setcolor(15);
  line(minx, maxy, maxx, maxy);
  line(minx, miny, minx, maxy);
  line(minx, miny, minx-3, miny+15);
  line(minx, miny, minx+3, miny+15);
  line(maxx, maxy, maxx-15, maxy-3);
  line(maxx, maxy, maxx-15, maxy+3);

  yscale[plotn] = maxy-miny;
  if(plottype == 1) {
    xscale[plotn] = (DOUBLE)(maxx-minx)/ (RD.max-RD.min);
    dr = 0.1*(RD.max-RD.min);
    rmin = RD.min;
  }
  else if(plottype == 2) {
    xscale[plotn] = (DOUBLE)(maxx-minx)/ (PD.max-PD.min);
    dr = 0.1*(PD.max-PD.min);
    rmin = PD.min;
  }
  else if(plottype == 3) {
    xscale[plotn] = (DOUBLE)(maxx-minx)/ (OBDM.max-OBDM.min);
    dr = 0.1*(OBDM.max-OBDM.min);
    rmin = OBDM.min;
  }
  else if(plottype == 4) {
    xscale[plotn] = (DOUBLE)(maxx-minx)/ (RDz.max-RDz.min);
    dr = 0.1*(RDz.max-RDz.min);
    rmin = RDz.min;
  }

  for(i=0; i<10; i++) {
    setcolor(WHITE);
    F = 0.1 * (DOUBLE) i;
    line(minx-5, maxy-(int)(F*yscale[plotn]), minx+5, maxy-(int)(F*yscale[plotn]));
    _setlinestyle(DASHED);
    setcolor(150);
    if(i) line(minx+5, maxy-(int)(F*yscale[plotn]), maxx, maxy-(int)(F*yscale[plotn]));
    _setlinestyle(SOLID);
    setcolor(LIGHTCYAN);
    sprintf(text, "%.2" LF "", 0.1*(DOUBLE) i);
    outtextxy(minx-5*width, maxy-F*yscale[plotn]-height/2, text);
  }
  for(i=1; i<10; i++) {
    line(minx+(int)(xscale[plotn]*(rmin+(DOUBLE)i*dr)), maxy-4, minx+(int)(xscale[plotn]*(rmin+(DOUBLE)i*dr)), maxy+4);
    _setlinestyle(DASHED);
    setcolor(150);
    line(minx+(int)(xscale[plotn]*(rmin+(DOUBLE)i*dr)), maxy-4, minx+(int)(xscale[plotn]*(rmin+(DOUBLE)i*dr)), miny);
    _setlinestyle(SOLID);
    setcolor(LIGHTCYAN);
    sprintf(text, "%.2" LF "", i*dr+rmin);
    outtextxy(minx+xscale[plotn]*(rmin+(DOUBLE)i*dr)-width/2, maxy+4+height/2, text);
  }

  if(!initialize) {
    if(plottype == 1) {
      rold = 0;
      fold = 0.;
      dr = 1./ (4.*PI*RD.step * RD.width *(DOUBLE) (RD.times_measured * N));
      size = RD.size;
    }
    else if(plottype == 2) {
      rold = 0;
      fold = 0.;
      size = PD.size-1;

#ifdef BC_ABSENT
      normalization = (1+1./N)/ ((DOUBLE)(PD.times_measured*N)*PD.step);
#else
      normalization = 2./ ((DOUBLE)(PD.times_measured*N)*n*4*PI*PD.step);
#endif

    }
    else if(plottype == 3) {
      rold = 0;
      fold = 1.;
      size = OBDM.size;
    }
    else if(plottype == 4) { // RDz
      rold = 0;
      fold = 0.;
      dr = 1./(2.*RDz.step*(DOUBLE) (RDz.times_measured * N));
      size = RD.size;
    }
    xold = (int) (minx+rold*xscale[plotn]);
    yold = (int) (maxy-fold*yscale[plotn]);

    for(i=0; i<size; i++) {
      if(plottype == 1) {
        r = RD.min + (DOUBLE) (i+0.5) * RD.step;
        f = Sqrt(dr*(DOUBLE)RD.N[i]/r);
      }
      else if(plottype == 2) {
        r = PD.min + (DOUBLE) (i+1) * PD.step;
#ifdef BC_ABSENT
        f = normalization*(DOUBLE)PD.N[i];
#else
        f = normalization*(DOUBLE)PD.N[i]/(r*r);
#endif

      }
      else if(plottype == 3) {
        r = (DOUBLE) (i) * OBDM.step;
        if(OBDM.N[i] > 0) {
          f = OBDM.f[i]/(DOUBLE) OBDM.N[i];
        }
      }
      else if(plottype == 4) { // RDz
        r = (DOUBLE) (i+0.5) * RDz.step;
        f = dr*(DOUBLE)RDz.N[i];

#ifdef TRIAL_1D // 1D
#ifdef BC_ABSENT // trap
       normalization = 1./(2*RDz.step*(DOUBLE) (RDz.times_measured * N));
       f = normalization*(DOUBLE)RDz.N[i];
#else // 1D PBC
       normalization = 1./(RDz.step*(DOUBLE) (RDz.times_measured * N));
       f = normalization*(DOUBLE)RDz.N[i]*N/n;
#endif
#endif
      }
      xc = (int) (minx+r*xscale[plotn]);
      yc = (int) (maxy-f*yscale[plotn]);
      if(xc<minx) xc = minx;
      if(xc>maxx) xc = maxx;
      if(yc<miny) yc = miny;
      if(yc>maxy) yc = maxy;

      if(i) line(xc, yc, xold, yold);

      xold = xc;
      yold = yc;
    }
  }
}

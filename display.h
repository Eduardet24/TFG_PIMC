/*display.h*/
#ifndef _DISPLAY_H_
#define _DISPLAY_H_

void CloseGraph(void);
void Picture(void);
void PictureCoordinates(char absciss, char ordinate, int minx, int maxx, 
  int miny, int maxy) ;
void PicturePlot(int plotn, char *title, int picturetype, int t, DOUBLE f, int minxi, 
  int maxxi, int minyi, int maxyi, int initialize, int colour);
void PictureDistribution(int plotn, int plottype, int minxi, int maxxi, 
  int minyi, int maxyi, int initialize, int colour);

#endif

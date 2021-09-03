/*********************************************************************
 * libconf: A small library to parse config files                    *
 *                                                                   *
 * Author: Andreas Heiss                                             *
 * Date  : March 2000                                                *
 *                                                                   *
 *********************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "libconf.h"

static int conf_buf_pointer;
static char conf_fn[500];
static char conf_buffer[500][100];

int conf_read(char *fn)
{
  FILE *cfile;
  char buf[100];

  cfile = fopen(fn,"r");
  if (cfile == NULL) return(1);
  
  conf_buf_pointer = -1;
  while(1) {
    if (fgets(buf,100,cfile) == NULL) break;
    if (! strncmp(buf,"#",1)) continue; /* ignore comments */
    while(buf[0]==32) strcpy(buf,&buf[1]); /* get rid of leading blanks */
    if (strlen(buf)<2) continue;
    strcpy(conf_buffer[++conf_buf_pointer],buf);
  }
  fclose(cfile);
       
  return(0);
}

int conf_readkeyword(char *name, char *fmt, ...)
{
  int i,l,ok,*d;
  char lbuf[100],*s;
  float *f;
  double *lf;
  va_list ap;
  
  l = strlen(name);
  ok = 0;
  for (i=0; i<=conf_buf_pointer; i++) {
    if (!strncmp(name,conf_buffer[i],l)) {
      ok = 1;
      strcpy(lbuf,conf_buffer[i]);
      while (lbuf[0]!=32)  strcpy(lbuf,&lbuf[1]); 
      while (lbuf[0]==32)  strcpy(lbuf,&lbuf[1]); 
      strcat(lbuf," ");
    }
  }
  if (!ok) return(1);

  va_start(ap, fmt);
  while(*fmt) {
    switch(*fmt++) {
    case 's':           /* string */
      s = va_arg(ap, char *);
      sscanf(lbuf,"%s",s);
      while (lbuf[0]!=32)  strcpy(lbuf,&lbuf[1]); 
      while (lbuf[0]==32)  strcpy(lbuf,&lbuf[1]); 
      break;
    case 'd':           /* int */
      d = va_arg(ap, int *);
      sscanf(lbuf,"%d",d);
      while (lbuf[0]!=32)  strcpy(lbuf,&lbuf[1]); 
      while (lbuf[0]==32)  strcpy(lbuf,&lbuf[1]); 
      break;
    case 'f':           /* float */
      f = va_arg(ap, float *);
      sscanf(lbuf,"%f",f);
      while (lbuf[0]!=32)  strcpy(lbuf,&lbuf[1]); 
      while (lbuf[0]==32)  strcpy(lbuf,&lbuf[1]); 
      break;
    case 'F':
      lf = va_arg(ap, double *);
      sscanf(lbuf,"%lf",lf);
      while (lbuf[0]!=32)  strcpy(lbuf,&lbuf[1]); 
      while (lbuf[0]==32)  strcpy(lbuf,&lbuf[1]); 
      break;
    }
  }
  va_end(ap);
  return(0);
}

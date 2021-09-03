/**********************************************
 * A simple example program using 'libconf.a' *
 **********************************************/

#include <stdio>
#include "libconf.h"

/* If you are using a c++ compile, do something like 
extern "C" {   
#include "libconf.h"
}
*/

main()
{
  int ierr;
  int a;
  float p[5];
  char string[100];
  double r;

  /* Read the config file */
  ierr = conf_read("example.conf");
  if (ierr != 0) {
    printf("ERROR reading config file !\n");
    exit(0);
  }
  
  /* Read information into variables */
  ierr = conf_readkeyword("A_NOT_DEFINED_KEYWORD","%d",&a);
  if (ierr) {
    printf("Keyword A_NOT_DEFINED_KEYWORD not found !\n");
  }
  
  /* Important: The format descriptor for 'double' is 'F' ! */
  ierr = conf_readkeyword("OTHERPAR","%d %s %F",&a,string,&r);
  if (!ierr) 
    printf("Other parameters: %d %s %21.17lf\n",a,string,r);

  ierr = conf_readkeyword("STARTPAR","%f %f %f %f %f",&p[0],&p[1],&p[2],&p[3],&p[4]);
  if (!ierr) 
    printf("Startpar = %f %f %f %f %f\n",p[0],p[1],p[2],p[3],p[4]);

}

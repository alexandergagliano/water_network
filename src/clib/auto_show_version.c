#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"The Grackle Version 3.1\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision 9060cbe0f2dc4a5047cbaf5e184c6217fc78963b\n");
   fprintf (fp,"\n");
}

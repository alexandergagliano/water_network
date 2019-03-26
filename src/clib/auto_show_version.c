#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"The Grackle Version 3.1\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision 40684201766f1f83b5312497fae0b2ee577d029b\n");
   fprintf (fp,"\n");
}

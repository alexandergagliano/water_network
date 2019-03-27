#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"The Grackle Version 3.1\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision ac7125ec2b32271b6160be620770f4a5a418901b\n");
   fprintf (fp,"\n");
}

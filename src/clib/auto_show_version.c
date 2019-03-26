#include <stdio.h>
void auto_show_version(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"The Grackle Version 3.1\n");
   fprintf (fp,"Git Branch   master\n");
   fprintf (fp,"Git Revision 9827f499045d29f29d5ef08811083fe3ddd2ab81\n");
   fprintf (fp,"\n");
}

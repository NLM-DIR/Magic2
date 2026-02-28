/*    dirsync_nfs.c   -- a c-program which will create/lock/unlock/remove a file in a directory
                         this will effectively flush the inode-cache for that directory only.
                         it allows the next stat of a file to be directly from the NFS server.
                         Very useful for trigger or controll files to coordinate multiple concurrent jobs in a cluster.
                         it prevents the otherwise necessary need to have an abusively low actimeo/acdirmax/acregmax mount options.
     This example file was found on the net more than 15 years ago, original author not credited only because I can't find them.
     the same technique works in python, perl, tcl, and even Bash with flock !
*/
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>


#define MAXSTR 1024
#define FILENAME "/@$%@#%@@"

int lock_unlock(int fd, char *path);

int main(int argc, char *argv[])
{
  int status;
  char filepath[MAXSTR];
  struct stat statbuf;
  int fd;

  /* must be 1 argument - a regular or directory file */
  if (argc == 2)
     strcpy(filepath, argv[1]);
  else {
    fprintf(stderr, "usage: %s <file> \n", argv[0]);
    exit(1);
  }

  /* make sure file/directory exists and stat for mode */
  if ((status= stat(filepath, &statbuf)) < 0){
      perror("stat");
      fprintf(stderr, "cannot stat %s\n", filepath);
      exit(1);
  }

 /* IF the argument is a directory, we create/lock/unlock/remove a file in it */
  if ( S_ISDIR(statbuf.st_mode) ) {
      strcat(filepath, FILENAME);
     /* CREATE the file */
      if((fd=open(filepath, O_RDWR | O_CREAT, 0777)) < 0) {
         fprintf(stderr, "could not create %s \n", filepath);
         exit(1);
      }
      lock_unlock(fd, filepath);  /* whether it works or fails, we unlink the fisle*/
      unlink(filepath);
  }

 /* IF the argument is a regular file, we lock/unlock it */
  if ( S_ISREG(statbuf.st_mode) ) {
      if((fd=open(filepath, O_RDWR)) < 0) {
         fprintf(stderr, "could not open %s \n", filepath);
         exit(1);
      }
      lock_unlock(fd,filepath);
  }
 /* IF the argument is not a regular or dir file, silently exit */
}

int lock_unlock(int fd, char *path)
{
  struct flock flockbuf;
  int status;

  flockbuf.l_type=F_WRLCK;
  flockbuf.l_whence=0;
  flockbuf.l_start=0;
  flockbuf.l_len=1;

  status = fcntl(fd, F_SETLK, &flockbuf);
  if (status < 0) {
       perror("fcntl");
       fprintf(stderr, "could not lock %s \n", path);
       return(-1);
  }
  status = fcntl(fd, F_SETLKW, &flockbuf);
  if (status < 0) {
       perror("fcntl");
       fprintf(stderr, "could not unlock %s \n", path);
       return(-1);
  }
  return(0);
}

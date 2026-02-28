/*  
 * misc files from Ace-Conn library concatenated together and hacked up so
 * as to not depend on libraries I don't have.  In the course of 
 * importing it, I also hacked it up to clean out features I don't use 
 * and otherwise clarify the code (to the extent that it is possible to do
 * that while concatenating somebody else's library into a single file).
*
* THE CURRENT STATE OF THIS CODE:  I have used this code in the wac
* library to talk to the Sanger Acedb 4.9 socket server.  I have not
* confirmed that I did not create any memory leaks or similar bugs.
*
 * 
 * Mark S. 3/2003
 *
 * Original information on this code is:
 *
 *  Author: Ed Griffiths (edgrif@sanger.ac.uk)
 *  Copyright (c) Sanger Institute, 2002
 *-------------------------------------------------------------------
 * Acedb is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * or see the on-line version at http://www.gnu.org/copyleft/gpl.txt
 *-------------------------------------------------------------------
 * Description: Internal utility functions for AceConn package.
 *              
 * Exported functions: None
 * HISTORY:
 * Last edited: Aug 27 20:22 2002 (edgrif)
 * Created: Thu Mar  7 09:34:46 2002 (edgrif)
 * CVS info:   $Id: acclient_socket_lib.c,v 1.14 2016/01/26 04:02:42 mieg Exp $
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#if defined(NEXT) || defined(HP)
extern void* malloc (mysize_t size) ;
#elif !defined(WIN32)  && !defined(MACINTOSH) && ! defined(MAC_X)
#include <malloc.h>   /* normal machines  */
#endif

#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>

#include "../wmd5/global.h"					    /* MD5 */
#include "../wmd5/md5.h"					    /* MD5 */


/* Must be set in header so server can detect clients byte order.            */
#define ACECONN_BYTE_ORDER 0x12345678

/* We strings rather than enums, easier for others to interface to...        */
/* No ACESERV_MSGnnn string should be longer than                            */
/*                             (ACECONN_MSGTYPE_BUFLEN - 1)                  */
/* you must not lightly change this length, this is the size clients will    */
/* be expecting, to change it means changing all clients.                    */
/*                                                                           */
enum {ACECONN_MSGTYPE_BUFLEN = 30} ;

/* Client only, its messages are either normal requests (commands) or are    */
/* ace data that needs to be parsed in.                                      */
#define ACECONN_MSGREQ    "ACESERV_MSGREQ"
#define ACECONN_MSGDATA   "ACESERV_MSGDATA"

/* Server only, it may just be sending or a reply or it may be sending an    */
/* instruction, such as "operation refused".                                 */
#define ACECONN_MSGOK     "ACESERV_MSGOK"
#define ACECONN_MSGENCORE "ACESERV_MSGENCORE"
#define ACECONN_MSGFAIL   "ACESERV_MSGFAIL"
#define ACECONN_MSGKILL   "ACESERV_MSGKILL"


/* The message header struct, this has to be packed into the the first       */
/* ACECONN_HEADER_BYTES of a message sent between client & server.           */
/* You can use these enums to help you in writing code to put bytes into the */
/* message.                                                                  */

/* Number of fields in the header that may need to be byte swopped.          */
enum {ACECONN_HEADER_SWAPFIELDS = 5} ;

/* Because we want a CONSTANT size for the buffer, we define this very       */
/* tediously...we have 5 ints of 4 bytes each and a char buffer to hold the  */
/* message type which is ACECONN_MSGTYPE_MAXLEN long.                        */
enum {ACECONN_HEADER_BYTES = ((ACECONN_HEADER_SWAPFIELDS * 4) + ACECONN_MSGTYPE_BUFLEN)} ;

typedef struct AceConnHeaderRecStruct
{
  int byte_swap ;					    /* Set to ACECONN_BYTE_ORDER */
  int length ;						    /* Length of data being sent in bytes. */
  int server_version ;					    /* Server sets this. */
  int client_id ;					    /* Unique client identifier, set by
							       server. */
  int max_bytes ;					    /* Maximum request size, set by server. */
  char msg_type[ACECONN_MSGTYPE_BUFLEN] ;		    /* See the msgs defined above. */
} AceConnHeaderRec, *AceConnHeader ;


/* On connection to the server, the client must send this text.              */
#define ACECONN_CLIENT_HELLO "bonjour"

/* The server will reply with this string when password verification is      */
/* complete.                                                                 */
#define ACECONN_SERVER_HELLO "et bonjour a vous"

/* Server can send a reply in slices, it will return "encore" when this is   */
/* true and client must reply with "encore" to get further slices.           */
#define ACECONN_SERVER_CLIENT_SLICE "encore"

/* To disconnect the server must sent this text.                             */
#define ACECONN_CLIENT_GOODBYE "quit"


/* Opaque handle to a connection to an Acedb socket server.                  */
typedef struct _AceConnRec *AceConnection ;

/* Return codes for the AceConn calls.                                       */
typedef enum {ACECONN_OK = 0,
	      ACECONN_QUIT,				    /* Command to server was "quit" so
							       connection is now closed. */
	      ACECONN_INVALIDCONN,			    /* Connection points to invalid memory. */
	      ACECONN_BADARGS,				    /* caller has supplied bad args. */
	      ACECONN_NOALLOC,				    /* could not allocate  */
	      ACECONN_NOSOCK,				    /* socket creation problem. */
	      ACECONN_UNKNOWNHOST,			    /* server host not known. */
	      ACECONN_NOCONNECT,			    /* Could not connect to host. */
	      ACECONN_NOSELECT,				    /* select on socket failed. */
	      ACECONN_HANDSHAKE,			    /* Handshake to server failed. */
	      ACECONN_READERROR,			    /* Error in reading from socket. */
	      ACECONN_WRITEERROR,			    /* Error in writing to socket. */
	      ACECONN_SIGSET,				    /* Problem with signal setting. */
	      ACECONN_NONBLOCK,				    /* Non-blocking for socket failed. */
	      ACECONN_TIMEDOUT,				    /* Connection timed out. */
	      ACECONN_NOCREATE,				    /* Could note create connection
							       control block. */
	      ACECONN_INTERNAL				    /* Dire, internal package error. */
} AceConnStatus ;


/* Default time (seconds) to wait for reply from server. */
enum {ACECONN_DEFAULT_TIMEOUT = 120} ;


/* Creates a connection struct. */
static AceConnStatus AceConnCreate(AceConnection *connection, const char *server_netid, int server_port,
			    const char *userid, const char *passwd, int timeout) ;


/* Open a connection to the database. */
static AceConnStatus AceConnConnect(AceConnection connection) ;

/* Send a request to the database and get the reply. */
static AceConnStatus AceConnRequest(AceConnection connection,
			     const char *request,
			     void **reply, int *reply_len) ;

/* Close down the connection. */
static AceConnStatus AceConnDisconnect(AceConnection connection) ;


/* Free the connection struct. */
static void AceConnDestroy(AceConnection connection) ;



/* Utility functions.                                                        */


/* Opaque handle to an message passed to server.                             */
typedef struct AceConnMsgRec *AceConnMsg ;

/* Connection can be in one of a number of states.                           */
typedef enum
{
  ACECONNSTATE_NONE,					    /* not connected at all to server. */
  ACECONNSTATE_OPEN,					    /* socket open. */
  ACECONNSTATE_READY					    /* connection ready for requests. */
} AceConnState ;


/* Holds all the data for a particular connection to the ace server.         */
typedef struct _AceConnRec
{
  char **valid ;					    /* Used for checking valid block. */

  AceConnState state ;

  /* Server details                                                          */
  const char *host ;
  int  port ;
  int  socket ;
  int  timeout ;					    /* command timeout in secs. */

  /* User details                                                            */
  const char *userid ;
  const char *passwd ;
  char *passwd_hash ;
  char *nonce_hash ;

  /* Current request.                                                        */
  AceConnMsg msg ;

  /* Misc.                                                                   */
  char *last_errmsg ;					    /* string describing last error. */

} AceConnRec ;


/* Basic creation/destruction of an AceConnRec.                              */
static AceConnection aceconnCreate(const char *server_netid, int server_port,
			    const char *userid, const char *passwd, int timeout) ;
static void aceconnDestroy(AceConnection connection) ;


/* Server requests.                                                          */
static AceConnStatus aceconnServerOpen(AceConnection connection) ;
static AceConnStatus aceconnServerHandshake(AceConnection connection) ;
static AceConnStatus aceconnServerRequest(AceConnection connection,
				   const char *request, void **reply, int *reply_len) ;
static AceConnStatus aceconnServerClose(AceConnection connection) ;


/* Utilities.                                                                */

/* This string is used as a unique symbol for checking that a connection     */
/* struct is valid. There can only be one of these pointers with this name   */
/* in this namespace so its address is unique and its this address that we   */
/* use for checking. Its most unlikely that any bit of random memory will    */
/* have this bit pattern. Using a string like this also helps in debugging   */
/* as its possible to see the string.                                        */
/* Note that the struct is initialised with this address and then when freed */
/* this part of the struct is cleared to zeros, this is vital because it is  */
/* possible for the memory that a dangling pointer addresses to remain       */
/* unaltered for some time and hence look like a valid pointer.              */

#define ACECONN_CONTROL_BLOCK_STRING "connection_control_block"

static char* CCB_ptr = ACECONN_CONTROL_BLOCK_STRING ;


/* Create a basic control block for an acedb server connection.              */
/*                                                                           */
static AceConnection aceconnCreate(const char *server_netid, int server_port,
			    const char *userid, const char *passwd, int timeout)
{
  AceConnection connection = NULL ;

  connection = (AceConnection)malloc(sizeof(AceConnRec)) ;
  if (connection)
    {
      connection->valid = &CCB_ptr ;
      connection->state = ACECONNSTATE_NONE ;
      connection->host = server_netid ;
      connection->port = server_port ;
      connection->socket = -1 ;
      if (timeout < 0)
	connection->timeout = ACECONN_DEFAULT_TIMEOUT ;
      else
	connection->timeout = timeout ;
      connection->userid = userid ;
      connection->passwd = passwd ;
      connection->passwd_hash = NULL ;
      connection->nonce_hash = NULL ;
      connection->msg = NULL ;
      connection->last_errmsg = NULL ;
    }

  return connection ;
}

/* Assumed never to fail.                                                    */
static void aceconnFreeErrMsg(AceConnection connection)
{
  if (connection->last_errmsg)
    free(connection->last_errmsg) ;

  return ;
}


static void aceconnDestroy(AceConnection connection)
{

  if (connection)
    {
      if (connection->passwd_hash)
	free(connection->passwd_hash) ;

      aceconnFreeErrMsg(connection) ;

      free(connection) ;
    }

  return ;
}


static void aceconnSetErrMsg(AceConnection connection, const char *msg1, const char *msg2 )
{
/*
* if the malloc fails, we crash.  not nice, but that was what we were
* going to do as soon as we return to the caller anyway.
*/
connection->last_errmsg = (char *) malloc(strlen(msg1) + strlen(msg2) + 2);
strcpy(connection->last_errmsg, msg1);
strcat(connection->last_errmsg, " ");
strcat(connection->last_errmsg, msg1);
}



/* Messages:                                                                 */
/*                                                                           */

struct cheapstring
	{
	unsigned char *str;
	int len;
	};

struct AceConnMsgRec
{
  /* header handling */
  AceConnHeaderRec proto_hdr ;
  int done_header ;
  char header_buf[ACECONN_HEADER_BYTES] ;
  int hBytesRequested ;
  int hBytesPending ;

  /* message handling */
  int done_message ;
  char *message ;
  struct cheapstring encore_message ;
  int mBytesRequested ;
  int mBytesPending ;
} ;


/* Messages are in one of these states while being sent, this is because even*/
/* without non-blocking I/O you can't guarantee sending the message in one   */
/* go if it is very big.                                                     */
typedef enum {ACECONNMSG_WAIT, ACECONNMSG_DONE, ACECONNMSG_ERROR} AceConnMsgState ;


/* Sockets can be accessed for READ or WRITE, not both currently.            */
typedef enum {ACECONNACC_READ, ACECONNACC_WRITE} AceConnSockAccess ;


/* The MD5 code returns an array of unsigned char of size 16, the value 16   */
/* has no symbolic constant in md5.h so I define one, plus the size of the   */
/* string required to hold the hexadecimal string version of the array.      */
/*                                                                           */
enum {MD5_HASHLEN = 16, MD5_HEX_HASHLEN = ((16 * 2) + 1)} ;


static AceConnStatus accessSocket(AceConnection connection, AceConnSockAccess access) ;
static AceConnStatus aceconnSocketRead(AceConnection connection) ;
static AceConnMsgState readSocket(int fd, char *readbuf, int *bytes_pending_inout,
				  char **errmsg_out) ;
static AceConnStatus aceconnSocketWrite(AceConnection connection) ;
static AceConnMsgState writeSocket(int fd, const char *writebuf, int *bytes_pending_inout,
				   char **errmsg_out) ;

static AceConnMsg messageCreate(void) ;
static void messageSetForRead(AceConnMsg msg, int encore) ;
static void messageSetForWrite(AceConnMsg msg, const char *message) ;
static void messageGetReply(AceConnMsg msg, void **reply, int *reply_len) ;
static int messageIsEncore(AceConnMsg msg_in) ;
static void messageGetEncore(AceConnMsg msg_in) ;
static int messageDone(AceConnMsg msg_in) ;
static void messageDestroy(AceConnMsg msg) ;

static void buf2hdr(char *buf, AceConnHeader hdr) ;
static void hdr2buf(AceConnHeader hdr, char *buf) ;
static void setMsgType(char buffer[], char *msgType) ;
static int testMsgType(char buffer[], char *msgType) ;

static AceConnStatus selectSocket(AceConnection connection, AceConnSockAccess access) ;
static AceConnStatus resetBlockingSocket(AceConnection connection, int *oldflags_out, int *newflags) ;
static AceConnStatus shutdownSocket(AceConnection connection) ;

static char *hashAndHexStrings(const char *strings[], int num_strings) ;
static char *makeHash(const char *userid, const char *passwd) ;
static char *convertMD5toHexStr(unsigned char digest[]) ;

/*                                                                           */
/* External package routines.                                                */
/*                                                                           */


/* Set up a socket connection to the server.                                 */
/*                                                                           */
static AceConnStatus aceconnServerOpen(AceConnection connection)
{ 
  AceConnStatus status = ACECONN_OK ;
  int sock ;
  struct sockaddr_in sname ;
  struct hostent *hp ;
  int n;

  /* create a socket for writing/reading */
  if (status == ACECONN_OK)
    {
      sock = socket(AF_INET, SOCK_STREAM, 0) ;
      if (sock < 0)
	{
	  status = ACECONN_NOSOCK ;
	  aceconnSetErrMsg(connection, "Failed to create socket: ", strerror(errno)) ;
	}
      else
	connection->socket = sock ;
    }

  n = 1;
  setsockopt( sock, 6, TCP_NODELAY , &n, sizeof(n));

  n=49152;
  setsockopt(sock, SOL_SOCKET, SO_SNDBUF, &n, sizeof(n));
  setsockopt(sock, SOL_SOCKET, SO_RCVBUF, &n, sizeof(n));


  if (status == ACECONN_OK)
    {
      /* create absolute socket name */
      memset((void *)&sname, 0, sizeof(sname)) ;
      sname.sin_family = AF_INET ;
      sname.sin_port = htons(connection->port) ;

      hp = gethostbyname(connection->host) ;
      if (!hp)
	{
	  status = ACECONN_UNKNOWNHOST ;
	  aceconnSetErrMsg(connection, "Unknown host: ", connection->host) ;
	}
    }


  /* OK, connect to the socket.                                              */
  if (status == ACECONN_OK)
    {
      memcpy((void*)&sname.sin_addr, (void*)hp->h_addr, hp->h_length) ;

      if (connect(sock, (struct sockaddr *)&sname, sizeof(sname)) == 0)
	connection->state = ACECONNSTATE_OPEN ;
      else
	{
	  status = ACECONN_NOCONNECT ;
	  aceconnSetErrMsg(connection, "Could not connect to socket: ", strerror(errno)) ;
	}
    }

  return status ;
}

/* Send an initial message to the client, we will then get a reply which     */
/* contains a 'nonce' or key with which we must encode the users userid and  */
/* password. Then get the users userid/password and do the encryption and    */
/* send it back to the server for verification. We should then get another   */
/* reply from the server to OK this.                                         */
/* If this routine fails then we exit.                                       */
/*                                                                           */
static AceConnStatus aceconnServerHandshake(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;
  char *reply = NULL ;
  int reply_len = 0 ;
  enum {HASH_STRING_NUM = 2} ;
  const char *hash_strings[HASH_STRING_NUM] ;
  char *hash_nonce ;
  char *request = NULL ;

  /* This is the first time we need to send a message to the server so       */
  /* create the message handling struct.                                     */
  connection->msg = messageCreate() ;

  /* We need to be able to add stuff to our error message to say we failed   */
  /* in handshake....                                                        */

  /* Say an initial hello to the server, the server returns a nonce string   */
  /* for us to hash our password hash with.                                  */
  if (status == ACECONN_OK)
    {
      status = aceconnServerRequest(connection, ACECONN_CLIENT_HELLO,
				    (void **)&reply, &reply_len) ;
      if (status == ACECONN_OK)
	{
	  /* WE SHOULD CHECK THAT THE SERVER ONLY RETURNED ONE WORD HERE....     */
	  /* where is the unix func to do that ??                                */

	  /* Create a hash of our userid and passwd and then hash this with the  */
	  /* nonce and convert to a hex string to pass back to server.           */
	  connection->passwd_hash = makeHash(connection->userid, connection->passwd) ;
	  hash_strings[0] = connection->passwd_hash ;
	  hash_strings[1] = reply ;
	  hash_nonce = hashAndHexStrings(hash_strings, HASH_STRING_NUM) ;
      
	  /* Now make a string containing the userid and hash to send back.  */
	  {
	  int n = strlen(connection->userid) + strlen(hash_nonce) + 2;
	  request = (char *) malloc(n);
  	  sprintf(request, "%s %s", connection->userid, hash_nonce) ;
	  }

	  free(reply) ;
	}
    }

  /* Send our nonce-hashed hash back, if all is ok the server will reply     */
  /* with its standard message.                                              */
  if (status == ACECONN_OK)
    {
      status = aceconnServerRequest(connection, request, (void **)&reply, &reply_len) ;
      if (status == ACECONN_OK)
	{
	  if (strcmp(reply, ACECONN_SERVER_HELLO) != 0)
	    {
	      status = ACECONN_HANDSHAKE ;
	      aceconnSetErrMsg(connection, reply, "") ;
	    }
	  else
	    { /* Return the userid and passwd hash for later reuse.   */ 
	      connection->nonce_hash = hash_nonce ;
	      connection->state = ACECONNSTATE_READY ;
	    }

	  free(reply) ;
	}
    }

  free(request) ;

  return status ;
}


/* Send a request to the server and get the reply.                           */
/*                                                                           */
static AceConnStatus aceconnServerRequest(AceConnection connection,
				   const char *request, void **reply, int *reply_len)
{
  AceConnStatus status = ACECONN_OK ;
  int encore = FALSE, finished = FALSE ;


  while (!finished)
    {
      if (status == ACECONN_OK)
	{
	  messageSetForWrite(connection->msg, request) ;

	  status = accessSocket(connection, ACECONNACC_WRITE) ;
	  if (status != ACECONN_OK)
	    finished = TRUE ;
	}
      
      if (status == ACECONN_OK)
	{
	  messageSetForRead(connection->msg, encore) ;

	  status = accessSocket(connection, ACECONNACC_READ) ;
	  if (status != ACECONN_OK)
	    finished = TRUE ;
	  else
	    {
	      if (messageIsEncore(connection->msg))
		{
		  encore = TRUE ;
		  request = ACECONN_SERVER_CLIENT_SLICE ;
		  messageGetEncore(connection->msg) ;
		}
	      else
		{
		  finished = TRUE ;
		  messageGetReply(connection->msg, reply, reply_len) ;
		}
	    }
	}
    }

  return status ;
}




static AceConnStatus aceconnServerClose(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;
  char *reply = NULL ;
  int reply_len = 0 ;

  if (connection->state == ACECONNSTATE_OPEN
      || connection->state == ACECONNSTATE_READY)
    {
      if (connection->state == ACECONNSTATE_READY)
	{
	  status = aceconnServerRequest(connection, ACECONN_CLIENT_GOODBYE,
					(void **)&reply, &reply_len) ;
	}

      /* Always close the socket.                                                */
      status =  shutdownSocket(connection) ;
    }

  messageDestroy(connection->msg) ;
  connection->msg = NULL ;
  connection->state = ACECONNSTATE_NONE ;

  return status ;
}



/*                                                                           */
/* Internal routines.                                                        */
/*                                                                           */


/* Read or write a message to a socket as a more or less atomic operation,   */
/* but allowing timeout via the use of select() to block for a specified time*/
/* if we have to wait on the socket.                                         */
/*                                                                           */
static AceConnStatus accessSocket(AceConnection connection, AceConnSockAccess access)
{
  AceConnStatus status = ACECONN_OK ;
  int finished ;

  finished = FALSE ;
  while (!finished)
    {
      /* Is the socket ready ?                                               */
      if ((status = selectSocket(connection, access)) != ACECONN_OK)
	finished = TRUE ;

      /* If so, process the rest of the message.                             */
      if (!finished)
	{
	  if (access == ACECONNACC_WRITE)
	    status = aceconnSocketWrite(connection) ;
	  else
	    status = aceconnSocketRead(connection) ;

	  /* We only finish if there is a problem (e.g. timeout) OR the msg  */
	  /* has been successfully processed.                                */
	  if (status != ACECONN_OK
	      || (status == ACECONN_OK && (messageDone(connection->msg))))
	    finished = TRUE ;
	}
    }

  return status ;
}


/**************************************************************/


static AceConnStatus aceconnSocketRead(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;
  AceConnMsgState state = ACECONNMSG_ERROR ;
  int sock_flags ;
  char *errmsg = NULL ;
  int fd = connection->socket ;
  AceConnMsg req = (AceConnMsg)connection->msg ;


  /* We can't just do a recv with MSG_DONTWAIT and MSG_PEEK because          */
  /* the former is not supported on many systems yet. So instead we make the */
  /* whole socket non-blocking, do the recv and then reset the socket.       */
  if (status == ACECONN_OK)
    {
      status = resetBlockingSocket(connection, &sock_flags, NULL) ;
    }

  if (status == ACECONN_OK)
    {
      /* Try to read the header, probably do this in one go.                     */
      if (req->done_header == FALSE)
	{
	  /* First time for this message so allocate a buffer.                   */
	  if (!req->hBytesPending)
	    {
	      req->hBytesPending = req->hBytesRequested = ACECONN_HEADER_BYTES ;
	    }

	  state = readSocket(fd, req->header_buf + (req->hBytesRequested - req->hBytesPending),
			     &(req->hBytesPending), &errmsg) ;

	  if (state == ACECONNMSG_DONE)
	    {
	      buf2hdr(req->header_buf, &(req->proto_hdr)) ;
	      req->done_header = TRUE ;
	    }
	  else if (state == ACECONNMSG_ERROR)
	    {
	      status = ACECONN_READERROR ;
	      aceconnSetErrMsg(connection, "Error reading message header from socket:", errmsg) ;
	      free(errmsg) ;
	    }
	}
    }

  /* We have a header so try to read buffer (this may take several calls to  */
  /* this function).                                                         */
  if (status == ACECONN_OK)
    {
      if (req->done_header == TRUE)
	{
	  /* First time for this message so allocate a buffer.                   */
	  if (!req->message)
	    {
	      req->message = (char *) malloc(req->proto_hdr.length) ;
	      req->mBytesRequested = req->mBytesPending = req->proto_hdr.length ;
	    }

	  state = readSocket(fd, req->message + req->mBytesRequested - req->mBytesPending,
			     &(req->mBytesPending), &errmsg) ;
	  if (state == ACECONNMSG_DONE)
	    {
	      req->done_message = TRUE ;
	    }
	  else if (state == ACECONNMSG_ERROR)
	    {
	      status = ACECONN_READERROR ;
	      aceconnSetErrMsg(connection, "Error reading message body from socket: ", errmsg) ;
	      free(errmsg) ;
	    }
	}
    }


  /* Reset the socket to blocking, try to do this even if we have an error.  */
  /* perhaps not worth it ?                                                  */
  if (status == ACECONN_OK)
    {
      status = resetBlockingSocket(connection, NULL, &sock_flags) ;
    }

  return status ;
}


/* Common read routine for handling reading from a socket and any resulting  */
/* errors. Note that it updates the the number of bytes left to read given   */
/* by bytesPending. If it fails it will give the reason in err_msg_out.      */
/*                                                                           */
/* n.b. should not be called with zero bytes to send, null readbuf etc.      */
/*                                                                           */
static AceConnMsgState readSocket(int fd, char *readbuf, int *bytes_pending_inout,
				  char **errmsg_out)
{
  AceConnMsgState state = ACECONNMSG_ERROR ;
  int bytes ;

  bytes = read(fd, readbuf, *bytes_pending_inout) ;
  if (bytes < 0)
    {
      /* This is potentially not strict enough, we may need to put a load of */
      /* #defines in here to check for SYS_V, POSIX etc. etc. to make sure   */
      /* we check for the correct errno's...aaaggghhh....                    */
      if (errno == EWOULDBLOCK || errno == EINTR || errno == EAGAIN)
	state = ACECONNMSG_WAIT ;
      else
	{
	  state = ACECONNMSG_ERROR ;
	  *errmsg_out = strdup("Read error trying to read client data on socket");
	}
    }
  else if (bytes == 0)
    {
      /* CHECK THIS ERROR CONDITION IN STEVENS....                           */
      state = ACECONNMSG_ERROR ;
      *errmsg_out = strdup("Zero bytes returned from read.") ;
    }
  else
    {
      *bytes_pending_inout -= bytes ;
      
      if (*bytes_pending_inout < 0)
	{
	  state = ACECONNMSG_ERROR ;
	  *errmsg_out = strdup("Logic error, read more bytes than number requested.") ;
	}
      else  if (*bytes_pending_inout == 0)
	state = ACECONNMSG_DONE ;
      else
	state = ACECONNMSG_WAIT ;
    }

  return state ;
}



/**************************************************************/

/* Jean I added code here to store any current handler for SIGPIPE, we don't */
/* know if the rest of the application will have installed one. We then      */
/* turn off signal handling for SIGPIPE so that we get EPIPE for a write to  */
/* a broken connection which we handle here. Then we turn back on any        */
/* existing signal handler.                                                  */
static AceConnStatus aceconnSocketWrite(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;
  AceConnMsgState state = ACECONNMSG_ERROR ;
  int sock_flags = -1 ;
  struct sigaction oursigpipe, oldsigpipe ;
  char *errmsg = NULL ;
  int fd = connection->socket ;
  AceConnMsg req = (AceConnMsg)connection->msg ;


/* bug: fix to work in my context */

  /* writes can deliver a SIGPIPE if the socket has been disconnected, by    */
  /* ignoring we will just receive errno = EPIPE.                            */
  oursigpipe.sa_handler = SIG_IGN ;
  sigemptyset(&oursigpipe.sa_mask) ;
  oursigpipe.sa_flags = 0 ;
  if (sigaction(SIGPIPE, &oursigpipe, &oldsigpipe) < 0)
    {
      status = ACECONN_SIGSET ;
      aceconnSetErrMsg(connection,
		       "Cannot set SIG_IGN for SIGPIPE for socket write operations: ", strerror(errno)) ;
    }

  /* We can't just do a recv with MSG_DONTWAIT and MSG_PEEK because          */
  /* the former is not supported on many systems yet. So instead we make the */
  /* whole socket non-blocking, do the recv and then reset the socket.       */
  if (status == ACECONN_OK)
    {
      status = resetBlockingSocket(connection, &sock_flags, NULL) ;
    }


  /* OK, Let's try to write the header                                       */

  if (status == ACECONN_OK)
    {
      if (req->done_header == FALSE)
	{
	  /* Do any byte swapping.                                               */
	  hdr2buf(&(req->proto_hdr), req->header_buf) ;

	  /* Now do the writes...                                                */
	  state = writeSocket(fd, req->header_buf + req->hBytesRequested - req->hBytesPending,
			      &(req->hBytesPending), &errmsg) ;
	  if (state == ACECONNMSG_DONE)
	    {
	      req->done_header = TRUE ;
	    }
	  else if (state == ACECONNMSG_ERROR)
	    {
	      status = ACECONN_WRITEERROR ;
	      aceconnSetErrMsg(connection, "Error reading message header from socket: ",
			       errmsg) ;
	      free(errmsg) ;
	    }
	}
    }

  /* OK, header written, so try to write buffer.                             */
  if (status == ACECONN_OK)
    {
      if (req->done_header == TRUE)
	{
	  state = writeSocket(fd,
			      req->message + (req->mBytesRequested - req->mBytesPending),
			      &(req->mBytesPending), &errmsg) ;
	  if (state == ACECONNMSG_DONE)
	    {
	      req->done_message = TRUE ;
	    }
	  else if (state == ACECONNMSG_ERROR)
	    {
	      status = ACECONN_WRITEERROR ;
	      aceconnSetErrMsg(connection, "Error reading message body from socket: ",
			       errmsg) ;
	      free(errmsg) ;
	    }
	}
    }

  if (status == ACECONN_OK)
    {
      status = resetBlockingSocket(connection, NULL, &sock_flags) ;
    }

  /* Reset the old signal handler, do this whatever our state.               */
  if (sigaction(SIGPIPE, &oldsigpipe, NULL) < 0)
    {
      status = ACECONN_SIGSET ;
      aceconnSetErrMsg(connection,
		       "Cannot reset previous handler for SIGPIPE after "
		       "socket write operations: ", strerror(errno)) ;
    }

  return status ;
}


/* Common write routine for handling writing to a socket and any resulting   */
/* errors. Note that it updates the the number of bytes left to write given  */
/* by bytesPending.                                                          */
static AceConnMsgState writeSocket(int fd, const char *writebuf, int *bytes_pending_inout,
				   char **errmsg_out)
{
  AceConnMsgState state = ACECONNMSG_ERROR ;
  int bytes ;

  bytes = write(fd, writebuf, *bytes_pending_inout) ;
  if (bytes < 0)
    {
      if (errno == EINTR || errno == EAGAIN)
	state = ACECONNMSG_WAIT ;
      else
	{
	  state = ACECONNMSG_ERROR ;
	  *errmsg_out = strdup("Write error trying to write client data to socket");
	}
    }
  else
    {
      *bytes_pending_inout -= bytes ;

      if (*bytes_pending_inout < 0)
	{
	  state = ACECONNMSG_ERROR ;
	  *errmsg_out = strdup("Logic error, wrote more bytes than number requested.") ;
	}
      else  if (*bytes_pending_inout == 0)
	state = ACECONNMSG_DONE ;
      else 
	state = ACECONNMSG_WAIT ;
    }

  return state ;
}


/*                                                                           */
/* Message handling code.                                                    */
/*                                                                           */

/* Allocate and initialise to zero.                                          */
/*                                                                           */
static AceConnMsg messageCreate(void)
{
  AceConnMsg msg = (AceConnMsg) calloc(sizeof(struct AceConnMsgRec),1) ;

  msg->proto_hdr.byte_swap = ACECONN_BYTE_ORDER ;		    /* Always the same. */

  return (AceConnMsg)msg ;
}


/* Set up a message to receive data from a socket.                           */
static void messageSetForRead(AceConnMsg msg_in, int encore)
{
  AceConnMsg msg = (AceConnMsg)msg_in ;

  /* Nothing to set in the protocol header.                                  */

  /* Set fields for controlling write of header and message.                 */
  msg->done_header = FALSE ;
  memset(msg->header_buf, 0, ACECONN_HEADER_BYTES) ;
  msg->hBytesRequested =  msg->hBytesPending = ACECONN_HEADER_BYTES ;

  msg->done_message = FALSE ;
  msg->message = NULL ;
  if (!encore && msg->encore_message.str)
    {
      free(msg->encore_message.str);
      msg->encore_message.str = NULL ;
    }
  msg->mBytesRequested =  msg->mBytesPending = 0 ;

  return ;
}


/* Set up a message to be written out to a socket.                           */
/*                                                                           */
static void messageSetForWrite(AceConnMsg msg_in, const char *message)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;

  /* Set fields in protocol header, may be special "encore" request meaning  */
  /* get next slice.                                                         */
  if (strcmp(message, ACECONN_SERVER_CLIENT_SLICE) == 0)
    setMsgType(msg->proto_hdr.msg_type, ACECONN_MSGENCORE) ;
  else
    setMsgType(msg->proto_hdr.msg_type, ACECONN_MSGREQ) ;
  msg->proto_hdr.length = strlen(message) + 1 ;

  /* Set fields for controlling write of header and message.                 */
  msg->done_header = FALSE ;
  memset(msg->header_buf, 0, ACECONN_HEADER_BYTES) ;
  msg->hBytesRequested = msg->hBytesPending = ACECONN_HEADER_BYTES ;

  msg->done_message = FALSE ;
  msg->message = (char *)message ;
  msg->mBytesPending = msg->mBytesRequested = msg->proto_hdr.length ;

  return ;
}

/* Has a message been completely sent ?                                      */
static int messageDone(AceConnMsg msg_in)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;
  int result = FALSE ;

  if (msg->done_header && msg->done_message)
    result = TRUE ;

  return result ;
}


/* Was message an "encore" message ?                                         */
static int messageIsEncore(AceConnMsg msg_in)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;
  int result = FALSE ;

  if (testMsgType(msg->proto_hdr.msg_type, ACECONN_MSGENCORE))
    result = TRUE ;

  return result ;
}

static void messageGetEncore(AceConnMsg msg_in)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;

  if (msg->message != NULL)
    {
      if (msg->encore_message.str == NULL)
	{
	msg->encore_message.str = (unsigned char *) malloc(msg->mBytesRequested);
	memcpy(msg->encore_message.str, msg->message, msg->mBytesRequested) ;
	}
      else
	{
	abort();
#if 0
 msg->encore_message = g_string_append(msg->encore_message, msg->message) ;
#endif
	}
      
      free((void *)msg->message) ;
      msg->message = NULL ;
    }

  return ;
}


static void messageGetReply(AceConnMsg msg_in, void **reply, int *reply_len)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;

  if (msg->encore_message.str)
    {
      *reply = msg->encore_message.str ;
      *reply_len = msg->encore_message.len ;
    }
  else
    {
      *reply = (void*) msg->message ;
      *reply_len = msg->mBytesRequested ;
    }

  return ;
}


/* Clean up a message of allocated resources                                 */
/* N.B. we DO NOT free any message buffer, that is always the clients        */
/* responsibility.                                                           */
/*                                                                           */
static void messageDestroy(AceConnMsg msg_in)
{
  AceConnMsg msg = (AceConnMsg) msg_in ;

  if (msg->encore_message.str)
    free(msg->encore_message.str) ;

  free(msg) ;

  return ;
}




/* Utilities to fiddle with protocol header.                                 */
/*                                                                           */
static void buf2hdr(char *buf, AceConnHeader hdr)
{
  memcpy(&(hdr->byte_swap), buf, 4) ;      buf += 4 ;
  memcpy(&(hdr->length), buf, 4) ;         buf += 4 ;
  memcpy(&(hdr->server_version), buf, 4) ; buf += 4 ;
  memcpy(&(hdr->client_id), buf, 4) ;      buf += 4 ;
  memcpy(&(hdr->max_bytes), buf, 4) ;      buf += 4 ;
  memcpy(&(hdr->msg_type), buf, ACECONN_MSGTYPE_BUFLEN) ;

  return ;
}

static void hdr2buf(AceConnHeader hdr, char *buf)
{
  memcpy(buf, &(hdr->byte_swap), 4) ;      buf += 4 ;
  memcpy(buf, &(hdr->length), 4) ;         buf += 4 ;
  memcpy(buf, &(hdr->server_version), 4) ; buf += 4 ;
  memcpy(buf, &(hdr->client_id), 4) ;      buf += 4 ;
  memcpy(buf, &(hdr->max_bytes), 4) ;      buf += 4 ;
  memcpy(buf, &(hdr->msg_type), ACECONN_MSGTYPE_BUFLEN) ;

  return ;
}


/* These two routines set the type and test it, they could be macros but     */
/* performance is not the problem here.                                      */
/*                                                                           */
/* They cope with caller supplying msgType which itself points to buffer.    */
/*                                                                           */
static void setMsgType(char buffer[], char *msgType)
{
  char *msg_copy = NULL ;

  /* Ugly bug here...what if msgType points to buffer ? We'll do belt and    */
  /* braces and clean the buffer anyway.                                     */
  if (msgType == &(buffer[0]))
    {
      msg_copy = strdup(msgType) ;
    }

  /* Reset the message type section of the header to be zeroed so there is   */
  /* no extraneous bumpf if the previous message was longer than this one.   */
  /* This is important for perl and other non-C languages that have to parse */
  /* this bit out of the message buffer.                                     */

  /* 2025_03_25: a false positive warning generated by gcc , avoid it */ 
/* 2025-compatible way: silence only GCC's false positive, ignored silently by Clang */
#ifndef __clang__
  #pragma GCC diagnostic push
  #pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

  memset(buffer, 0, ACECONN_MSGTYPE_BUFLEN) ;

#ifndef __clang__
  #pragma GCC diagnostic pop
#endif
  
  if (msg_copy != NULL)
    msgType = msg_copy ;

  strcpy(buffer, msgType) ;

  if (msg_copy != NULL)
    free(msg_copy) ;

  return ;
}


static int testMsgType(char buffer[], char *msgType)
{
  int result = FALSE ;

  if (msgType == &(buffer[0]))
    result = TRUE ;
  else if (strcmp(buffer, msgType) == 0)
    result = TRUE ;

  return result ;
}



/*                                                                           */
/*                    Socket utilities                                       */
/*                                                                           */

/* Currently we only wait for READ or WRITE, not both, could alter interface */
/* to allow both quite easily.                                               */
/* On success returns ACECONN_OK                                             */
/*                                                                           */
static AceConnStatus selectSocket(AceConnection connection, AceConnSockAccess access)
{
  AceConnStatus status = ACECONN_OK ;
  int rc ;
  struct timeval tv, *tv_ptr ;
  fd_set readset, writeset ;
  fd_set *readset_ptr = NULL, *writeset_ptr = NULL, *exceptset_ptr = NULL ;
  int maxfd ;

  if (access == ACECONNACC_WRITE)
    {
      FD_ZERO(&writeset) ;
      FD_SET(connection->socket, &writeset) ;
      writeset_ptr = &writeset ;
    }
  else
    {
      FD_ZERO(&readset) ;
      FD_SET(connection->socket, &readset) ;
      readset_ptr = &readset ;
    }

  maxfd = connection->socket + 1 ;			    /* "+ 1" is usual array offset gubbins */

  /* We allow the "never timeout" option, i.e. only return if the socket is  */
  /* ready.                                                                  */
  if (connection->timeout == 0)
    tv_ptr = NULL ;
  else
    {
      tv.tv_sec = connection->timeout ;
      tv.tv_usec = 0 ;
      tv_ptr = &tv ;
    }

  rc = select(maxfd, readset_ptr, writeset_ptr, exceptset_ptr, tv_ptr) ;
  if (rc < 0)
    { 
      status = ACECONN_NOSELECT ;
      aceconnSetErrMsg(connection, "select() on socket failed: ", strerror(errno)) ;
    }
  else if (rc == 0)
    {
      status = ACECONN_TIMEDOUT ;
      aceconnSetErrMsg(connection, "", "select() on socket timed out.") ;
    }
  else
    status = ACECONN_OK ;

  return status ;
}


/* Set/reset non-blocking/blocking on socket.                                */
/* If oldflags_out is non-NULL the current flags will be returned in it.     */
/* If newflags is non-NULL the socket will be set to those flags, otherwise  */
/* the existing socket flags will be OR'd with O_NONBLOCK.                   */
/*                                                                           */
static AceConnStatus resetBlockingSocket(AceConnection connection,
					 int *oldflags_out, int *newflags)
{
  AceConnStatus status = ACECONN_OK ;
  int fd = connection->socket ;
  int flags ;

  if (newflags)
    flags = *newflags ;
  else
    {
      flags = fcntl(fd, F_GETFL, 0) ;
      if (flags < 0)
	{
	  status = ACECONN_NONBLOCK ;
	  aceconnSetErrMsg(connection, "Failed to get socket flags: ",
			   strerror(errno)) ;
	}
      else
	{
	  if (oldflags_out)				    /* Return old flags ? */
	    *oldflags_out = flags ;
	  flags |= O_NONBLOCK ;
	}
    }

  if (status == ACECONN_OK)
    {
      flags = fcntl(fd, F_SETFL, flags) ;
      if (flags < 0)
	{
	  status = ACECONN_NONBLOCK ;
	  aceconnSetErrMsg(connection, "Failed to set socket flags: ",
			   strerror(errno)) ;
	}
    }

  return status ;
}


/* We ignore ENOTCONN which means that the socket has probably been closed   */
/* by the server already, other errors imply a coding error by us.           */
/*                                                                           */
static AceConnStatus shutdownSocket(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;

  if (shutdown(connection->socket, 2) < 0)
    {
      if (errno != ENOTCONN)
	{
	  status = ACECONN_INTERNAL ;
	  aceconnSetErrMsg(connection, 
			   "Internal code error, invalid option to shutdown call: ",
			   strerror(errno)) ;
	}
    }

   return status ;
}



/* Encryption and hashing code.                                              */
/*                                                                           */
/* We use MD5 for encryption and this can be found in the wmd5 directory.    */
/*                                                                           */
/* Note that MD5 produces as its output a 128 bit value that uniquely        */
/* represents the input strings. We convert this output value into a hex     */
/* string version of the 128 bit value for several reasons:                  */
/*                                                                           */
/*     - the md5 algorithm requires strings as input and we need to use      */
/*       some of the md5 output as input into a new md5 hash.                */
/*     - it makes all handling of the encrypted data much simpler            */
/*     - using hex means that the string we produce consists entirely of the */
/*       digits 0-9 and the letters a-f (i.e. no unprintable chars)          */
/*     - the passwd hash can be kept in a plain text file in the database    */
/*                                                                           */

/* hash the userid and password together and convert into a hex string.      */
/*                                                                           */
static char *makeHash (const char *userid, const char *passwd)
{
  char *hash = NULL ;
  enum {HASH_STRING_NUM = 2} ;
  const char *hash_strings[HASH_STRING_NUM] ;

  hash_strings[0] = userid ;
  hash_strings[1] = passwd ;
  hash = hashAndHexStrings (hash_strings, HASH_STRING_NUM) ;

  return hash ;
}

/* Takes an array of strings, hashes then together using MD5 and then        */
/* produces a hexstring translation of the MD5 output.                       */
/*                                                                           */
static char *hashAndHexStrings(const char *strings[], int num_strings)
{
  char *hex_string = NULL ;
  MD5_CTX Md5Ctx ;
  unsigned char digest[MD5_HASHLEN] ;
  int i ;

  MD5Init(&Md5Ctx) ;

  for (i = 0 ; i < num_strings ; i++)
    {
      MD5Update(&Md5Ctx, (unsigned char *)strings[i], strlen(strings[i])) ;
    }

  MD5Final(digest, &Md5Ctx) ;

  hex_string = convertMD5toHexStr(&digest[0]) ;

  return hex_string ;
}


/************************************************************/

/* Takes the array of unsigned char output by the MD5 routines and makes a   */
/* hexadecimal version of it as a null-terminated string.                    */
/*                                                                           */
static char *convertMD5toHexStr(unsigned char digest[])
{
  char *digest_str ;
  int i ;
  char *hex_ptr ;

  digest_str = (char *) malloc(MD5_HEX_HASHLEN) ;
  for (i = 0, hex_ptr = digest_str ; i < MD5_HASHLEN ; i++, hex_ptr+=2)
    {
      sprintf(hex_ptr, "%02x", digest[i]) ;
    }

  return digest_str ;
}

static AceConnStatus doServerClose(AceConnection connection) ;


/**@name Interface Functions
 * The AceConn package has a simple set of interface functions that allow
 * you to open a connection to an acedb server, send requests and receive
 * answers. */

/*@{ Start of interface functions. */


/** Create the an unconnected connection control block.
 *
 * @return ACECONN_OK if connect to server successful, otherwise
 * status indicates what failure was.
 *
 * @param connection_out returned valid connection (opaque type),
 *        untouched unless create succeeds.
 * @param server_netid host network id of server
 * @param server_port servers port on host
 * @param userid user known to acedb server
 * @param passwd users password (will be encrypted before sending).
 * @param timeout time to wait for reply from server before returning
 *        control to application
 */
static AceConnStatus AceConnCreate(AceConnection *connection_out,
			    const char *server_netid, int server_port,
			    const char *userid, const char *passwd, int timeout)
{
  AceConnStatus status = ACECONN_OK ;
  AceConnection connection = NULL ;

  /* Could check port number for  server_port < 1024 || server_port > 65535,
   * but probably over kill. */
  if (server_netid == NULL || server_port <= 0
      || userid == NULL || passwd == NULL)
    status = ACECONN_BADARGS ;
  else
    {
      if (!(connection = aceconnCreate(server_netid, server_port, userid, passwd, timeout)))
	status = ACECONN_NOCREATE ;
      else
	*connection_out = connection ;
    }

  return status ;
}

/** Open a connection to the server. Does all the hand shaking necessary to
 * say "hello" to the server, after this acedb requests can be sent to the
 * server.
 *
 * @return ACECONN_OK if connect to server successful, otherwise
 * status indicates what failure was.
 *
 * @param connection_inout server connection
 */
static AceConnStatus AceConnConnect(AceConnection connection_inout)
{
  AceConnStatus status = ACECONN_OK ;

  if (connection_inout == NULL) 
    {
      status = ACECONN_BADARGS ;
    }

  if (status == ACECONN_OK)
    {
      status = aceconnServerOpen(connection_inout) ;
    }

  if (status == ACECONN_OK)
    {
      status = aceconnServerHandshake(connection_inout) ;
    }

  /* If there was a problem, then clear up. */
  if (status != ACECONN_OK)
    {
      doServerClose(connection_inout) ;
    }

  return status ;
}



/** Sends an acedb request to the server and then waits for and returns
 * the reply from the server.
 *
 * @return ACECONN_OK if request sent and reply retrieved successfully, otherwise
 * status indicates what failure was.
 *
 * @param connection_out server connection
 * @param request, a C string containing the acedb request, e.g. "find sequence *"
 * @param reply_out, the reply received from the server, <B>note</B> this may be either
 * a C string or binary data as when a postscript image is returned.
 * @param reply_len_out the length of the data returned, <B>note</B> that for a C string,
 * this <b>includes</b> the terminating null char of the string.
 */
static AceConnStatus AceConnRequest(AceConnection connection,
			     const char *request,
			     void **reply_out, int *reply_len_out)
{
  AceConnStatus status = ACECONN_OK ;

  if (status == ACECONN_OK)
    {
      status = aceconnServerRequest(connection, request, reply_out, reply_len_out) ;

      /* For a normal request we just send it, but if the request was "quit" */
      /* then we close the connection to the server as well, note that if    */
      /* the close was successful we return a status of "quit" so the caller */
      /* can detect this situation.                                          */
      if (status == ACECONN_OK && (strcmp(request, ACECONN_CLIENT_GOODBYE) == 0))
	{
	  connection->state = ACECONNSTATE_OPEN ;
	  status = aceconnServerClose(connection) ;
	  if (status == ACECONN_OK)
	    {
	      status = ACECONN_QUIT ;
	      aceconnSetErrMsg(connection, "",
			       "Connection to server closed by \"quit\" request.") ;
	    }
	}
    }

  return status ;
}


/** Closes the connection to an acedb server. If close is
 * successful then connection is freed, otherwise it is not
 * to allow application to get the error message.
 *
 * @return ACECONN_OK if disconnected successfully, otherwise
 * status indicates what failure was.
 *
 * @param connection server connection (opaque type)
 */
static AceConnStatus AceConnDisconnect(AceConnection connection)
{
  AceConnStatus status;
  status = doServerClose(connection) ;
  return status ;
}


/** Simple free of a connection struct.
 *
 * @return nothing
 *
 * @param connection server connection (opaque type)
 */
static void AceConnDestroy(AceConnection connection)
{
  aceconnDestroy(connection) ;

  return ;
}

/*@} End of interface functions. */


/**@name Utility Functions
 * Following are all basically utility functions that supplement the Interface
 * functions. */


/*                                                                           */
/* Internal functions.                                                       */
/*                                                                           */

static AceConnStatus doServerClose(AceConnection connection)
{
  AceConnStatus status = ACECONN_OK ;

  if (status == ACECONN_OK && connection->state != ACECONNSTATE_NONE)
    {
      status = aceconnServerClose(connection) ;
    }

  return status ;
}


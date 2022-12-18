/* parser.c */

/*
	$Id: parser.c 8838 2003-11-30 07:04:15Z katayama $
	$Log$
	Revision 1.5  2003/11/30 07:04:14  katayama
	check return value of malloc

	Revision 1.4  2002/02/23 19:05:26  katayama
	Not sys/fcntl.h
	
	Revision 1.3  1998/07/28 23:58:42  katayama
	Less warnings with g++ -Wall
	
	Revision 1.2  1998/06/18 07:47:11  higuchit
	*** empty log message ***

	Revision 1.1  1998/06/17 08:42:59  higuchit
	*** empty log message ***

*/

#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <malloc.h>
#include <string.h>
#include <unistd.h>

#include "parser.h"

static void prsnull(buf)
char *buf;
{
	char *ptr;

	for(ptr=buf;*ptr;ptr++){
		switch(*ptr){
			case '\t':
				*ptr = ' ';
				break;

			case '\n':
				if(ptr!=buf && *(ptr-1) == '\\')
					*ptr = '\2', *(ptr-1) = ' ';
				else
					*ptr = '\1';
				break;
			
			/*case '#' :*/
			case '!':
			case ';':
				for(;*ptr;ptr++){
					if(*ptr == '\n')
						break;
					else
					if(*ptr == '\\' && *(ptr+1) == '\n')
						*ptr = '\2', *(ptr+1) = ' ';
					else
						*ptr = ' ';
				}
				*ptr = '\1';
		}
	}
}

PRS *prsopen(file)
char *file;
{
	PRS *prs;
	int fd;
	struct stat stbuf;

	prs = (PRS*)malloc(sizeof(PRS));
	if(!prs) {
       	  perror("$Source:prs:malloc");
	  return NULL;
	}
	
	fd = -1;
	prs->master = prs->buf = NULL;

	if((fd=open(file,O_RDONLY))==-1) goto ERROR;
	if(fstat(fd,&stbuf)==-1) goto ERROR;
	if(!(prs->master=(char*)malloc(stbuf.st_size+1))) goto ERROR;
	if(!(prs->buf   =(char*)malloc(stbuf.st_size+1))) goto ERROR;
	if(read(fd,prs->master,stbuf.st_size)!=stbuf.st_size) goto ERROR;
	*(prs->master+stbuf.st_size) = '\0';
	close(fd);

	prsrewind(prs);

	return prs;

	ERROR:
	close(fd);
	free(prs->master);
	free(prs->buf);
	free(prs);
	return NULL;
}

void prsrewind(prs)
PRS *prs;
{
	prs->pline = 1;
	prs->cline = 1;
	strcpy(prs->buf,prs->master);
	prsnull(prs->buf);
	prs->ptr = prs->buf;
}

int prstell(prs)
PRS *prs;
{
	return prs->pline;
}

int prsnext(prs,argc,argv)
PRS *prs;
int *argc;
char ***argv;
{
	static char *v[64];
	char *new = NULL;
	int has_ent = 0;

	*argc = 0;

	prs->pline = prs->cline;
	new = prs->ptr;

	for(;;prs->ptr++)
		switch(*(prs->ptr)){
			case '\0':
				return 0;

			case '\1':
				*(prs->ptr) = ' ';
				prs->cline++;
				if(has_ent){ *(prs->ptr++) = '\0'; goto NL_END;}
				prs->pline = prs->cline;
				continue;

			case '\2':
				*(prs->ptr) = ' ';
				prs->cline++;
				continue;

			case ' ' :
				continue;

			default  :
				has_ent = 1;
		}
	NL_END:

	if(!has_ent) return 0;

	v[(*argc)++] = strtok(new," ");
	do{
		char *tmp = v[*argc-1];
		while(*tmp && *tmp!=' ') tmp++;
		*tmp = '\0';
	}while((v[(*argc)++] = strtok(NULL," ")));

	(*argc)--;
	*argv = v;

	return 1;
}

void prsclose(prs)
PRS *prs;
{
	free(prs->master);
	free(prs->buf);
	free(prs);
}


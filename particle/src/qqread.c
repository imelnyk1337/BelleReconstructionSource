/* qqread.c */

/*
	$Id: qqread.c 8838 2003-11-30 07:04:15Z katayama $
	$Log$
	Revision 1.8  2003/11/30 07:04:15  katayama
	check return value of malloc

	Revision 1.7  2001/12/05 12:35:12  katayama
	For gcc-3
	
	Revision 1.6  1998/10/14 12:50:24  jtanaka
	uses qq2panther to get Ptype information.
	
	Revision 1.5  1998/10/05 10:55:56  jtanaka
	bug fix in qqread.c and modified static objects in Relation.cc

	Revision 1.4  1998/09/22 06:14:46  jtanaka
	bug fix : bug = no user.dec -> can not use Ptype(int PDGID).

	Revision 1.3  1998/09/09 08:07:00  jtanaka
	added Ptype(idhep).

	Revision 1.2  1998/07/28 23:58:43  katayama
	Less warnings with g++ -Wall

	Revision 1.1  1998/07/16 07:03:06  higuchit
	QQ2PDT -> QQREAD
	
*/

/* killed on 10/8/98 */
#if 0

/**********************/
/****** INCLUDES ******/
/**********************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string.h>
#include <values.h>
#include <sys/types.h>
#include "qqread.h"
#include "parser.h"

/*******************************************/
/****** PARTICLE PROPERTIES FUNCTIONS ******/
/*******************************************/

static size_t nPDT;
static struct pdtent  PDTarray[MAXNPDT];
static struct pdtent *PDTarray2[MAXNPDT]; /* copy pointers of PDTarray */

static struct pdtent *pdtent_build(argc,argv)
int argc;
char **argv;
{
	static struct pdtent pd;

	strncpy(pd.name,argv[1],sizeof(pd.name)-1);
	pd.pdgid  = 0;
	pd.stable = atoi(argv[3]) != -1;
	pd.mass   = atof(argv[4]);
	pd.charge = atof(argv[5]);
	pd.spin   = atof(argv[6]);
	pd.ctau   = atof(argv[7]);

	pd.width    = (argc==11) ? atof(argv[ 8]) : -1.0;
	pd.mass_min = (argc==11) ? atof(argv[ 9]) :  0.0;
	pd.mass_max = (argc==11) ? atof(argv[10]) :  0.0;

	return &pd;
}

static int pdtent_append(pd)
struct pdtent *pd;
{
	if(nPDT>=MAXNPDT) return 0;

	bcopy(pd,&(PDTarray[nPDT++]),sizeof(struct pdtent));

	return 1;
}

static int pdtent_cmp_with_name(pd1,pd2)
struct pdtent *pd1, *pd2;
{
	/* comparison of names only */
	return strcmp(pd1->name,pd2->name);
}

static int pdtent_cmp_with_pdgid(pd1,pd2)
struct pdtent **pd1, **pd2;
{
	/* comparison of pdgid only */
	return (*pd1)->pdgid - (*pd2)->pdgid;
}

static void pdtent_sort_with_name()
{
	qsort(PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp_with_name);
}

static void pdtent_sort_with_pdgid()
{
	qsort(PDTarray2,nPDT,sizeof(struct pdtent*),pdtent_cmp_with_pdgid);
}

static void pdtent_do(prs)
PRS *prs;
{
	int argc;
	char **argv;

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){
		struct pdtent *pd;

		if((argc!=8 && argc!=11) || strcasecmp("PARTICLE",argv[0]))
			continue;

		pd = pdtent_build(argc,argv);
		pdtent_append(pd);
	}

	pdtent_sort_with_name();

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){
		struct pdtent *pd, pdtmp;

		if(argc!=3 || strcasecmp("PDG",argv[0]))
			continue;

		strncpy(pdtmp.name,argv[1],sizeof(pdtmp.name)-1);
		pd = (struct pdtent*)bsearch(
			&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp_with_name
		);

		if(pd) pd->pdgid = atoi(argv[2]);
	}
}


/***********************************/
/****** DECAY TABLE FUNCTIONS ******/
/***********************************/

static size_t nDECAY;
static struct decay DECAYarray[MAXNDECAY];

#define DECAY_CMP_MOTHER_ONLY MAXINT

static int pname_cmp(pname1,pname2)
pname *pname1, *pname2;
{
	return strcmp((char*)pname1,(char*)pname2);
}

static struct decay *decay_build(argc,argv)
int argc;
char **argv;
{
	int i;
	static struct decay dec;
	pname daughter_tmp[8];

	dec.matrix    = atoi(argv[1]);
	dec.br        = atof(argv[2]);
	dec.ndaughter = argc-3;

	for(i=3;i<argc;i++) strcpy(daughter_tmp[i-3],argv[i]);
	qsort(daughter_tmp,dec.ndaughter,sizeof(pname),pname_cmp);
	for(i=0;i<dec.ndaughter;i++) strcpy(dec.daughter[i],daughter_tmp[i]);

	return &dec;
}

static int decay_append(dec)
struct decay *dec;
{
	if(nDECAY>=MAXNDECAY) return 0;

	bcopy(dec,&(DECAYarray[nDECAY++]),sizeof(struct decay));

	return 1;
}

static int decay_cmp(dec1,dec2)
struct decay *dec1, *dec2;
{
	int cmp, ncheck, i;

	cmp = strcmp(dec1->mother,dec2->mother);
	if(
		cmp                                        ||
		(dec1->ndaughter == DECAY_CMP_MOTHER_ONLY) ||
		(dec2->ndaughter == DECAY_CMP_MOTHER_ONLY)
	)
		return cmp;

	cmp = dec1->ndaughter - dec2->ndaughter;
	if(cmp) return cmp;

	ncheck = dec1->ndaughter;
	for(i=0;i<ncheck;i++){
		cmp = strcmp(dec1->daughter[i],dec2->daughter[i]);
		if(cmp) return cmp;
	}

	return 0;
}

static void decay_sort()
{
	qsort(DECAYarray,nDECAY,sizeof(struct decay),decay_cmp);
}

static void decay_normalize()
{
	int i=0;

	while(i<nDECAY){
		pname *mother;

		mother = &(DECAYarray[i].mother);

		/* stability check */
		{
			struct pdtent *pd, pdtmp;

			strncpy(pdtmp.name,(char*)mother,sizeof(pdtmp.name)-1);
			pd = (struct pdtent*)bsearch(
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp_with_name
			);

			if(!pd || pd->stable){
				/* partcile was not found or is stable */
				while( i<nDECAY &&
					!strcmp((char*)mother,DECAYarray[i].mother) ) i++;
				continue;
			}
		}

		{
			int istart=i, iend = 0, j;
			float totbr=0.0;

			for(;i<nDECAY;i++){
				if(strcmp((char*)mother,DECAYarray[i].mother)){
					iend = i-1;
					break;
				}
				totbr += DECAYarray[i].br;
			}

			for(j=istart;j<=iend;j++)
				DECAYarray[j].br = totbr>0.0 ? DECAYarray[j].br/totbr : 0.0;
		}
	}
}

static void decay_delete(prs)
PRS *prs;
{
	int argc;
	char **argv;
	int ndelete = 0;
	struct decay *decstart = &(DECAYarray[     0]),
							 *decend   = &(DECAYarray[nDECAY]);

	prsrewind(prs);

	while(prsnext(prs,&argc,&argv)){
		struct decay *dec, *dec_i, dectmp;

		/* keyword DECAY */
		if(argc!=2 || strcasecmp("DECAY",argv[0])) continue;

		strcpy((char*)dectmp.mother,argv[1]);
		dectmp.ndaughter = DECAY_CMP_MOTHER_ONLY;

		dec = (struct decay*)bsearch(
			&dectmp,DECAYarray,nDECAY,sizeof(struct decay),decay_cmp
		);

		if(!dec) continue;

		for(dec_i=dec  ;dec_i>=decstart;dec_i--){
			if(strcmp((char*)dectmp.mother,(char*)dec_i->mother)) break;
			((char*)dec_i->mother)[0] = '\177'; /* sorted to be last */
			((char*)dec_i->mother)[1] = '\0';
			ndelete++;
		}
		for(dec_i=dec+1;dec_i<=decend  ;dec_i++){
			if(strcmp((char*)dectmp.mother,(char*)dec_i->mother)) break;
			((char*)dec_i->mother)[0] = '\177'; /* sorted to be last */
			((char*)dec_i->mother)[1] = '\0';
			ndelete++;
		}
	}

	decay_sort();

	/* truncate "\177" particles */
	nDECAY -= ndelete;
}

static void decay_do(prs)
PRS *prs;
{
	int argc;
	char **argv;
	pname mother;
	int decay_block = 0;

	/* preparation for overriding of new decay table */
	decay_delete(prs);

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){

		/* keyword STABLE */
		if(argc==2 && !strcasecmp("STABLE",argv[0])){
			/* NOTE: override pdtent.stable */

			struct pdtent *pd, pdtmp;

			strncpy(pdtmp.name,argv[1],sizeof(pdtmp.name)-1);
			pd = (struct pdtent*)bsearch(
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp_with_name
			);

			if(pd) pd->stable = 1;
			continue;
		}

		/* keyword DECAY */
		if(argc==2 && !strcasecmp("DECAY",argv[0])){
			strcpy((char*)mother,argv[1]);
			decay_block = 1;
			continue;
		}

		/* keyword ENDDECAY */
		if(argc==1 && !strcasecmp("ENDDECAY",argv[0])){
			decay_block = 0;
			continue;
		}

		/* out of DECAY - ENDDECAY block */
		if(!decay_block) continue;

		/* keyword CHANNEL ? */
		if(argc<3 || strcasecmp("CHANNEL",argv[0])) continue;

		{
			struct decay *dec;
			dec = decay_build(argc,argv);
			strcpy((char*)dec->mother,(char*)mother);
			decay_append(dec);
		}
	}

	decay_sort();
	decay_normalize();
}

/*****************************/
/****** QQ TABLE READER ******/
/*****************************/

int qqread(file_qq,nfile_qq)
char *file_qq[];
size_t nfile_qq;
{
	int i;

	/* parser */
	for(i=0;i<nfile_qq;i++){
		PRS *prs;
		
		/* if(!(prs=prsopen(file_qq[i]))) return 0; */
		if(i == 0 && !(prs=prsopen(file_qq[i]))) return 0;
		if(!(prs=prsopen(file_qq[i]))) goto for_pdgid;
	
		pdtent_do(prs);
		decay_do (prs);

		prsclose(prs);
	}

for_pdgid:
	for(i=0;i<nPDT;i++){
	  PDTarray2[i] = &PDTarray[i];
	}

	pdtent_sort_with_pdgid();

	return 1;
}


struct pdtent *qqgetpdt_with_name(pname name)
{
      struct pdtent pdtmp;

	strcpy(pdtmp.name,name);
	
	return (struct pdtent*)bsearch(
		&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp_with_name
	);
}


struct pdtent **qqgetpdt_with_pdgid(int pdgid)
{
  struct pdtent *pdtmp;
  struct pdtent **tmp;
  
  pdtmp = (struct pdtent*)malloc(sizeof(struct pdtent));
  if(pdtmp==NULL) return NULL;
  
  pdtmp->pdgid = pdgid;
  
  tmp = (struct pdtent**)bsearch(&pdtmp, PDTarray2, nPDT,
					   sizeof(struct pdtent*), pdtent_cmp_with_pdgid);

  free(pdtmp);
  
  return tmp;
}

#endif

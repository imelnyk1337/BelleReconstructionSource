/* qq2panther.c */

/**********************/
/****** INCLUDES ******/
/**********************/
#include "belle.h"

#include <stdio.h>
#include <stdlib.h>
/* for non-ansi string functions */
/* strcasecmp */
#if defined(__sparc)
#  if defined(__EXTENSIONS__)
#    include <string.h>
#  else
#    define __EXTENSIONS__
#    include <string.h>
#    undef __EXTENSIONS__
#  endif
#elif defined(__GNUC__)
#  if defined(_XOPEN_SOURCE)
#    include <string.h>
#  else
#    define _XOPEN_SOURCE
#    include <string.h>
#    undef _XOPEN_SOURCE
#  endif
#endif
#include <limits.h>
#include <stdarg.h>
#include <sys/types.h>
#include "qq2panther.h"
#include "parser.h"
#include QQ_H

typedef char pname[8];

#define MAXNPDT 1024
#define MAXNDECAY 4096

struct pdtent{
	pname name;
	int pdgid;
	int stable;
	float mass /* in GeV */, charge, spin, ctau;
	float width /* FWHM in GeV */, mass_min, mass_max;
};

struct decay{
	pname mother;
	int matrix;
	float br;
	size_t ndaughter;
	pname daughter[8];
};

int qqread();
struct pdtent *qqgetpdt();

/*******************************************/
/****** PARTICLE PROPERTIES FUNCTIONS ******/
/*******************************************/

static size_t nPDT;
static struct pdtent PDTarray[MAXNPDT];

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

	memcpy(&(PDTarray[nPDT++]),pd,sizeof(struct pdtent));

	return 1;
}

static int pdtent_cmp(pd1,pd2)
struct pdtent *pd1, *pd2;
{
	/* comparison of names only */
	return strcmp(pd1->name,pd2->name);
}

typedef int (*sortfunc)(const void*, const void*);

static void pdtent_sort()
{
	qsort(PDTarray,nPDT,sizeof(struct pdtent),(sortfunc)pdtent_cmp);
}

static void pdtent_do(prs)
PRS *prs;
{
	int argc;
	char **argv;

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){
		struct pdtent *pd;

		if(argc==3 && !strcasecmp("QQBAR",argv[0])){

			static struct pdtent _pd;
			pd = &_pd;

			strncpy(pd->name,argv[1],sizeof(pd->name)-1);
			pd->pdgid  = 0;
			pd->stable = 0;
			pd->mass   = 0.0;
			pd->charge = 0.0;
			pd->spin   = 0.0;
			pd->ctau   = 0.0;

			pd->width    = 0.0;
			pd->mass_min = 0.0;
			pd->mass_max = 0.0;

		}else
		if((argc==8 || argc==11) && !strcasecmp("PARTICLE",argv[0])){

			pd = pdtent_build(argc,argv);

		}else continue;

		pdtent_append(pd);
	}

	pdtent_sort();

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){
	  if(strcasecmp("PDG",argv[0]) == 0){
	    struct pdtent *pd, pdtmp;
	    
	    /* if(argc!=3 || strcasecmp("PDG",argv[0]))continue; */
	    if(argc >= 3){
	      /* 3 is standard */
	      /* but there is exception!! */
	      strncpy(pdtmp.name,argv[1],sizeof(pdtmp.name)-1);
	      pd = (struct pdtent*)bsearch(
					   &pdtmp,PDTarray,nPDT,sizeof(struct pdtent),(sortfunc)pdtent_cmp
					   );
	      
	      if(pd) pd->pdgid = atoi(argv[2]);
	    }
	  }
	}
}

/***********************************/
/****** DECAY TABLE FUNCTIONS ******/
/***********************************/

static size_t nDECAY;
static struct decay DECAYarray[MAXNDECAY];

#define DECAY_CMP_MOTHER_ONLY INT_MAX

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
	qsort(daughter_tmp,dec.ndaughter,sizeof(pname),(sortfunc)pname_cmp);
	for(i=0;i<(int)dec.ndaughter;i++) strcpy(dec.daughter[i],daughter_tmp[i]);

	return &dec;
}

static int decay_append(dec)
struct decay *dec;
{
	if(nDECAY>=MAXNDECAY) return 0;

	memcpy(&(DECAYarray[nDECAY++]),dec,sizeof(struct decay));

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
	qsort(DECAYarray,nDECAY,sizeof(struct decay),(sortfunc)decay_cmp);
}

static void decay_normalize()
{
	size_t i=0;

	while(i<nDECAY){
		pname *mother;

		mother = &(DECAYarray[i].mother);

		/* stability check */
		{
			struct pdtent *pd, pdtmp;

			strncpy(pdtmp.name,(char*)mother,sizeof(pdtmp.name)-1);
			pd = (struct pdtent*)bsearch(
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),(sortfunc)pdtent_cmp
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
			&dectmp,DECAYarray,nDECAY,sizeof(struct decay),(sortfunc)decay_cmp
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
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),(sortfunc)pdtent_cmp
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

/****************************/
/****** PANTHER FILLER ******/
/****************************/

static void qqpanther()
{
	size_t i;
	pname moname;

	/* QQ_particle */
	for(i=0;i<nPDT;i++){
		struct qq_particle *qqpdt;

		qqpdt = (struct qq_particle*)BsNewEnt(QQ_PARTICLE);

		strcpy(qqpdt->m_name,PDTarray[i].name);
		qqpdt->m_pdgid = PDTarray[i].pdgid;
		qqpdt->m_charge = PDTarray[i].charge;
		qqpdt->m_mass = PDTarray[i].mass;
		qqpdt->m_ctau = PDTarray[i].ctau;
		qqpdt->m_spin = PDTarray[i].spin;
		qqpdt->m_decay[0] = qqpdt->m_decay[1] = 0;
	}

	/* QQ_decay */
	i = 0;
	while(i<nDECAY){
		struct pdtent *mopdt;
		int moqqpdtid;
		struct qq_particle *moqqpdt;

		strcpy(moname,DECAYarray[i].mother);

		mopdt = qqgetpdt(moname);
		moqqpdtid = mopdt - PDTarray + 1 /* +1 is for PANTHER */;

		moqqpdt =
			(struct qq_particle*)BsGetEnt(QQ_PARTICLE,moqqpdtid,BBS_No_Index);
		moqqpdt->m_decay[0] = i + 1;

		for(;i<nDECAY;i++){
			int j;
			struct qq_decay *qqdec;
			
			if(strcmp(moname,DECAYarray[i].mother)) break;

			qqdec = (struct qq_decay*)BsNewEnt(QQ_DECAY);
			qqdec->m_matrix = DECAYarray[i].matrix;
			qqdec->m_ratio = DECAYarray[i].br;
			qqdec->m_mother = moqqpdtid;
			qqdec->m_numdg = DECAYarray[i].ndaughter;
			for(j=0;j<qqdec->m_numdg;j++){
				int qqpdtid;
				struct pdtent *pdt;
				pdt = qqgetpdt(DECAYarray[i].daughter[j]);
				qqpdtid = pdt - PDTarray + 1;
				qqdec->m_dg[j] = qqpdtid;
			}
		}
		moqqpdt->m_decay[1] = i;
	}
}

/*****************************/
/****** QQ TABLE READER ******/
/*****************************/

int qq2panther_(const char *file_qq, ...)
{
	va_list var;
	char *file_user;
	PRS *prs;

	BsNoClrT(QQ_PARTICLE);
	BsNoClrT(QQ_DECAY);

	va_start(var,file_qq);

	/* decay.dec */
	if(!(prs=prsopen(file_qq))) return 0;

	pdtent_do(prs);
	decay_do (prs);

	prsclose(prs);

	/* user.dec */
	while((file_user = va_arg(var,char*))){
	        if(file_user==0 || *file_user==0) break;
		if(!(prs=prsopen(file_user))) continue;
	
		pdtent_do(prs);
		decay_do (prs);

		prsclose(prs);
	}

	va_end(var);

	qqpanther();

	return 1;
}

struct pdtent *qqgetpdt(name)
pname name;
{
	struct pdtent pdtmp;

	strcpy(pdtmp.name,name);
	
	return (struct pdtent*)bsearch(
		&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),(sortfunc)pdtent_cmp
	);
}


/* qq2pdt.c */

/*
	$Id: qq2pdt.c 8838 2003-11-30 07:04:15Z katayama $
	$Log$
	Revision 1.4  2003/11/30 07:04:14  katayama
	check return value of malloc

	Revision 1.3  1998/07/16 07:03:34  higuchit
	obsolete
	
	Revision 1.2  1998/06/18 07:47:21  higuchit
	*** empty log message ***

	Revision 1.1  1998/06/17 08:43:30  higuchit
	*** empty log message ***
*/

/* obsoleted on 7/16/98 */
#if 0

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include "parser.h"

/* type of particle name */
typedef char pname[16];

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

#define MAXNPDT 1024
static size_t nPDT;
static struct pdtent *PDTarray;

#define MAXNDECAY 4096
static size_t nDECAY;
static struct decay *DECAYarray;

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

	bcopy(pd,PDTarray+(nPDT++),sizeof(struct pdtent));

	return 1;
}

static int pdtent_cmp(pd1,pd2)
struct pdtent *pd1, *pd2;
{
	/* comparison of names only */
	return strcmp(pd1->name,pd2->name);
}

static void pdtent_sort()
{
	qsort(PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp);
}

static void qq2pdt_pdtent(prs)
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

	pdtent_sort();

	prsrewind(prs);
	while(prsnext(prs,&argc,&argv)){
		struct pdtent *pd, pdtmp;

		if(argc!=3 || strcasecmp("PDG",argv[0]))
			continue;

		strncpy(pdtmp.name,argv[1],sizeof(pdtmp.name)-1);
		pd = (struct pdtent*)bsearch(
			&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp
		);

		if(pd) pd->pdgid = atoi(argv[2]);
	}
}

static struct decay *decay_build(argc,argv)
int argc;
char **argv;
{
	int i;
	static struct decay dec;

	dec.matrix    = atoi(argv[1]);
	dec.br        = atof(argv[2]);
	dec.ndaughter = argc-3;
	for(i=3;i<argc;i++) strcpy(dec.daughter[i-3],argv[i]);

	return &dec;
}

static int decay_append(dec)
struct decay *dec;
{
	if(nDECAY>=MAXNDECAY) return 0;

	bcopy(dec,DECAYarray+(nDECAY++),sizeof(struct decay));

	return 1;
}

static int decay_cmp(dec1,dec2)
struct decay *dec1, *dec2;
{
	/* comparison of mother names only */
	return strcmp(dec1->mother,dec2->mother);
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

		mother = &((DECAYarray+i)->mother);

		/* stability check */
		{
			struct pdtent *pd, pdtmp;

			strncpy(pdtmp.name,(char*)mother,sizeof(pdtmp.name)-1);
			pd = (struct pdtent*)bsearch(
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp
			);

			if(!pd || pd->stable){
				/* partcile was not found or is stable */
				while( i<nDECAY &&
					!strcmp((char*)mother,(DECAYarray+i)->mother) ) i++;
				continue;
			}
		}

		{
			int istart=i, iend, j;
			float totbr=0.0;

			for(;i<nDECAY;i++){
				if(strcmp((char*)mother,(DECAYarray+i)->mother)){
					iend = i-1;
					break;
				}
				totbr += (DECAYarray+i)->br;
			}

			for(j=istart;j<=iend;j++)
				(DECAYarray+j)->br = totbr>0.0 ? (DECAYarray+j)->br/totbr : 0.0;
		}
	}
}

static void decay_delete(prs)
PRS *prs;
{
	int argc;
	char **argv;
	int ndelete = 0;
	struct decay *decstart = DECAYarray, *decend = DECAYarray + nDECAY;

	prsrewind(prs);

	while(prsnext(prs,&argc,&argv)){
		struct decay *dec, *dec_i, dectmp;

		/* keyword DECAY */
		if(argc!=2 || strcasecmp("DECAY",argv[0])) continue;

		strcpy((char*)dectmp.mother,argv[1]);

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

	nDECAY -= ndelete;
}

static void qq2pdt_decay(prs)
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
				&pdtmp,PDTarray,nPDT,sizeof(struct pdtent),pdtent_cmp
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

int qq2pdt(file_qq,nfile_qq,file_pdt)
size_t nfile_qq;
char *file_qq[], *file_pdt;
{
	int i;

	PDTarray   = malloc(sizeof(struct pdtent)*MAXNPDT  );
	if(PDTarray==NULL) {
	  perror("$Source$:PDTarray:malloc");
	  return 0;
	}
	DECAYarray = malloc(sizeof(struct decay )*MAXNDECAY);
	if (DECAYarray==NULL) {
       	  perror("$Source:DECAYarray:malloc");
	  return 0;
	}
	bzero(PDTarray  ,sizeof(struct pdtent)*MAXNPDT  );
	bzero(DECAYarray,sizeof(struct decay )*MAXNDECAY);

	/* parser */
	for(i=0;i<nfile_qq;i++){
		PRS *prs;
		
		if(!(prs=prsopen(file_qq[i]))) return 0;
	
		qq2pdt_pdtent(prs);
		qq2pdt_decay (prs);

		prsclose(prs);
	}

	/* dump */
	{
		FILE *fp_pdt;
		int i;

		if(!(fp_pdt=fopen(file_pdt,"w"))) return 0;

		for(i=0;i<nPDT;i++){
			struct pdtent *pd = PDTarray+i;
			fprintf(
				fp_pdt,
				"P %d %s %c %f %f %f %f 0.0 0.0 %f 0.0 0.0 %f %f 0.0 0 JETSET Default 0\n",
				pd->pdgid, pd->name, pd->stable?'y':'n',
				pd->mass, pd->mass_min, pd->mass_max,
				pd->width, pd->ctau, pd->charge, pd->spin
			);
		}

#if 0
		for(i=0;i<nDECAY;i++){
			int j;
			char decay_to[128];
			struct decay *dec = DECAYarray+i;

			fprintf(fp_pdt,"D %s:",dec->mother);

			decay_to[0] = '\0';
			for(j=0;j<dec->ndaughter;j++){
				strcat(decay_to,",");
				strcat(decay_to,dec->daughter[j]);
			}
			fprintf(fp_pdt,"%s; ",&(decay_to[1]));

			fprintf(fp_pdt,"%lf 0 0 1 JETSET %lf LUDECY000\n",dec->br,dec->br);
		}
#endif

		fclose(fp_pdt);
	}

	free(PDTarray  );
	free(DECAYarray);

	return 1;
}
#endif


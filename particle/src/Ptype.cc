#include <string.h>
// Ptype.cc
//
// $Log$
// Revision 1.28  2003/12/06 10:56:34  katayama
// Read from share/bin
//
// Revision 1.27  2003/11/30 07:04:14  katayama
// check return value of malloc
//
// Revision 1.26  2003/09/14 00:05:28  katayama
// for evtgen
//
// Revision 1.25  2001/12/24 12:03:32  katayama
// gcc 3.0 and headers are cleaned up
//
// Revision 1.24  2001/12/05 12:35:12  katayama
// For gcc-3
//
// Revision 1.23  2001/04/09 04:54:34  katayama
// Do not overwrite when a new qq_particle table is read from panther
//
// Revision 1.22  2000/04/13 12:41:57  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.21  2000/03/09 06:44:54  jtanaka
// bug fix.
// The bug is due to "additional comments in decay.dec"
//
// Revision 1.20  2000/03/07 11:14:04  katayama
// compatibility with CC5.0
//
// Revision 1.19  2000/01/05 06:09:39  jtanaka
// major updates: please see BELLE whiteboard.
//
// Revision 1.18  1999/07/20 13:57:21  jtanaka
// bug: user.dec
//
// Revision 1.17  1999/07/20 10:57:19  jtanaka
// default : no "user.dec"
//
// Revision 1.16  1999/07/20 10:48:04  jtanaka
// added some methods to use decay.dec and user.dec
//
// Revision 1.15  1999/06/07 07:10:45  higuchit
// higuchit
// getenv("BELLE_TOP_DIR")
//
// Revision 1.14  1999/04/14 10:32:34  higuchit
// getenv("QQ_USER_TABLE")
//
// Revision 1.13  1999/01/16 10:31:41  katayama
// clean up includes
//
// Revision 1.12  1999/01/09 12:51:36  jtanaka
// Bugs fix, added some members in Ptype Class, removed some members from Particle.cc because of functions for setting member of other classes.
//
// Revision 1.11  1998/12/03 11:30:53  katayama
// added belle.h
//
// Revision 1.10  1998/10/14 12:50:24  jtanaka
// uses qq2panther to get Ptype information.
//
// Revision 1.9  1998/09/09 08:07:00  jtanaka
// added Ptype(idhep).
//
// Revision 1.8  1998/09/08 09:54:01  jtanaka
// Modify Mdst_charged constructor in Particle and some functions in Momentum.
//
// Revision 1.7  1998/09/03 23:55:15  jtanaka
// We need to include "Particle.h" only. Comment out some members of Ptype.
//
// Revision 1.6  1998/08/03 06:06:43  katayama
// Read decay.dec/udecay.dec used in qq (temp. fix for now)
//
// Revision 1.5  1998/07/16 07:05:00  higuchit
// gave up of using PDT/CLHEP
//
// Revision 1.4  1998/07/03 12:08:30  jtanaka
// modified some parts
//
// Revision 1.3  1998/07/02 09:28:23  higuchit
// null flag -> `usable' flag
//
// Revision 1.2  1998/07/01 11:54:06  jtanaka
// add m_null.
//
// Revision 1.1  1998/06/14 14:11:32  higuchit
// *** empty log message ***
//
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "belle.h"

#include "particle/Ptype.h"
#include "panther/panther.h"
#include QQ_H

#include "qq2panther.h"
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


extern "C" {
	void qqgfil_(const char *, char *, int *, int, int);
	int qq2panther_(const char *, ...);
}


static bool use_EVTGEN;

static int nQQParticle;
static struct qq_particle **SortName, **SortPDGid;
static struct qq_particle *m_QQ_PARTICLE_LIST;


static int
cmp_name(const void *a, const void *b)
{
  const struct qq_particle **aa((const struct qq_particle **)a);
	const struct qq_particle **bb((const struct qq_particle **)b);
	return strcmp((*aa)->m_name,(*bb)->m_name);
}

static int
cmp_pdgid(const void *a, const void *b)
{
  const struct qq_particle **aa((const struct qq_particle **)a);
	const struct qq_particle **bb((const struct qq_particle **)b);
	return (*aa)->m_pdgid - (*bb)->m_pdgid;
}


static char *
get_decay_dec(void)
{
	char *env = getenv("QQ_DECAY_TABLE");
	char *belle_top_dir = getenv("BELLE_TOP_DIR");
	static char file[256] = {0};

	if( env ){
		strcpy(file, env);
	} else if ( belle_top_dir ){
		strcpy(file, belle_top_dir);
		strcat(file, "/share/data-files/qq98/decay.dec");
	} else {
		strcpy(file, "decay.dec");
	}
	dout(Debugout::INFO,"Ptype") <<"(Particle Class)\"decay.dec\": " << file<< std::endl;

	return file;
}


static char *
get_user_dec(void)
{
	char *env = getenv("QQ_USER_TABLE");
	static char file[256] = {0};

	if( env ){
		strcpy(file, env);
		dout(Debugout::INFO,"Ptype") <<"(Particle Class)\"user.dec\": " << file<< std::endl;
	} else {
	  dout(Debugout::INFO,"Ptype") <<"(Particle Class)\"user.dec\": (none)"<< std::endl;
	}

	return file;
}


static char *
get_evt_pdl(void)
{
	char *env = getenv("EVTGEN_PDL");
	char *belle_top_dir = getenv("BELLE_TOP_DIR");
	static char file[256] = {0};

	if( env ){
		strcpy(file, env);
	} else if ( belle_top_dir ){
		strcpy(file, belle_top_dir);
		strcat(file, "/share/data-files/BelEvtGen/evt.pdl");
	} else {
		strcpy(file, "evt.pdl");
	}
	dout(Debugout::INFO,"Ptype") <<"(Particle Class)\"evt.pdl\": " <<file<< std::endl;

	return file;
}


static void
qqlist_qq(void)
{
  dout(Debugout::INFO,"Ptype") <<"(Particle Class)creating from \"decay.dec\" and \"user.dec\" ..."<< std::endl;

	char *qq_sys_tab     = get_decay_dec();
	char *qq_user_tab    = get_user_dec();

	if( !qq2panther_( qq_sys_tab, qq_user_tab, NULL) ){
	  dout(Debugout::ERR,"Ptype") <<"(Particle Class)failure in creation"<< std::endl;
		exit(1);
	}
}


static void
qqlist_evtgen(void)
{
  dout(Debugout::INFO,"Ptype") <<"(Particle Class)creating from \"decay.dec\", \"user.dec\", and \"evt.pdl\" ..."<< std::endl;

	char *qq_sys_tab     = get_decay_dec();
	char *qq_user_tab    = get_user_dec();
	char *evtgen_sys_tab = get_evt_pdl();

	char command[1024];
	
	char *belle_top_dir = getenv("BELLE_TOP_DIR");
	if( belle_top_dir ){
	std::sprintf(command, "perl %s/share/bin/qqlist.pl -sortbyname ",
			belle_top_dir);
	} else {
		strcpy(command, "perl qqlist.pl -sortbyname ");
	}

	if( qq_sys_tab[0] ){
		char tmp[256];
	std::sprintf(tmp, "-qq=%s ", qq_sys_tab);
		strcat(command, tmp);
	}

	if( qq_user_tab[0] ){
		char tmp[256];
	std::sprintf(tmp, "-qq=%s ", qq_user_tab);
		strcat(command, tmp);
	}

	if( evtgen_sys_tab[0] ){
		char tmp[256];
	std::sprintf(tmp, "-evtgen=%s ", evtgen_sys_tab);
		strcat(command, tmp);
	}

	dout(Debugout::INFO,"Ptype") <<"(Particle Class)command: " << command<< std::endl;
	FILE *fp = popen(command, "r");

	char buf[256];
	while( fgets(buf,256,fp) ){
		char name[256], gen[256];
		int pdgid, stable;
		double mass, charge, spin, ctau, width, mass_min, mass_max;
		sscanf(buf, "%s %d %d %lf %lf %lf %lf %lf %lf %lf %s\n",
			name, &pdgid, &stable, &mass, &charge, &spin, &ctau,
			&width, &mass_min, &mass_max, gen);
		struct qq_particle *qqpdt;

		qqpdt = (struct qq_particle*)BsNewEnt(QQ_PARTICLE);

		strcpy(qqpdt->m_name, name);
		qqpdt->m_pdgid  = pdgid;
		qqpdt->m_charge = charge;
		qqpdt->m_mass   = mass;
		qqpdt->m_ctau   = ctau;
		qqpdt->m_spin   = spin;
		qqpdt->m_decay[0] = qqpdt->m_decay[1] = !stable;
	}

	fclose(fp);
}


static void
initHepPDTable(void)
{
	static bool flgFirst = true;
	if( !flgFirst ){
	  return;
	} else {
	  flgFirst = false;
		use_EVTGEN = !!getenv("BELLE_USE_EVTGEN");
	}
	

	nQQParticle = BsCouTab(QQ_PARTICLE);

	if( nQQParticle == 0 ){

	  dout(Debugout::INFO,"Ptype") <<"(Particle Class)PANTHER table for QQ not found"<< std::endl;

		if( use_EVTGEN )
			qqlist_evtgen();
		else
			qqlist_qq();

		dout(Debugout::INFO,"Ptype") <<"(Particle Class)creation completed"<< std::endl;

		nQQParticle = BsCouTab(QQ_PARTICLE);

	} else {

	  dout(Debugout::INFO,"Ptype") <<"(Particle Class)Panther table \"QQ_PARTICLE\" exists."<< std::endl;
	  dout(Debugout::INFO,"Ptype") <<"(Particle Class)This table is used to create Ptype data."<< std::endl;

	}


	// 2001/04/07 jtanaka
	// not use Panther:QQ_PARTICLE but m_QQ_PARTICLE_LIST.
	// note: QQ_PARTICLE is not always same with m_QQ_PARTICLE_LIST.
	//
	// I see.  2002/11/01, T.Higuchi
	if (NULL == (m_QQ_PARTICLE_LIST =
		     (struct qq_particle*)malloc(sizeof(struct qq_particle)*nQQParticle))) {
	  perror("$Source$:m_QQ_PARTICLE_LIST:malloc");
	  exit(1);
	}
	for( int i=1; i<=nQQParticle; i++ ){
	  struct qq_particle *qq_particle
			= (struct qq_particle*)BsGetEnt(QQ_PARTICLE, i, BBS_No_Index);
		m_QQ_PARTICLE_LIST[i-1] = *qq_particle;
	}


	if (NULL == (SortName  = (struct qq_particle**)
		     malloc(sizeof(struct qq_particle*)*nQQParticle))) {
	  perror("$Source$:SortName:malloc");
	  exit(1);
	}
	if (NULL == (SortPDGid = (struct qq_particle**)
		     malloc(sizeof(struct qq_particle*)*nQQParticle))) {
	  perror("$Source$:SortPDGid:malloc");
	  exit(1);
	}

	// save pointers
	for(int i=0; i<nQQParticle; ++i){
	  SortName[i] = SortPDGid[i] = &(m_QQ_PARTICLE_LIST[i]);
	}

	// sort by name/pdgid
	qsort(SortName ,nQQParticle,sizeof(struct qq_particle*),cmp_name );
	qsort(SortPDGid,nQQParticle,sizeof(struct qq_particle*),cmp_pdgid);
}


Ptype::Ptype(const char *name)
{
	initHepPDTable();

	struct qq_particle *key;
	struct qq_particle **tmp;

	if (NULL == (key = (struct qq_particle*)malloc(sizeof(struct qq_particle)))) {
	  perror("$Source$:key:malloc");
	  exit(1);
	}

	strcpy(key->m_name,name);

	tmp = (struct qq_particle**)bsearch(&key, SortName ,nQQParticle, sizeof(struct qq_particle*), cmp_name);
	
	free(key);

	if(!tmp){
	  dout(Debugout::ERR,"Ptype") << "(Particle Class)Particle Name = " << name << " does not exist." << std::endl;
	  m_qq = NULL;
	  m_usable = UNUSABLE;
	  return;
	}

	m_qq = *tmp;
	m_usable = USABLE;
}


Ptype::Ptype(const int pdgid)
{
	initHepPDTable();

	struct qq_particle *key;
	struct qq_particle **tmp;

	if (NULL == (key = (struct qq_particle*)malloc(sizeof(struct qq_particle)))) {
	  perror("$Source$:key:mallock");
	  exit(1);
	}

	key->m_pdgid = pdgid;
	
	tmp = (struct qq_particle**)bsearch(
			&key,SortPDGid,nQQParticle,sizeof(struct qq_particle*),cmp_pdgid);
	
	free(key);
	
	if(!tmp){
	  dout(Debugout::ERR,"Ptype") << "(Particle Class)Particle ID = " << pdgid << " does not exist." << std::endl;
	  m_qq = NULL;
	  m_usable = UNUSABLE;
	  return;
	}

	m_qq = *tmp;
	m_usable = USABLE;
}


const unsigned int 
Ptype::nDecay(void) const
{
	if( use_EVTGEN ){
		return 911;
	} else {
		if(!m_qq->m_decay[0])return 0;
		else return abs(m_qq->m_decay[1]-m_qq->m_decay[0]+1);
	}
}


Ptype 
Ptype::chargeConjugation(void) const
{
  Ptype cc(-1*(m_qq->m_pdgid));

  if(cc)return cc;
  else return *this;
}


void 
Ptype::dump(const std::string & keyword,
            const std::string & prefix) const
{
  if(keyword.find("nameonly") != std::string::npos){
    dout(Debugout::DUMP,"Ptype") << name();
    return;
  }

  bool full = false;
  if (keyword.find("full") != -1) full = true;

  dout(Debugout::DUMP,"Ptype") << prefix;
  dout(Debugout::DUMP,"Ptype") << "Ptype:";
  if(full || keyword.find("name")   != std::string::npos)dout(Debugout::DUMP,"Ptype") << " name="   << name();
  if(full || keyword.find("lund7")  != std::string::npos)dout(Debugout::DUMP,"Ptype") << " lund7="  << lund();
  if(full || keyword.find("mass")   != std::string::npos)dout(Debugout::DUMP,"Ptype") << " mass="   << mass();
  if(full || keyword.find("ctau")   != std::string::npos)dout(Debugout::DUMP,"Ptype") << " ctau="   << cTau();
  if(full || keyword.find("charge") != std::string::npos)dout(Debugout::DUMP,"Ptype") << " charge=" << charge();
  if(full || keyword.find("spin")   != std::string::npos)dout(Debugout::DUMP,"Ptype") << " spin="   << spin();
  if(!use_EVTGEN && (full || keyword.find("ndecay") != std::string::npos))dout(Debugout::DUMP,"Ptype") << " nDecay=" << nDecay();
  if(full || keyword.find("stable") != std::string::npos)dout(Debugout::DUMP,"Ptype") << " stable=" << (int)stable();
  if(full || keyword.find("return") != std::string::npos)dout(Debugout::DUMP,"Ptype") << std::endl;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

// File: Muid_mdst.cc
//
// Read MDST tables and setup Muid_mdst class.
//
// Creation: 02-Mar-2000
// Version: 19-Apr-2000
//
// $Id: Muid_mdst.cc 9944 2006-11-29 07:36:07Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.4  2001/12/13 15:31:51  katayama
// MDST_OBS
//
// Revision 1.3  2000/04/19 07:47:46  katayama
// removed message
//
// Revision 1.2  2000/04/13 12:36:22  katayama
// Added std:: to cout,cerr,endl etc.
//
// Revision 1.1  2000/04/07 21:41:47  katayama
// Muid_mdst is here
//
// Revision 1.5  2000/04/05 11:05:21  teramoto
// Modify the implementation of Muid_mdst, so that it can also read MDST
// without MDST_KLM_Mu_EX.
//
// Revision 1.4  2000/03/17 11:35:41  teramoto
// Added a status member variable/function.
//
// Revision 1.3  2000/03/15 11:15:23  teramoto
// Added a protection for the very rare cases that MDST_charged exists but
// there's no MDST_MuId.
//
// Revision 1.2  2000/03/09 01:41:28  katayama
// Compatibility with CC5.0
//
// Revision 1.1  2000/03/02 10:37:49  teramoto
// Added the Muid_mdst class, which is a user interface for accessing MDST_MuId,
// etc panther tables.
//

#include "belle.h"
#include        "mdst/Muid_mdst.h"

#include        "panther/panther.h"
#include	MDST_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


// Constructor by specifying the MDST_charged.

Muid_mdst::Muid_mdst( const Mdst_charged& mdst_charged )
{

  if( mdst_charged.muid() ){

    Mdst_muid mdst_muid( mdst_charged.muid() );

    m_level = mdst_muid.muon();
    m_quality = mdst_muid.quality();

    Mdst_klm_mu klm_mu( mdst_muid.klm() );

    m_absp_len_trk = klm_mu.absp_len_pas();
    m_absp_len_hit = klm_mu.absp_len_hit();
    m_layer_trk = klm_mu.layer_pas();
    m_layer_hit = klm_mu.layer_hit();

// Get managers.

    Mdst_klm_mu_ex_Manager &ex_mgr = Mdst_klm_mu_ex_Manager::get_manager();

    if( ex_mgr.begin() != ex_mgr.end() ){

      Mdst_klm_mu_ex_Index ex_idx = ex_mgr.index("pMDST_Charged");
      ex_idx.update();
      std::vector<Mdst_klm_mu_ex> ex_v 
	= point_from( mdst_charged.get_ID(), ex_idx );
      Mdst_klm_mu_ex &ex = ex_v.front();	// ex_v should have only one element.

      m_muon_likelihood	= ex.Muon_likelihood();
      m_pion_likelihood	= ex.Pion_likelihood();
      m_kaon_likelihood	= ex.Kaon_likelihood();
      m_miss_likelihood	= ex.Miss_likelihood();
      m_junk_likelihood	= ex.Junk_likelihood();

      m_outcome		= ex.Outcome();
      m_chi_2		= ex.Chi_2();
      m_n_hits		= ex.N_hits();
      m_n_strips	= ex.N_strips();
      m_layer_trk_brl	= ex.Layer_trk_brl();
      m_layer_hit_brl	= ex.Layer_hit_brl();
      m_layer_trk_end	= ex.Layer_trk_end();
      m_layer_hit_end	= ex.Layer_hit_end();
      m_n_layer_trk_brl	= ex.N_layer_trk_brl();
      m_n_layer_hit_brl	= ex.N_layer_hit_brl();
      m_n_layer_trk_end	= ex.N_layer_trk_end();
      m_n_layer_hit_end	= ex.N_layer_hit_end();
      m_trk_pattern	= ex.Trk_pattern();
      m_hit_pattern	= ex.Hit_pattern();

      m_status		= 1;		// Success.

    } else {

      m_muon_likelihood	= -1;
      m_pion_likelihood	= -1;
      m_kaon_likelihood	= -1;
      m_miss_likelihood	= -1;
      m_junk_likelihood	= -1;

      m_outcome		= -1;
      m_chi_2		= -1;
      m_n_hits		= -1;
      m_n_strips	= -1;
      m_layer_trk_brl	= -1;
      m_layer_hit_brl	= -1;
      m_layer_trk_end	= -1;
      m_layer_hit_end	= -1;
      m_n_layer_trk_brl	= -1;
      m_n_layer_hit_brl	= -1;
      m_n_layer_trk_end	= -1;
      m_n_layer_hit_end	= -1;
      m_trk_pattern	= 0;
      m_hit_pattern	= 0;

      m_status		= 2;		// Success.
    }

  } else {

    m_level = -2;
    m_quality = 0;

    m_absp_len_trk = -1;
    m_absp_len_hit = -1;
    m_layer_trk = -2;
    m_layer_hit = -2;

    m_muon_likelihood	= -1;
    m_pion_likelihood	= -1;
    m_kaon_likelihood	= -1;
    m_miss_likelihood	= -1;
    m_junk_likelihood	= -1;

    m_outcome		= -1;
    m_chi_2		= -1;
    m_n_hits		= -1;
    m_n_strips	= -1;
    m_layer_trk_brl	= -1;
    m_layer_hit_brl	= -1;
    m_layer_trk_end	= -1;
    m_layer_hit_end	= -1;
    m_n_layer_trk_brl	= -1;
    m_n_layer_hit_brl	= -1;
    m_n_layer_trk_end	= -1;
    m_n_layer_hit_end	= -1;
    m_trk_pattern	= 0;
    m_hit_pattern	= 0;

// dout(Debugout::DDEBUG,"Muid_mdst") << "%ERROR: Muid_mdst: failed to access the MDST_MuId table." <<std::endl;
    m_status		= 0;		// Failed.

  }

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

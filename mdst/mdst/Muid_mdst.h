// File: Muid_mdst.h
//
// Muid_mdst class.
//
// Creation: 02-Mar-2000
// Version: 09-Mar-2002
//
// $Id: Muid_mdst.h 9932 2006-11-12 14:26:53Z katayama $
//
// Revision history
//
// $Log$
// Revision 1.4  2002/03/09 10:46:07  hitoshi
// added 2 functions, Shared_KLM_hit, Muid_mdst::Shared_KLM_hit (by Teramoto).
//
// Revision 1.3  2001/12/13 15:31:51  katayama
// MDST_OBS
//
// Revision 1.2  2000/10/27 02:00:06  katayama
// check mx layer
//
// Revision 1.1  2000/04/07 21:41:45  katayama
// Muid_mdst is here
//
// Revision 1.5  2000/03/21 11:53:07  teramoto
// Fix a comment sentence, which effect doc++.
//
// Revision 1.4  2000/03/21 00:23:18  teramoto
// Likelihood for Endcap_MX_layer = 13(def) and 11, separately calculated by
// Piilonen-san. In addition, Prerejection() member function, which is actually
// the chi2==0 check, is added in the Muid_mdst class.
//
// Revision 1.3  2000/03/17 11:35:27  teramoto
// Added a status member variable/function.
//
// Revision 1.2  2000/03/04 04:29:47  teramoto
// A minor fix for making doc++ run.
//
// Revision 1.1  2000/03/02 10:37:56  teramoto
// Added the Muid_mdst class, which is a user interface for accessing MDST_MuId,
// etc panther tables.
//

#ifndef _Muid_Mdst_Flag_
#define _Muid_Mdst_Flag_
#include "belle.h"
#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



class Mdst_charged;


/**
  MDST\_MuId, MDST\_KLM\_Mu, MDST\_KLM\_Mu\_EX  Handler.
*/

class Muid_mdst {

public:

/// Constructor by specifying the MDST\_charged.
  Muid_mdst( const Mdst_charged& mdst_charged );

/// Destructor.
  virtual ~Muid_mdst(){}

// Accessors

/// Get MDST\_MuId:muon.
  inline int Level() const;

/// Get MDST\_KLM\_Mu\_EX:Muon\_likelihood.
  inline double Muon_likelihood() const;

/// Get MDST\_KLM\_Mu\_EX:Pion\_likelihood.
  inline double Pion_likelihood() const;

/// Get MDST\_KLM\_Mu\_EX:Kaon\_likelihood.
  inline double Kaon_likelihood() const;

/// Get MDST\_KLM\_Mu\_EX:Miss\_likelihood.
  inline double Miss_likelihood() const;

/// Get MDST\_KLM\_Mu\_EX:Junk\_likelihood.
  inline double Junk_likelihood() const;

/// Get MDST\_KLM\_Mu\_EX:Outcome.
  inline int Outcome() const;

/// Get MDST\_MuId:quality.
  inline int Quality() const;

/// Get MDST\_KLM\_Mu\_EX:Chi\_2.
  inline double Chi_2() const;

/// Get MDST\_KLM\_Mu\_EX:N\_hits.
  inline int N_hits() const;

/// Get MDST\_KLM\_Mu\_EX:N\_strips.
  inline int N_strips() const;

/// Get MDST\_KLM\_Mu:Absp\_len\_pas.
  inline double Absp_len_trk() const;

/// Get MDST\_KLM\_Mu:Absp\_len\_hit.
  inline double Absp_len_hit() const;

/// Get MDST\_KLM\_Mu:Layer\_pas.
  inline int Layer_trk() const;

/// Get MDST\_KLM\_Mu:Layer\_hit.
  inline int Layer_hit() const;

/// Get MDST\_KLM\_Mu\_EX:Layer\_trk\_brl.
  inline int Layer_trk_brl() const;

/// Get MDST\_KLM\_Mu\_EX:Layer\_hit\_brl.
  inline int Layer_hit_brl() const;

/// Get MDST\_KLM\_Mu\_EX:Layer\_trk\_end.
  inline int Layer_trk_end() const;

/// Get MDST\_KLM\_Mu\_EX:Layer\_hit\_end.
  inline int Layer_hit_end() const;

/// Get N\_layer\_pas.
  inline int N_layer_trk() const;

/// Get N\_layer\_hit.
  inline int N_layer_hit() const;

/// Get MDST\_KLM\_Mu\_EX:N\_layer\_trk\_brl.
  inline int N_layer_trk_brl() const;

/// Get MDST\_KLM\_Mu\_EX:N\_layer\_hit\_brl.
  inline int N_layer_hit_brl() const;

/// Get MDST\_KLM\_Mu\_EX:N\_layer\_trk\_end.
  inline int N_layer_trk_end() const;

/// Get MDST\_KLM\_Mu\_EX:N\_layer\_hit\_end.
  inline int N_layer_hit_end() const;

/// Get MDST\_KLM\_Mu\_EX:Trk\_pattern.
  inline unsigned int Trk_pattern() const;

/// Get MDST\_KLM\_Mu\_EX:Hit\_pattern.
  inline unsigned int Hit_pattern() const;

/// Get MDST access status.
  inline bool Status() const;

/// Get prerejection information: Y/N(=1/0).
  inline bool Prerejection() const;

/// Get Endcap_MX_layer=11 flag: Y/N(=1/0).
  inline bool Endcap_MX_layer_11() const;

/// Get Shared KLM hit track flag: Y/N(=1/0).
  inline bool Shared_KLM_hit() const;

private:

  int		m_level;
  double	m_muon_likelihood;
  double	m_pion_likelihood;
  double	m_kaon_likelihood;
  double	m_miss_likelihood;
  double	m_junk_likelihood;
  int		m_outcome;
  int		m_quality;
  double	m_chi_2;
  int		m_n_hits;
  int		m_n_strips;
  double	m_absp_len_trk;
  double	m_absp_len_hit;
  int		m_layer_trk;
  int		m_layer_hit;
  int		m_layer_trk_brl;
  int		m_layer_hit_brl;
  int		m_layer_trk_end;
  int		m_layer_hit_end;
  int		m_n_layer_trk_brl;
  int		m_n_layer_hit_brl;
  int		m_n_layer_trk_end;
  int		m_n_layer_hit_end;
  int		m_trk_pattern;
  int		m_hit_pattern;

  bool		m_status;

};

// Inline functions.

inline int Muid_mdst::Level() const { return m_level; }

inline double Muid_mdst::Muon_likelihood() const { return m_muon_likelihood; }

inline double Muid_mdst::Pion_likelihood() const { return m_pion_likelihood; }

inline double Muid_mdst::Kaon_likelihood() const { return m_kaon_likelihood; }

inline double Muid_mdst::Miss_likelihood() const { return m_miss_likelihood; }

inline double Muid_mdst::Junk_likelihood() const { return m_junk_likelihood; }

inline int Muid_mdst::Outcome() const { return m_outcome; }

inline int Muid_mdst::Quality() const { return m_quality; }

inline double Muid_mdst::Chi_2() const { return m_chi_2; }

inline int Muid_mdst::N_hits() const { return m_n_hits; }

inline int Muid_mdst::N_strips() const { return m_n_strips; }

inline double Muid_mdst::Absp_len_trk() const { return m_absp_len_trk; }

inline double Muid_mdst::Absp_len_hit() const { return m_absp_len_hit; }

inline int Muid_mdst::Layer_trk() const { return m_layer_trk; }

inline int Muid_mdst::Layer_hit() const { return m_layer_hit; }

inline int Muid_mdst::Layer_trk_brl() const { return m_layer_trk_brl; }

inline int Muid_mdst::Layer_hit_brl() const { return m_layer_hit_brl; }

inline int Muid_mdst::Layer_trk_end() const { return m_layer_trk_end; }

inline int Muid_mdst::Layer_hit_end() const { return m_layer_hit_end; }

inline int Muid_mdst::N_layer_trk() const { 
	return m_n_layer_trk_brl + m_n_layer_trk_end; }

inline int Muid_mdst::N_layer_hit() const { 
	return m_n_layer_hit_brl + m_n_layer_hit_end; }

inline int Muid_mdst::N_layer_trk_brl() const { return m_n_layer_trk_brl; }

inline int Muid_mdst::N_layer_hit_brl() const { return m_n_layer_hit_brl; }

inline int Muid_mdst::N_layer_trk_end() const { return m_n_layer_trk_end; }

inline int Muid_mdst::N_layer_hit_end() const { return m_n_layer_hit_end; }

inline unsigned int Muid_mdst::Trk_pattern() const { return m_trk_pattern; }

inline unsigned int Muid_mdst::Hit_pattern() const { return m_hit_pattern; }

inline bool Muid_mdst::Status() const { return m_status; }

inline bool Muid_mdst::Prerejection() const {
	return ( ( m_chi_2 == 0 ) ? 1: 0 ); }

inline bool Muid_mdst::Endcap_MX_layer_11() const {
	return ( ( m_quality & 0x400000 ) ? 1: 0 ); }

inline bool Muid_mdst::Shared_KLM_hit() const {
	return ( ( m_quality & 0x10000000 ) ? 1: 0 ); }

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
#endif /* _Muid_Mdst_Flag_ */

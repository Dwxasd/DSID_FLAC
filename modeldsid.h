// modeldsid.cpp
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/07  Hao Xu   First release 

#ifndef MODELDSID_H_
#define MODELDSID_H_

#pragma once

#include "C:\Program Files\Itasca\Flac3d500\pluginfiles\models\src\conmodel.h"

namespace models {
  class Modelbexmc : public ConstitutiveModel  {
    public:
      // Creators
      Modelbexmc();
      virtual String        getName() const { return L"bexmc"; }
      virtual String        getFullName() const { return L"bexmc"; }
      virtual UInt          getMinorVersion() const;
      virtual String        getProperties() const;
      virtual String        getStates() const;
      virtual Variant       getProperty(UInt index) const;
      virtual void          setProperty(UInt index, const Variant &p,UInt restoreVersion=0);

      //virtual Modelbexmc *clone() const { return(NEW("Modelbexmc") Modelbexmc()); }
      virtual Modelbexmc    *clone() const { return new Modelbexmc(); }
      virtual Double        getConfinedModulus() const { return(BulkB_ + Shear_*4.0/3.0); }
      virtual Double        getShearModulus() const { return(Shear_); }
      virtual Double        getBulkModulus() const { return(Bulk_T_); }
      virtual void          copy(const ConstitutiveModel *mod);
      virtual void          initialize(UByte dim,State *s);
      virtual void          run(UByte dim,State *s); 

      // Optional
      virtual bool          supportsHystereticDamping() const { return(false); }
    private:
      Double Bulk_,BulkB_,Shear_,Poisson_,Kappa_PS0_,Lambda_PS0_,MM_,MPC_,MP1_;
      Double MV_L_,MV_,EV_,EV_P_,MP_,MQ_,T_now_,Suct_,Alpha_PS_,Alpha_0_,Alpha_1_;
	  Double Alpha_2_,Alpha_3_,Beta_L_,r_L_,T_ref_,MPS0_,k_PS_,row_PS_,Kappa_PS_,Lambda_PS_;
	  Double Kappa_SP_,Kappa_SP0_,Alpha_SP_,Alpha_SS_,P_ref_,Bulk_H_,DSuct_,MPC_ST_,MPC_T_;
	  Double MPS_,Nonass_;
	  Double Bulk_m_,Kappa_m_,MV_m_,MV_Ma_,EV_m_,EV_Ma_,P0_m_,dP_,Bulk_T_;
	  Double f_SD_,f_SD0_,f_SD1_,n_SD_,f_SI_,f_SI0_,f_SI1_,n_SI_,sd_,si_;
	  Double EV_PLC_,EV_PmM_,dEV_PmM_,MSC_,S11old_,dp_LC_,dP_m_,a_1_,a_2_,alpha_m_,betha_m_;
	  Double ftanh_,f_SD2_,f_SD3_,f_SI2_,f_SI3_,cyc_,Bulk_s_,Suct_m_,Suct_O_,sigM_,sigH_,sigT_, dEVs_m_,dEVs_mtr_,PP0_,f_sc_,
		     dds11_,dds22_,dds33_,de11_,de22_,de33_,Debug_,toughTime_,dFS_copy_,fs_frac_,chem_in_,sig_miswell_,sig_tswell_,
			 dP_m_tr_,dEV_m_;
  };
}
//EOF

#endif 

// modeldsid.cpp
// subroutine:
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/07  Hao Xu   First release 

#include "modeldsid.h"
#include "C:\Program Files\Itasca\Flac3d500\pluginfiles\models\src\state.h"
#include "C:\Program Files\Itasca\Flac3d500\pluginfiles\models\src\convert.h"
#include "version.txt"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <assert.h>


#ifdef DSID_EXPORTS
int __stdcall DllMain(void *,unsigned, void *) {
    return(1);
}

extern "C" EXPORT_TAG const char *getName() {
    return "modeldsid";
}

extern "C" EXPORT_TAG unsigned getMajorVersion() {
    return MAJOR_VERSION;
}

extern "C" EXPORT_TAG unsigned getMinorVersion() {
    return MINOR_VERSION;
}

extern "C" EXPORT_TAG void *createInstance() {
    models::Modeldsid *m = new models::Modeldsid();
    return((void *)m);
}
#endif // DSID_EXPORTS

namespace models {
    static const unsigned long mShearNow    = 0x01;
    // static const unsigned long mTensionNow  = 0x02;
    static const unsigned long mShearPast   = 0x04;
    // static const unsigned long mTensionPast = 0x08;
    static const Double dC1D3 = 1.0 / 3.0;
    static const Double dC2D3 = 2.0 / 3.0;
    static const Double dC4D3 = 4.0 / 3.0;

    Modeldsid::Modeldsid(): Bulk_(0.0),BulkB_(0.0),Shear_(0.0),Poisson_(0.0),
              Kappa_PS0_(0.0),Lambda_PS0_(0.0),MM_(0.0),MPC_(0.0),MP1_(0.0),
              MV_L_(0.0),MV_(0.0),EV_(0.0),EV_P_(0.0),MP_(0.0),MQ_(0.0),
              T_now_(0.0),Suct_(0.0),Alpha_PS_(0.0),Alpha_0_(0.0),Alpha_1_(0.0),
              Alpha_2_(0.0),Alpha_3_(0.0),Beta_L_(0.0),r_L_(0.0),T_ref_(0.0),
              MPS0_(0.0),k_PS_(0.0),row_PS_(0.0),Kappa_PS_(0.0),Lambda_PS_(0.0),
              Kappa_SP_(0.0),Kappa_SP0_(0.0),Alpha_SP_(0.0),Alpha_SS_(0.0),
              P_ref_(0.0),Bulk_H_(0.0),DSuct_(0.0),MPC_ST_(0.0),MPC_T_(0.0),
              MPS_(0.0),Nonass_(1.0),
              Bulk_m_(0.0),Kappa_m_(0.0),MV_m_(0.0),MV_Ma_(0.0),
              EV_m_(0.0),EV_Ma_(0.0),P0_m_(0.0),dP_(0.0),Bulk_T_(0.0),
              f_SD_(0.0),f_SD0_(0.0),f_SD1_(0.0),n_SD_(0.0),
              f_SI_(0.0),f_SI0_(0.0),f_SI1_(0.0),n_SI_(0.0),sd_(0.0),si_(0.0),
              EV_PLC_(0.0),EV_PmM_(0.0),dEV_PmM_(0.0),MSC_(0.0),S11old_(0.0),dp_LC_(0.0),
              dP_m_(0.0),a_1_(0.0),a_2_(0.0),alpha_m_(0.0),betha_m_(1.0),
              ftanh_(0.0),f_SD2_(0.0),f_SD3_(0.0),f_SI2_(0.0),f_SI3_(0.0),cyc_(1.0),
              Bulk_s_(0.0),Suct_m_(0.0),Suct_O_(0.0),sigM_(0.0),sigH_(0.0),
              sigT_(0.0),dEVs_m_(0.0),dEVs_mtr_(0.0),PP0_(0.0),f_sc_(0.0),
              dds11_(0.0),dds22_(0.0),dds33_(0.0),de11_(0.0),de22_(0.0),de33_(0.0),
              Debug_(0.0),toughTime_(0.0),dFS_copy_(0.0),fs_frac_(1.0),chem_in_(0.0),
              sig_miswell_(0.0),sig_tswell_(0.0),dP_m_tr_(0.0),dEV_m_(0.0){
    }

    String Modeldsid::getProperties(void) const {
        return L"bulk,shear,bulk_bound,poisson,"
               L"kappa_PS0,lambda_PS0,mm,mpc,mp1,"
               L"mv_l,cv,cam_ev,camev_p,cam_cp,cq,"
               L"temperature,suction,alpha_PS,alpha_0,alpha_1,"
               L"alpha_2,alpha_3,beta_L,r_L,temp_ref,"
               L"mps0,k_PS,row_PS,kkkappa_PS,lllambda_PS,"
               L"kappa_SP,kappa_SP0,alpha_SP,alpha_SS,"
               L"p_ref,bulk_h,dsuct,mpc_ST,mpc_T,"
               L"mpmps,nonass,"
               L"Bulk_m,kappa_m,MV_m,MV_Ma,"
               L"EV_m,EV_Ma,P0_m,dP,Bulk_T,"
               L"f_SD,f_SD0,f_SD1,n_SD,"
               L"f_SI,f_SI0,f_SI1,n_SI,sd,si,"
               L"EV_PLC,EV_PmM,dEV_PmM,MSC,S11old,dp_LC,"
               L"dP_m,a_1,a_2,alpha_m,betha_m,"
               L"ftanh,f_SD2,f_SD3,f_SI2,f_SI3,cyc,"
               L"Bulk_s,Suct_m,Suct_O,sigM,sigH,"
               L"sigT,dEVs_m,dEVs_mtr,PP0,f_sc,"
               L"dds11,dds22,dds33,de11,de22,de33,"
               L"Debug,toughTime,dFS_copy,fs_frac,chem_in,"
               L"sig_miswell,sig_tswell,dP_m_tr,dEV_m";

    }

    UInt Modeldsid::getMinorVersion() const {
        return MINOR_VERSION;
    }

    String Modeldsid::getStates(void) const {
        return L"shear-n,shear-p";
    }

    Variant Modeldsid::getProperty(UInt ul) const {
        switch (ul) {
            case 1:  return(Bulk_);
            case 2:  return(Shear_);
            case 3:  return(BulkB_);
            case 4:  return(Poisson_);
//
            case 5:  return(Kappa_PS0_);
            case 6:  return(Lambda_PS0_);
            case 7:  return(MM_);
            case 8:  return(MPC_);
            case 9:  return(MP1_);
//
            case 10: return(MV_L_);
            case 11: return(MV_);
            case 12: return(EV_);
            case 13: return(EV_P_);
            case 14: return(MP_);
            case 15: return(MQ_);
//
            case 16: return(T_now_); 
            case 17: return(Suct_);
            case 18: return(Alpha_PS_);
            case 19: return(Alpha_0_);
            case 20: return(Alpha_1_);
//
            case 21: return(Alpha_2_);
            case 22: return(Alpha_3_);
            case 23: return(Beta_L_);
            case 24: return(r_L_);
            case 25: return(T_ref_);
//
            case 26: return(MPS0_);
            case 27: return(k_PS_);
            case 28: return(row_PS_);
            case 29: return(Kappa_PS_);
            case 30: return(Lambda_PS_);
//
            case 31: return(Kappa_SP_); 
            case 32: return(Kappa_SP0_);
            case 33: return(Alpha_SP_);
            case 34: return(Alpha_SS_);
//
            case 35: return(P_ref_); 
            case 36: return(Bulk_H_);
            case 37: return(DSuct_);
            case 38: return(MPC_ST_);
            case 39: return(MPC_T_);
//
            case 40: return(MPS_); 
            case 41: return(Nonass_);
//
            case 42: return(Bulk_m_);
            case 43: return(Kappa_m_);
            case 44: return(MV_m_);
            case 45: return(MV_Ma_);
//
            case 46: return(EV_m_);
            case 47: return(EV_Ma_);
            case 48: return(P0_m_);
            case 49: return(dP_);
            case 50: return(Bulk_T_);
//
            case 51: return(f_SD_);
            case 52: return(f_SD0_);
            case 53: return(f_SD1_);
            case 54: return(n_SD_);
//
            case 55: return(f_SI_);
            case 56: return(f_SI0_);
            case 57: return(f_SI1_);
            case 58: return(n_SI_);
            case 59: return(sd_);
            case 60: return(si_);
//
            case 61: return(EV_PLC_);
            case 62: return(EV_PmM_);
            case 63: return(dEV_PmM_);
            case 64: return(MSC_);
            case 65: return(S11old_);
            case 66: return(dp_LC_);
//
            case 67: return(dP_m_);
            case 68: return(a_1_);
            case 69: return(a_2_);
            case 70: return(alpha_m_);
            case 71: return(betha_m_);
//
            case 72: return(ftanh_);
            case 73: return(f_SD2_);
            case 74: return(f_SD3_);
            case 75: return(f_SI2_);
            case 76: return(f_SI3_);
            case 77: return(cyc_);
//
            case 78: return(Bulk_s_);
            case 79: return(Suct_m_);
            case 80: return(Suct_O_);
            case 81: return(sigM_);
            case 82: return(sigH_);
//
            case 83: return(sigT_);
            case 84: return(dEVs_m_);
            case 85: return(dEVs_mtr_);
            case 86: return(PP0_);
            case 87: return(f_sc_);
//
            case 88: return(dds11_);
            case 89: return(dds22_);
            case 90: return(dds33_);
            case 91: return(de11_);
            case 92: return(de22_);
            case 93: return(de33_);
//
            case 94: return(Debug_);
            case 95: return(toughTime_);
            case 96: return(dFS_copy_);
            case 97: return(fs_frac_);
            case 98: return(chem_in_);
//
            case 99: return(sig_miswell_);
            case 100: return(sig_tswell_);
            case 101: return(dP_m_tr_);
            case 102: return(dEV_m_);
        }
        return(0.0);
    }

    void Modeldsid::setProperty(UInt ul,const Variant &p,UInt restoreVersion) {
        ConstitutiveModel::setProperty(ul,p,restoreVersion);
        switch (ul) {
            case 1:  Bulk_       = p.toDouble();  break;
            case 2:  Shear_      = p.toDouble();  break;
            case 3:  BulkB_      = p.toDouble();  break;
            case 4:  Poisson_    = p.toDouble();  break;
            case 5:  Kappa_PS0_  = p.toDouble();  break;
            case 6:  Lambda_PS0_ = p.toDouble();  break;
            case 7:  MM_         = p.toDouble();  break;
            case 8:  MPC_        = p.toDouble();  break;
            case 9:  MP1_        = p.toDouble();  break;
            case 10: MV_L_       = p.toDouble();  break;
            case 11: MV_         = p.toDouble();  break;
            case 12: EV_         = p.toDouble();  break;
            case 13: EV_P_       = p.toDouble();  break;
            case 14: MP_         = p.toDouble();  break;
            case 15: MQ_         = p.toDouble();  break;
            case 16: T_now_      = p.toDouble();  break;
            case 17: Suct_       = p.toDouble();  break; 
            case 18: Alpha_PS_   = p.toDouble();  break;
            case 19: Alpha_0_    = p.toDouble();  break;
            case 20: Alpha_1_    = p.toDouble();  break;
            case 21: Alpha_2_    = p.toDouble();  break;
            case 22: Alpha_3_    = p.toDouble();  break;
            case 23: Beta_L_     = p.toDouble();  break;
            case 24: r_L_        = p.toDouble();  break;
            case 25: T_ref_      = p.toDouble();  break;
            case 26: MPS0_       = p.toDouble();  break;
            case 27: k_PS_       = p.toDouble();  break;
            case 28: row_PS_     = p.toDouble();  break;
            case 29: Kappa_PS_   = p.toDouble();  break;
            case 30: Lambda_PS_  = p.toDouble();  break;
            case 31: Kappa_SP_   = p.toDouble();  break; 
            case 32: Kappa_SP0_  = p.toDouble();  break;
            case 33: Alpha_SP_   = p.toDouble();  break;
            case 34: Alpha_SS_   = p.toDouble();  break;
            case 35: P_ref_      = p.toDouble();  break; 
            case 36: Bulk_H_     = p.toDouble();  break;
            case 37: DSuct_      = p.toDouble();  break; 
            case 38: MPC_ST_     = p.toDouble();  break; 
            case 39: MPC_T_      = p.toDouble();  break; 
            case 40: MPS_        = p.toDouble();  break; 
            case 41: Nonass_     = p.toDouble();  break;
            case 42: Bulk_m_     = p.toDouble();  break;
            case 43: Kappa_m_    = p.toDouble();  break;
            case 44: MV_m_       = p.toDouble();  break;
            case 45: MV_Ma_      = p.toDouble();  break;
            case 46: EV_m_       = p.toDouble();  break;
            case 47: EV_Ma_      = p.toDouble();  break;
            case 48: P0_m_       = p.toDouble();  break;
            case 49: dP_         = p.toDouble();  break;
            case 50: Bulk_T_     = p.toDouble();  break;
            case 51: f_SD_       = p.toDouble();  break;
            case 52: f_SD0_      = p.toDouble();  break;
            case 53: f_SD1_      = p.toDouble();  break;
            case 54: n_SD_       = p.toDouble();  break;
            case 55: f_SI_       = p.toDouble();  break;
            case 56: f_SI0_      = p.toDouble();  break;
            case 57: f_SI1_      = p.toDouble();  break;
            case 58: n_SI_       = p.toDouble();  break;
            case 59: sd_         = p.toDouble();  break;
            case 60: si_         = p.toDouble();  break;
            case 61: EV_PLC_     = p.toDouble();  break;
            case 62: EV_PmM_     = p.toDouble();  break;
            case 63: dEV_PmM_    = p.toDouble();  break;
            case 64: MSC_        = p.toDouble();  break;
            case 65: S11old_     = p.toDouble();  break;
            case 66: dp_LC_      = p.toDouble();  break;
            case 67: dP_m_       = p.toDouble();  break;
            case 68: a_1_        = p.toDouble();  break;
            case 69: a_2_        = p.toDouble();  break;
            case 70: alpha_m_    = p.toDouble();  break;
            case 71: betha_m_    = p.toDouble();  break;
            case 72: ftanh_      = p.toDouble();  break;
            case 73: f_SD2_      = p.toDouble();  break;
            case 74: f_SD3_      = p.toDouble();  break;
            case 75: f_SI2_      = p.toDouble();  break;
            case 76: f_SI3_      = p.toDouble();  break;
            case 77: cyc_        = p.toDouble();  break;
            case 78: Bulk_s_     = p.toDouble();  break;
            case 79: Suct_m_     = p.toDouble();  break;
            case 80: Suct_O_     = p.toDouble();  break;
            case 81: sigM_        = p.toDouble();  break;
            case 82: sigH_        = p.toDouble();  break;
            case 83: sigT_        = p.toDouble();  break;
            case 84: dEVs_m_       = p.toDouble();  break;
            case 85: dEVs_mtr_     = p.toDouble();  break;
            case 86: PP0_          = p.toDouble();  break;
            case 87: f_sc_         = p.toDouble();  break;
            case 88: dds11_         = p.toDouble();  break;
            case 89: dds22_         = p.toDouble();  break;
            case 90: dds33_         = p.toDouble();  break;
            case 91: de11_        = p.toDouble();  break;
            case 92: de22_        = p.toDouble();  break;
            case 93: de33_        = p.toDouble();  break;
            case 94: Debug_       = p.toDouble();  break;
            case 95: toughTime_   = p.toDouble();  break;
            case 96: dFS_copy_    = p.toDouble();  break;
            case 97: fs_frac_     = p.toDouble();  break;
            case 98: chem_in_     = p.toDouble();  break;
            case 99: sig_miswell_ = p.toDouble();  break;
            case 100: sig_tswell_ = p.toDouble();  break;
            case 101: dP_m_tr_    = p.toDouble();  break;
            case 102: dEV_m_      = p.toDouble();  break;
        }
    }

    void Modeldsid::copy(const ConstitutiveModel *m) {
        ConstitutiveModel::copy(m);
        const Modeldsid *em = dynamic_cast<const Modeldsid *>(m);
        if (!em) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        Bulk_       = em->Bulk_;
        Shear_      = em->Shear_;
        BulkB_      = em->BulkB_;
        Poisson_    = em->Poisson_;
        Kappa_PS0_  = em->Kappa_PS0_;
        Lambda_PS0_ = em->Lambda_PS0_;
        MM_         = em->MM_;
        MPC_        = em->MPC_;
        MP1_        = em->MP1_;
        MV_L_       = em->MV_L_;
        MV_         = em->MV_;
        EV_         = em->EV_;
        EV_P_       = em->EV_P_;
        MP_         = em->MP_;
        MQ_         = em->MQ_;
        T_now_      = em->T_now_;
        Suct_       = em->Suct_;
        Alpha_PS_   = em->Alpha_PS_;
        Alpha_0_    = em->Alpha_0_;
        Alpha_1_    = em->Alpha_1_;
        Alpha_2_    = em->Alpha_2_;
        Alpha_3_    = em->Alpha_3_;
        Beta_L_     = em->Beta_L_;
        r_L_        = em->r_L_;
        T_ref_      = em->T_ref_;
        MPS0_       = em->MPS0_;
        k_PS_       = em->k_PS_;
        row_PS_     = em->row_PS_;
        Kappa_PS_   = em->Kappa_PS_;
        Lambda_PS_  = em->Lambda_PS_;
        Kappa_SP_   = em->Kappa_SP_; 
        Kappa_SP0_  = em->Kappa_SP0_;
        Alpha_SP_   = em->Alpha_SP_;
        Alpha_SS_   = em->Alpha_SS_;
        P_ref_      = em->P_ref_;  
        Bulk_H_     = em->Bulk_H_;
        DSuct_      = em->DSuct_;
        MPC_ST_     = em->MPC_ST_;
        MPC_T_      = em->MPC_T_; 
        MPS_        = em->MPS_; 
        Nonass_     = em->Nonass_;
        Bulk_m_     = em->Bulk_m_;
        Kappa_m_    = em->Kappa_m_;
        MV_m_       = em->MV_m_;
        MV_Ma_      = em->MV_Ma_;
        EV_m_       = em->EV_m_;
        EV_Ma_      = em->EV_Ma_;
        P0_m_       = em->P0_m_;
        dP_         = em->dP_;
        Bulk_T_     = em->Bulk_T_;
        f_SD_       = em->f_SD_;
        f_SD0_      = em->f_SD0_;
        f_SD1_      = em->f_SD1_;
        n_SD_       = em->n_SD_;
        f_SI_       = em->f_SI_;
        f_SI0_      = em->f_SI0_;
        f_SI1_      = em->f_SI1_;
        n_SI_       = em->n_SI_;
        sd_         = em->sd_;
        si_         = em->si_;
        EV_PLC_     = em->EV_PLC_;
        EV_PmM_     = em->EV_PmM_;
        dEV_PmM_    = em->dEV_PmM_;
        MSC_        = em->MSC_;
        S11old_     = em->S11old_;
        dp_LC_     = em->dp_LC_;
        dP_m_       = em->dP_m_;
        a_1_        = em->a_1_;
        a_2_        = em->a_2_;
        alpha_m_    = em->alpha_m_;
        betha_m_    = em->betha_m_;
        ftanh_      = em->ftanh_;
        f_SD2_      = em->f_SD2_;
        f_SD3_      = em->f_SD3_;
        f_SI2_      = em->f_SI2_;
        f_SI3_      = em->f_SI3_;
        cyc_        = em->cyc_;
        Bulk_s_     = em->Bulk_s_;
        Suct_m_     = em->Suct_m_;
        Suct_O_     = em->Suct_O_;
        sigM_        = em->sigM_;
        sigH_        = em->sigH_;
        sigT_        = em->sigT_;
        dEVs_m_       = em->dEVs_m_;
        dEVs_mtr_     = em->dEVs_mtr_;
        PP0_          = em->PP0_;
        f_sc_         = em->f_sc_;
        dds11_          = em->dds11_;
        dds22_          = em->dds22_;
        dds33_          = em->dds33_;
        de11_          = em->de11_;
        de22_          = em->de22_;
        de33_          = em->de33_;
        Debug_         = em->Debug_;  
        toughTime_     = em->toughTime_;
        dFS_copy_      = em->dFS_copy_;
        fs_frac_       = em->fs_frac_;
        chem_in_       = em->chem_in_;
        sig_miswell_   = em->sig_miswell_;
        sig_tswell_    = em->sig_tswell_;
        dP_m_tr_       = em->dP_m_tr_;
        dEV_m_         = em->dEV_m_;
    }

    void Modeldsid::initialize(UByte dim,State *s) {
        ConstitutiveModel::initialize(dim,s);
        // initialize mean pressure
        if (MP1_ == 0.0) MP1_ = 1.0;
        if (MP1_ < 0.0 || MPC_ < 0.0) throw std::runtime_error("dsid: MP1_ < 0 or MPC_ < 0.0 ");
        if (MV_L_ < 1.0 || Kappa_PS0_ == 0.0) throw std::runtime_error("dsid: MV_L_ < 1 or Kappa_PS0_ = 0");
        if (Lambda_PS0_ == Kappa_PS0_) throw std::runtime_error("dsid: Lambda_PS0_ = Kappa_PS0_");
        if (BulkB_ == 0.0) throw std::runtime_error("dsid: BulkB_ = 0.0");
        //if (Poisson_ == 0.0) throw std::runtime_error("dsid: Poisson_ = 0.0, enter the value for Poisson's ratio");
        Double dP0 = MP_;           // initialized in dsid_ini_p
        if (dP0 < 0.0) {
            throw std::runtime_error("dsid: Mean net pressure is negative -- 0");
        }

        /* --- initialize specific volume --- */
        Suct_m_ = Suct_ + Suct_O_;
        if (Suct_m_ < 0.) Suct_m_ = 0.0;
        Kappa_PS_ = Kappa_PS0_ * ( 1.0 + Alpha_PS_ * Suct_m_ );    
        Lambda_PS_ = Lambda_PS0_ * ( ( 1.0 - r_L_ ) * exp( - Beta_L_ * Suct_m_ ) + r_L_ ); 
        Double DTemp_= T_now_ - T_ref_;
        MPC_T_ = MPC_ + 2.0 * ( Alpha_1_ * DTemp_ + Alpha_3_ * DTemp_ * abs( DTemp_ ) );
        Double P_exp_ = 0.;

        if ( Lambda_PS_ == Kappa_PS_ ) {
            throw std::runtime_error("dsid: Lambda_PS_ == Kappa_PS_");
        } else {
            P_exp_=(Lambda_PS0_-Kappa_PS_)/(Lambda_PS_-Kappa_PS_);
        }

        Double MPCRATIO_ = MPC_T_ / MP1_;
        if ( MPCRATIO_ > 0 ){
            MPC_ST_ = MP1_ * pow(MPCRATIO_,P_exp_);
            if ( MPC_ST_ ==0.0){
                std::cout << "  initial MPC_ST_ = " << MPC_ST_ << " MP1_ = " << MP1_ << " MPCRATIO_ = " << MPCRATIO_ << " P_exp_ = " << P_exp_ << std::endl;
                std::cout << " MPC_T_ = " << MPC_T_ << " MPC_ = " << MPC_ << std::endl;
            }
        } else {
            throw std::runtime_error("dsid: MPCRATIO_ <= 0 ");
        }

        if (MV_Ma_ == 0.0) {
                /* ==== initialize specific volumes ===== */
            MV_Ma_ = MV_L_ - Lambda_PS_ * log( MPC_ST_ / MP1_ ) + Kappa_PS_ * log( MPC_ST_ / dP0 ); // HX_12/02
        }

        if (MV_Ma_ < 0.0) throw std::runtime_error("dsid:  MV_Ma_ < 1.0, initial void ratio of macrostructure is negative");
        if (MV_m_ < 0.0) throw std::runtime_error("dsid:  MV_m_ < 1.0, initial void ratio of microstructure is negative");
        MV_ = MV_Ma_ + MV_m_ - 1.0; 
        if (MV_ < 0.0) throw std::runtime_error("dsid:  MV_ < 1.0, initial total specific volume is negative");

        /* --- initialize current bulk modulus --- */
        if (Bulk_ == 0.0) {
            //--Effective stress of the microstructure
            P0_m_ = MP_ + Suct_m_;
            Bulk_ = MV_Ma_ * dP0 / Kappa_PS_;
            //--Bulk modulus of the microstructure
                if (Kappa_m_ == 0.0 || betha_m_ == 0.0) {
                    throw std::runtime_error("dsid: Kappa_m_ = 0 or betha_m_ = 0");
                } else {
                    Bulk_m_ = a_1_ * MV_m_ * P0_m_ / Kappa_m_ + a_2_ * exp( alpha_m_ * P0_m_) / betha_m_;  // HX_11/30
                }
        }

        // mark std::cout << " ###############  Iitialization ############### " << std::endl;
        // mark std::cout << " MPC_ST_ = " <<  MPC_ST_ << " MPC_ = " <<  MPC_ << " MPC_T_ = " <<  MPC_T_ << std::endl;
        /* --- initialize hydraulic modulus --- */

        if ((P_ref_ > 0.0) && (Alpha_SP_ > 0.0)) {
            Kappa_SP_ = Kappa_SP0_*(1.0+Alpha_SP_*log(dP0 / P_ref_)) * exp(Alpha_SS_*Suct_m_) ; 
        } else {
            Kappa_SP_  = Kappa_SP0_;
        }

        if (Kappa_SP_ == 0.) {
            throw std::runtime_error("dsid: Kappa_SP_ = 0 ");
        } else {
            Bulk_H_ = MV_Ma_ * ( Suct_m_ + 0.1e6) / Kappa_SP_ ;
        }

        /* --- initialize total mechanical and hydraulic modulus --- */
        if (Bulk_m_ == 0.0) {
            if (Bulk_ == 0.0) {
                throw std::runtime_error("dsid: All Bulk Modulus = 0 ");
            } else {
                Bulk_T_ = Bulk_ ;
                Bulk_s_ = Bulk_H_ ;
            }
        } else {
            if (Bulk_ == 0.0) {
                Bulk_T_ = Bulk_m_* fs_frac_;
                Bulk_ = 1e3;
                Bulk_s_ = 1.0 / (1.0 / Bulk_H_ + 1.0 / Bulk_m_* fs_frac_);
            } else {
                Bulk_T_ = 1.0 / (1.0 / Bulk_ + 1.0 / Bulk_m_* fs_frac_);
                Bulk_s_ = 1.0 / (1.0 / Bulk_H_ + 1.0 / Bulk_m_* fs_frac_);
            }
        }

        /* --- initialize shear modulus --- */
        Shear_ = 1.5*Bulk_T_*(1.0-2.0*Poisson_)/(1.0+Poisson_);
    }

    static const UInt Pav    = 0;
    static const UInt Qav    = 1;
    static const UInt Evav   = 2;
    static const UInt EvdPav = 3;
    static const UInt Evavm  = 4;
    static const UInt EvavMa = 5;
    static const UInt EvdPmM = 6;
    static const UInt EvdPLC = 7;
    static const UInt dPav   = 8;
    static const UInt Pavbup = 9;
    static const UInt Pind   = 0;

    // The above are now indices into s->working_[] and s->iworkingp[].
    // This is necessary for thread safety.

    //--------------------------------------------------------------------------
    void Modeldsid::run(UByte dim,State *s) {
        ConstitutiveModel::run(dim,s);
        /* --- plasticity indicator:                                  */
        if (s->state_ & mShearNow) s->state_ |= mShearPast;
        s->state_ &= ~mShearNow;
        /* --- initialize stacks --- */
        if (!s->sub_zone_) {
            s->iworking_[Pind]  = 0;
            s->working_[Pav]    = 0.0;
            s->working_[Qav]    = 0.0;
            s->working_[Evav]   = 0.0;
            s->working_[EvdPav] = 0.0;
            s->working_[Evavm]  = 0.0;
            s->working_[EvavMa] = 0.0;
            s->working_[EvdPmM] = 0.0;
            s->working_[EvdPLC] = 0.0;
            s->working_[dPav]   = 0.0;
            s->working_[Pavbup]   = - dC1D3 * (s->stnS_.rs11() + s->stnS_.rs22() + s->stnS_.rs33());
        }

        /* --- trial elastic stresses --- */
        Double dA1 = Bulk_T_ + dC4D3 * Shear_;
        Double dA2 = Bulk_T_ - dC2D3 * Shear_;
        Double dG2 = 2.0 * Shear_;

        s->stnS_.rs11() += s->stnE_.s11() * dA1 + (s->stnE_.s22() + s->stnE_.s33()) * dA2;
        s->stnS_.rs22() += s->stnE_.s22() * dA1 + (s->stnE_.s11() + s->stnE_.s33()) * dA2;
        s->stnS_.rs33() += s->stnE_.s33() * dA1 + (s->stnE_.s11() + s->stnE_.s22()) * dA2;
        s->stnS_.rs12() += s->stnE_.s12() * dG2;
        s->stnS_.rs13() += s->stnE_.s13() * dG2;
        s->stnS_.rs23() += s->stnE_.s23() * dG2;

        Double dS11 = s->stnE_.s11() * dA1 + (s->stnE_.s22() + s->stnE_.s33()) * dA2;
        Double dS22 = s->stnE_.s22() * dA1 + (s->stnE_.s11() + s->stnE_.s33()) * dA2;
        Double dS33 = s->stnE_.s33() * dA1 + (s->stnE_.s11() + s->stnE_.s22()) * dA2;
        dds11_ = dS11;
        dds22_ = dS22;
        dds33_ = dS33;
        de11_ = s->stnE_.s11();
        de22_ = s->stnE_.s22();
        de33_ = s->stnE_.s33();

        dP_ = - (dS11 + dS22 + dS33) * dC1D3; // compression +, tensil -

        // default settings, altered below if model is found to be failing
        s->viscous_ = true;  // allow viscous strains

        if (!canFail()) return;

        Double sMeanOld = s->stnS_.s11() + s->stnS_.s22() + s->stnS_.s33();

        /* --- mean pressure --- */
        Double dPVal = - (sMeanOld) * dC1D3; //compression +, tensil -

        //dP_m_ += dP_ ; //compression +, tensil -
        dP_m_tr_ = dP_m_ + dP_; //compression +, tensil -
        //if (s->sub_zone_  == 0){
        //std::cout << std::endl << "=======================  " << s->sub_zone_ <<"   ========================================" << std::endl ;
        //std::cout << " dP_m_ = " << dP_m_  << " dP_ = " << dP_ << " dP_m_tr_ = "<<dP_m_tr_<< std::endl;
        //}
        /* --- deviatoric stresses --- */
        Double dDS11 = s->stnS_.s11() + dPVal;
        Double dDS22 = s->stnS_.s22() + dPVal;
        Double dDS33 = s->stnS_.s33() + dPVal;
        Double dDS12 = s->stnS_.s12();
        Double dDS13 = s->stnS_.s13();
        Double dDS23 = s->stnS_.s23();

        /* --- second deviatoric stress invariant --- */
        Double dJ2   = 0.5*(dDS11*dDS11+dDS22*dDS22+dDS33*dDS33)+s->stnS_.s12()*s->stnS_.s12()+s->stnS_.s13()*s->stnS_.s13()+s->stnS_.s23()*s->stnS_.s23();
        Double dQVal = sqrt(3.*dJ2);

        /* --- Suction and temperature dependent dsid parameters --- */
        Double DTemp_ = T_now_ - T_ref_;
        MPC_T_ = MPC_ + 2.0 * (Alpha_1_*DTemp_ + Alpha_3_*DTemp_*abs(DTemp_));
        Double P_exp_=(Lambda_PS0_-Kappa_PS_)/(Lambda_PS_-Kappa_PS_); // HX_12/02
        Double MPCRATIO_ = MPC_T_ / MP1_;
        MPC_ST_ = MP1_ * pow(MPCRATIO_,P_exp_);
        MPS_ = k_PS_ * Suct_m_; // HX_11/30

        /* --- dsid yield criterion --- */
        Double dM2 = MM_*MM_;
        Double dFS = dQVal*dQVal + dM2 * (dPVal+MPS_) * (dPVal - MPC_ST_);
        Double dEVVal  = s->stnE_.s11() + s->stnE_.s22() + s->stnE_.s33(); //compression -, tensile +

        // ------------------------------------DEBUG--------------------------------------------------------------------
        dEV_m_ = 0.0;

        Double dPNew   = dPVal;
        Double dPf     = 0.0;
        Double dPf2     = 0.0;
        Double dES = 0.;
        Double dPNew2   = 0.0;
        Double dQNew   = dQVal;
        Double dEvdPVal = 0.0;
        Double dEvdPLC = 0.0;
        Double dEV_Ma_ = - dEVVal;
        dEV_PmM_ = 0.0;

        Double dAlam  = 0.0;
        Double dSc_copy = 0.0;
        Double dEvdPVal_copy_1 = 0.0;

        //Double itest = 0.0;
        Double dES_copy = 0.0;

        /* --- detect small pc --- */
        if (MPC_ < MP1_*1.e-3) {
            s->iworking_[Pind] = 1;
            s->state_ |= mShearNow;
            s->stnS_.rs11() = 0.0;
            s->stnS_.rs22() = 0.0;
            s->stnS_.rs33() = 0.0;
            s->stnS_.rs12() = 0.0;
            s->stnS_.rs13() = 0.0;
            s->stnS_.rs23() = 0.0;
            dEvdPVal = - dEVVal;   //compression +, tensile -
        } else {
            if (dFS > 0.0) {
                /* ---    shear failure --- */
                s->iworking_[Pind] = 1;
                s->state_ |= mShearNow;
                if (Nonass_ <= 0.0) Nonass_ = 1.0;
                Double dSa    = 6.0 * Shear_ * dQVal * Nonass_;  // 3*G*cb*alpha=6*G*q*alpha
                Double dSc    = dM2 * (2.0 * dPVal - MPC_ST_ + MPS_);
                dSc_copy = dSc;
                Double dSb    = Bulk_T_ * dSc;
                Double dBa    = dSa*dSa + dM2 * dSb*dSb;
                Double dBb    = dSa * dQVal + 0.5 * dSb * dSc;
                Double dBc    = dFS;

                //Double dAlam  = 0.0;
                Double dAlam1 = 0.0;
                if (dBa != 0.0) {
                    Double dBoa   = dBb / dBa;
                    Double dVal   = dBoa*dBoa - dBc / dBa;

                    if (dVal < 0.0) {
                        std::cout << " " << std::endl;
                        std::cout << " dPVal = " << dPVal << std::endl;
                        std::cout << " MPC_ST_ = " << MPC_ST_ << std::endl;
                        std::cout << " MPC_ = " << MPC_ << std::endl;
                        throw std::runtime_error("dsid:  Yield envelope cannot be reached -- dVal < 0.0");
                    }

                    Double dVal1   = sqrt(dVal);
                    dAlam  = dBoa + dVal1;
                    dAlam1 = dBoa - dVal1;

                } else {
                    if (dBb != 0.0) {
                        dAlam  = 0.5 * dBc / dBb;
                        dAlam1 = dAlam;
                    } else throw std::runtime_error("dsid:  Yield envelope cannot be reached -- dBb = 0.0");
                } 

                if (fabs(dAlam1) < fabs(dAlam)) dAlam = dAlam1;
                if (dAlam < 0.0) dAlam = 0.0;
                dQNew  = dQVal - dSa * dAlam;
                dPNew  = dPVal - dSb * dAlam; //compression +, tensile - new mean pressure

                dEvdPVal = dAlam * dSc;     // HX_11/30 //compression +, tensile -
                dEvdPLC = dAlam * dSc;   // HX_12/01 //compression +, tensile -
                dEvdPVal_copy_1 = dEvdPVal;
                if (s->sub_zone_  == 0){
                std::cout << std::endl << "=============  yield surface  ================" << std::endl ;
                std::cout << "dPVal = " << dPVal << " dPNew = " << dPNew  << "dP_ = " << dP_ << "dP_m_ = "<<dP_m_<< std::endl;
                }
            }  
        }

        //========================================================================================

        //    MICRO-SCTRUCTURE  -1-

        //========================================================================================

        f_SI_ = f_SD_ = 0.0;

        if ( Bulk_m_ != 0. ) {
            //Double p_r = (dQVal * dQVal)/dM2/(dPVal + MPS_) + dPVal;○
            Double p_r = (dQNew * dQNew)/dM2/(dPNew + MPS_) + dPNew;
            Double PP0 = p_r / MPC_ST_;
            PP0_ = PP0;
            dES = dPNew - s->working_[Pavbup]+DSuct_;

            //--------------------------------------------------------------------------------------------
            if (dES > 0.) {
                /* --- microstructural compression --- fc*/
                s->iworking_[Pind] = 1;        
                if (PP0 > 0.0) {
                    f_SD_ = (1.0 - ftanh_) * (f_SD0_ + f_SD1_ * pow(PP0,n_SD_))+ ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                } else {
                    f_SD_ = (1.0 - ftanh_) * f_SD0_ + ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                }

                f_SI_ =0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_ = (dPNew - s->working_[Pavbup]) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_ += dEVs_m_;
                dEV_PmM_ = f_SD_ * dEV_m_;

            //--------------------------------------------------------------------------------------------
            } else if  (dES < 0.) {     
                /* --- microstructural swelling --- fs*/
                s->iworking_[Pind] = 1;
                if (PP0 >= 1.0) {
                    f_SI_ = (1.0 - ftanh_) * f_SI0_ + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                } else {
                    f_SI_ = (1.0 - ftanh_) * (f_SI0_ + f_SI1_ * pow(1-PP0,n_SI_)) + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                }
                f_SD_ = 0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_ = (dPNew-s->working_[Pavbup]) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_ += dEVs_m_;
                dEV_PmM_ = f_SI_ * dEV_m_;
            } else {
                dEV_m_ = 0.0;
                dEV_PmM_ = 0.0;
            }
            dPf = dEV_PmM_ * Bulk_T_;
            dPNew2 = dPNew - dPf;
            dES_copy = dES;
        }

        if (s->sub_zone_  == 0){
            std::cout << std::endl << "========================================" << std::endl ;
            std::cout << " ***************** 1 dES = " << dES  << std::endl;
            std::cout << " dPf = " << dPf << " dES1 = " << dPNew-s->working_[Pavbup]+DSuct_   << " dP_m_ = "<<dP_m_<< " dP_ = "<<dP_<<std::endl;
            std::cout << " dEV_m = " << dEV_m_ << " P =" << s->working_[Pavbup] << " dPNew = " << dPNew << std::endl;
            std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
            std::cout << " dEVs_m_ = " <<dEVs_m_ << " DSuct_ = " << DSuct_ << " Suct_ = " << Suct_ << std::endl;
        }

        Double dEV_m_2=0.0;
        Double dEV_PmM_2 =0.0;

        if ( Bulk_m_ != 0. ) {
            Double p_r = (dQNew * dQNew)/dM2/(dPNew2 + MPS_) + dPNew2;
            Double PP0 = p_r / MPC_ST_;
            std::cout << " dPf2 = " << dPf2 << " dES2 = " << dPNew2-s->working_[Pavbup]+DSuct_   << " dP_m_ = "<< dP_m_ << " dP_ = "<<dP_<< " PP0 = " << PP0_ <<std::endl;
            PP0_ = PP0;
            dES = dPNew2-s->working_[Pavbup]+DSuct_;

            //--------------------------------------------------------------------------------------------
            if (dES > 0.) {
                /* --- microstructural compression --- fc*/
                s->iworking_[Pind] = 1;        
                if (PP0 > 0.0) {
                    f_SD_ = (1.0 - ftanh_) * (f_SD0_ + f_SD1_ * pow(PP0,n_SD_))+ ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                } else {
                    f_SD_ = (1.0 - ftanh_) * f_SD0_ + ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                }

                f_SI_=0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_2 = ((dPNew2-s->working_[Pavbup])) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_2 += dEVs_m_;
                dEV_PmM_2 = f_SD_ * dEV_m_2;
                dPf2 = dEV_PmM_2 * Bulk_T_;

            //--------------------------------------------------------------------------------------------
            } else if  (dES < 0.) {
                /* --- microstructural swelling --- fs*/
                s->iworking_[Pind] = 1;

                if (PP0 >= 1.0) {
                    f_SI_ = (1.0 - ftanh_) * f_SI0_ + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                } else {
                    f_SI_ = (1.0 - ftanh_) * (f_SI0_ + f_SI1_ * pow(1-PP0,n_SI_)) + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                }

                f_SD_ = 0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_2 = ((dPNew2-s->working_[Pavbup])) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_2 += dEVs_m_;
                dEV_PmM_2 = f_SI_ * dEV_m_2;
                dPf2 = dEV_PmM_2 * Bulk_T_;//compression +, tensil -
            } else {
                dEV_m_2 = 0.0;
                dEV_PmM_2 = 0.0;
                dPf2 = 0.0;
            }

            dPNew -= (dPf+dPf2)/2.;

        }

        if (s->sub_zone_  == 0){
            std::cout << std::endl << "--------------------------------" << std::endl ;
            std::cout << " ***************** 2 dES = " << dES  << std::endl;
            std::cout << " dPf2 = " << dPf2 << " dES2 = " << dES_copy   << " dP_m_ = "<< dP_m_ << " dP_ = "<<dP_<< " PP0 = " << PP0_ <<std::endl;
            std::cout << " dEV_m2 = " << dEV_m_2  << " dPNew2=" << dPNew2  << std::endl;
            std::cout << " sigH_ = " << sigH_ << std::endl;
            std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;

        }

        f_SI_ = f_SD_ = 0.0;
        if ( Bulk_m_ != 0. ) {
            Double p_r = (dQNew * dQNew)/dM2/(dPNew + MPS_) + dPNew;
            Double PP0 = p_r / MPC_ST_;
            PP0_ = PP0;
            dES = dPNew-s->working_[Pavbup]+DSuct_;

            //--------------------------------------------------------------------------------------------
            if (dES > 0.) {
                /* --- microstructural compression --- fc*/
                s->iworking_[Pind] = 1;        
                if (PP0 > 0.0) {
                    f_SD_ = (1.0 - ftanh_) * (f_SD0_ + f_SD1_ * pow(PP0,n_SD_))+ ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                } else {
                    f_SD_ = (1.0 - ftanh_) * f_SD0_ + ftanh_ * (f_SD0_ + f_SD1_ * tanh(f_SD2_ * (PP0 - f_SD3_)));
                }
                f_SI_ =0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_ = (dPNew-s->working_[Pavbup]) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_ += dEVs_m_;
                dEV_PmM_ = f_SD_ * dEV_m_;

            //--------------------------------------------------------------------------------------------
            } else if  (dES < 0.) {
            //SI Yield Surface
                s->iworking_[Pind] = 1;

                if (PP0 >= 1.0) {
                    f_SI_ = (1.0 - ftanh_) * f_SI0_ + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                } else {
                    f_SI_ = (1.0 - ftanh_) * (f_SI0_ + f_SI1_ * pow(1-PP0,n_SI_)) + ftanh_ * (f_SI0_ - f_SI1_ * tanh(f_SI2_ * (PP0 - f_SI3_)));
                }
                f_SD_ = 0.0;
                std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
                dEV_m_ = (dPNew-s->working_[Pavbup]) / Bulk_m_* fs_frac_;
                if (dEVs_mtr_ == dEVs_m_) dEV_m_ += dEVs_m_;
                dEV_PmM_ = f_SI_ * dEV_m_;

            } else {
                dEV_PmM_ = 0.0;
            }
            dEV_Ma_ -= dEV_m_;
            dEvdPVal += dEV_PmM_; 
            dES_copy = dES;
        }

        if (abs(f_SD_) > 0.0) {
            f_sc_ = f_SD_;
        } else {
            f_sc_ = f_SI_;
        }

        if (dQVal == 0.0) {
            Double dAdd = dPNew - dPVal; // dAdd: compression -, tensile +
            s->stnS_.rs11() -= dAdd;
            s->stnS_.rs22() -= dAdd;
            s->stnS_.rs33() -= dAdd;
        } else {
            s->stnS_.rs11() = (dDS11 / dQVal) * dQNew - dPNew; //compression -, tensile +
            s->stnS_.rs22() = (dDS22 / dQVal) * dQNew - dPNew;
            s->stnS_.rs33() = (dDS33 / dQVal) * dQNew - dPNew;
            s->stnS_.rs12() = (dDS12 / dQVal) * dQNew;
            s->stnS_.rs23() = (dDS23 / dQVal) * dQNew;
            s->stnS_.rs13() = (dDS13 / dQVal) * dQNew;
        }

        if (s->sub_zone_  == 0){
            std::cout << std::endl << "--------------------------------" << std::endl ;
            std::cout << " ***************** 3 dES = " << dES  << std::endl;
            std::cout << " dPf3 = " << dPf2 << " dES3 = " <<dES_copy   << " dP_m_ = "<< dP_m_ << " dP_ = "<<dP_<< " PP0 = " << PP0_ <<std::endl;
            std::cout << " dEV_m3 = " << dEV_m_  << " dPNew3=" << dPNew  << std::endl;
            std::cout << " f_si  = " << f_SI_ << " f_sd = " << f_SD_ << std::endl;
            std::cout << "Press Enter key to continue..." << std::endl ;
            std::cin.ignore();
        }

        /* update zone stresses and strains --- */
        Double dVol = s->getSubZoneVolume();
        //std::cout << " dPNew = " << dPNew <<std::endl;
        s->working_[Pav]   += dPNew * dVol;    
        s->working_[Qav]   += dQNew * dVol;
        s->working_[Evav]  += dEVVal * dVol;
        s->working_[EvdPav] += dEvdPVal * dVol;
        s->working_[Evavm]  += dEV_m_ * dVol;
        s->working_[EvavMa] += dEV_Ma_ * dVol;  
        s->working_[EvdPmM] += dEV_PmM_ * dVol;    
        s->working_[EvdPLC] += dEvdPLC * dVol;
        s->working_[dPav]   += (dPNew - s->working_[Pavbup]) * dVol;

        dFS_copy_ = dFS;

        /* --- the last zone has been processed, update parameters: --- */
        if (s->sub_zone_==s->total_sub_zones_-1) {
            dVol = 1.0 / s->getZoneVolume();

            if (s->overlay_==2) dVol *= 0.5;
            s->working_[Pav]   = s->working_[Pav] * dVol;
            s->working_[Qav]   = s->working_[Qav] * dVol;
            s->working_[Evav]  = s->working_[Evav] * dVol;
            s->working_[EvdPav] = s->working_[EvdPav] * dVol;
            s->working_[Evavm] = s->working_[Evavm] * dVol;
            s->working_[EvavMa] = s->working_[EvavMa] * dVol;
            s->working_[EvdPmM] = s->working_[EvdPmM] * dVol;    
            s->working_[EvdPLC] = s->working_[EvdPLC] *dVol;
            s->working_[dPav]   = s->working_[dPav] *dVol;
            s->working_[Pavbup] = s->working_[Pav];

            MP_    = s->working_[Pav];
            MQ_    = s->working_[Qav];
            EV_   += s->working_[Evav];
            EV_P_ += s->working_[EvdPav];
            EV_m_ += s->working_[Evavm];
            EV_Ma_ += s->working_[EvavMa];
            EV_PmM_ += s->working_[EvdPmM];
            EV_PLC_ += s->working_[EvdPLC];
            dEVs_m_ += s->working_[Evavm];
            dP_m_ += s->working_[dPav];

            // ------------------ microstructural swelling stress ----------------------
            //sig_miswell_ = s->working_[EvdPmM] * Bulk_T_;
            sig_miswell_ += s->working_[EvdPmM]* Bulk_T_;
            S11old_ = EV_m_;
            sig_tswell_  = sig_miswell_ + sigH_;

            /* ---                  specific dVolume --- */
            Double dVold  = MV_;
            MV_    = dVold * (1.0 + s->working_[Evav]);
            if (MV_ < 0.0) throw std::runtime_error("dsid:  Current void ratio is negative, MV_ < 0 -- 1");
            if ( Bulk_m_ != 0. ) {
                Double dVmold = MV_m_;
                MV_m_ = dVmold - s->working_[Evavm] * dVold;

                if (MV_m_ < 0.0) {
                    std::cout << " " << std::endl;
                    throw std::runtime_error("dsid:  Specific volume of microstructure is < 0 -- 1");
                }
            }

            Double dVMaold = MV_Ma_;
            MV_Ma_ = dVMaold - s->working_[EvavMa] * dVold;

            if (MV_Ma_ < 0.0) {
                std::cout << " " << std::endl;
                std::cout << " dVMaold = " << dVMaold  << std::endl;
                std::cout << " s->working_[EvavMa] = " << s->working_[EvavMa] << std::endl;
                std::cout << " dEV_Ma_ = " << dEV_Ma_ << std::endl;
                throw std::runtime_error("dsid:  Specific volume of macrostructure is < 0 -- 1, 0");
            }

            /* ---                  bulk modulus --- */
            Bulk_  = MV_Ma_ * s->working_[Pav] / Kappa_PS_;

            //if (Bulk_ > BulkB_) throw std::runtime_error("dsid:  Bulk modulus bound must be increased for stability");
            /* ---                  update bulk modulus of microstructure --- */
            Suct_m_ = Suct_ + Suct_O_;
            if (Suct_m_ < 0.) Suct_m_ =0. ;
            P0_m_ = MP_ + Suct_m_;
            Bulk_m_  = a_1_ * MV_m_ * P0_m_ / Kappa_m_ + a_2_ * exp( alpha_m_ * P0_m_) / betha_m_; // HX_11/30
            Kappa_PS_ = Kappa_PS0_*(1.0+Alpha_PS_*Suct_m_);    // HX_11/30
            Lambda_PS_ = Lambda_PS0_*((1.0-r_L_)*exp(-Beta_L_*Suct_m_)+r_L_);
            /* ---                  cap pressure --- */
            if (s->iworking_[Pind] == 1) {
                //Double MPC_new = MPC_ * (1.0 + s->working_[EvdPav] * dVMaold/(Lambda_PS0_ - Kappa_PS_)); //lloret's equation
                Double MPC_new = MPC_ * (1.0 + s->working_[EvdPav] * dVold/(Lambda_PS0_ - Kappa_PS_));//sanchez's equation
                if (MPC_new >= MP1_ ) {
                    MPC_  = MPC_new;
                    //Bulk_ = 0.0;
                } else {
                    MPC_  = MP1_;
                }
            }
            ///* --- dsid update hydraulic bulk modulus for next step ---*/ 
            if ( s->working_[Pav] / P_ref_ >= 1.0) {
                Kappa_SP_ = Kappa_SP0_*(1.0+Alpha_SP_*log(s->working_[Pav] / P_ref_)) * exp(Alpha_SS_*Suct_m_) ; // HX_11/30
            } else {
                Kappa_SP_ = Kappa_SP0_ * exp(Alpha_SS_*Suct_m_) ; // HX_11/30
            }

            Bulk_H_ = MV_Ma_ * (Suct_m_+0.1e6) / Kappa_SP_;
        if (Bulk_m_ == 0.0) {
            if (Bulk_ == 0.0) {
                throw std::runtime_error("dsid: All Bulk Modulus = 0 ");
            } else {
                Bulk_T_ = Bulk_ ;
                Bulk_s_ = Bulk_H_ ;
            }
        } else {
            if (Bulk_ == 0.0) {
                Bulk_T_ = Bulk_m_ * fs_frac_;
                Bulk_ = 1e3;
                Bulk_s_ = 1.0 / (1.0 / Bulk_H_ + 1.0 / Bulk_m_* fs_frac_);
            } else {
                Bulk_T_ = 1.0 / (1.0 / Bulk_ + 1.0 / Bulk_m_* fs_frac_);
                Bulk_s_ = 1.0 / (1.0 / Bulk_H_ + 1.0 / Bulk_m_* fs_frac_);
            }
        }

            /* ---                  shear modulus --- */
            Shear_ = 1.5*Bulk_T_*(1.0-2.0*Poisson_)/(1.0+Poisson_);
        }

        if (s->state_ & mShearNow) {
            s->viscous_ = false; // inhibit viscous strains
        } else {
            s->viscous_ = true;  // allow viscous strains
        }
    }
}

// EOF

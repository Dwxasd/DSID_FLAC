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

    Modeldsid::Modeldsid(): Bulk_(0.0),BulkB_(0.0),Shear_(0.0),Poisson_(0.0),Bulk0_(0.0),
               Poisson0_(0.0),a1_(0.0),a2_(0.0),a3_(0.0),a4_(0.0),
               C0_(0.0),C1_(0.0),alpha_(0.0),Debug_(0.0) {   }

    String Modeldsid::getProperties(void) const {
        return L"bulk,shear,bulk_bound,poisson,bulk0"
               L"poisson0,a1,a2,a3,a4,"
               L"c0,c1,alpha,debug";

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
            case 5:  return(Bulk0_);
//
            case 6:  return(Poisson0_);
            case 7:  return(a1_);
            case 8:  return(a2_);
            case 9:  return(a3_);
            case 10: return(a4_);
//
            case 11: return(C0_);
            case 12: return(C1_);
            case 13: return(alpha_);
            case 14: return(Debug_);
        }
        return(0.0);
    }

    void Modeldsid::setProperty(UInt ul,const Variant &p,UInt restoreVersion) {
        ConstitutiveModel::setProperty(ul,p,restoreVersion);
        switch (ul) {
            case  1: Bulk_        = p.toDouble();  break;
            case  2: Shear_       = p.toDouble();  break;
            case  3: BulkB_       = p.toDouble();  break;
            case  4: Poisson_     = p.toDouble();  break;
            case  5: Bulk0_       = p.toDouble();  break;
            case  6: Poisson0_    = p.toDouble();  break;
            case  7: a1_          = p.toDouble();  break;
            case  8: a2_          = p.toDouble();  break;
            case  9: a3_          = p.toDouble();  break;
            case 10: a4_          = p.toDouble();  break;
            case 11: C0_          = p.toDouble();  break;
            case 12: C1_          = p.toDouble();  break;
            case 13: alpha_       = p.toDouble();  break;
            case 14: Debug_       = p.toDouble();  break;
        }
    }

    void Modeldsid::copy(const ConstitutiveModel *m) {
        ConstitutiveModel::copy(m);
        const Modeldsid *em = dynamic_cast<const Modeldsid *>(m);
        if (!em) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        Bulk_       = em->Bulk_    ;
        Shear_      = em->Shear_   ;
        BulkB_      = em->BulkB_   ;
        Poisson_    = em->Poisson_ ;
        Bulk0_      = em->Bulk0_   ;
        Poisson0_   = em->Poisson0_;
        a1_         = em->a1_      ;
        a2_         = em->a2_      ;
        a3_         = em->a3_      ;
        a4_         = em->a4_      ;
        C0_         = em->C0_      ;
        C1_         = em->C1_      ;
        alpha_      = em->alpha_   ;
        Debug_      = em->Debug_   ;
    }

    void Modeldsid::initialize(UByte dim,State *s) {
        ConstitutiveModel::initialize(dim,s);
        // initialize mean pressure

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

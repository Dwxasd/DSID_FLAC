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
#include <complex>
#include <assert.h>
#include "../mathLib/arithmetic.h"
#include "../mathLib/r1Tensor.h"
#include "../mathLib/r2Tensor.h"
#include "../mathLib/r3Tensor.h"
#include "../mathLib/errInfo.h"
#include "../mathLib/gaussj.h"
#include "../mathLib/eigen.h"


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

void effectiveStiffness(r2Tensor<double> &Matdom, const r1Tensor<double> &Omega, const double & E0_,
                        const double & Poisson0_, const double & a1_, const double & a2_,
                        const double & a3_, const double & a4_, const double & C0_,
                        const double & C1_, const double & alpha_);

void damageFunction(double & fd, const r1Tensor<double> &Sigma, const r1Tensor<double> &Omega, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_, const int & ioptfd);

void matP1(r2Tensor<double> &P1, const r1Tensor<double> &Sigma);

void matP2(r2Tensor<double> &P2, const r1Tensor<double> &Sigma);

void dY_dSigFunction(const r1Tensor<double> &Sigma, r2Tensor<double> &dY_dSig,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_);

// SUB FD_LAM
void cuttingPlaneMethod(const r1Tensor<double> &Omega, const r1Tensor<double> &Sigma, r2tensor<double> &Matdom,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_, double & H0,
             double & Hp, r1Tensor<double> &dG_dY, r1Tensor<double> &df_dSig,
             r1Tensor<double> &temp, const int & iopt);

void Mat_dS_dOmega(r3Tensor<double> &dS_dO, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_);



namespace models {
    static const unsigned long mShearNow    = 0x01;
    // static const unsigned long mTensionNow  = 0x02;
    static const unsigned long mShearPast   = 0x04;
    // static const unsigned long mTensionPast = 0x08;
    static const Double dC1D3 = 1.0 / 3.0;
    static const Double dC2D3 = 2.0 / 3.0;
    static const Double dC4D3 = 4.0 / 3.0;

    Modeldsid::Modeldsid(): Bulk_(0.0),BulkB_(0.0),Shear_(0.0),Poisson_(0.0),E0_(0.0),
               Poisson0_(0.0),a1_(0.0),a2_(0.0),a3_(0.0),a4_(0.0),
               C0_(0.0),C1_(0.0),alpha_(0.0),Debug_(0.0),Omega_00_(0.0),
               Omega_11_(0.0),Omega_22_(0.0),Omega_01_(0.0),Omega_12_(0.0),Omega_20_(0.0),
               Epsid_00_(0.0),Epsid_11_(0.0),Epsid_22_(0.0),Epsid_01_(0.0),Epsid_12_(0.0),
               Epsid_20_(0.0) {   }

    String Modeldsid::getProperties(void) const {
        return L"bulk,shear,bulk_bound,poisson,E0"
               L"poisson0,a1,a2,a3,a4,"
               L"c0,c1,alpha,debug,omega00,"
               L"omega11,omega22,omega01,omega12,omega20,"
               L"epsid00,epsid11,epsid22,epsid01,epsid12,"
               L"epsid20";

    }

    UInt Modeldsid::getMinorVersion() const {
        return MINOR_VERSION;
    }

    String Modeldsid::getStates(void) const {
        return L"shear-n,shear-p";
    }

    Variant Modeldsid::getProperty(UInt ul) const {
        switch (ul) {
            case 1:  return(Bulk_)    ;
            case 2:  return(Shear_)   ;
            case 3:  return(BulkB_)   ;
            case 4:  return(Poisson_) ;
            case 5:  return(E0_)   ;
//
            case 6:  return(Poisson0_);
            case 7:  return(a1_)      ;
            case 8:  return(a2_)      ;
            case 9:  return(a3_)      ;
            case 10: return(a4_)      ;
//
            case 11: return(C0_)      ;
            case 12: return(C1_)      ;
            case 13: return(alpha_)   ;
            case 14: return(Debug_)   ;
            case 15: return(Omega_00_);
//
            case 16: return(Omega_11_);
            case 17: return(Omega_22_);
            case 18: return(Omega_01_);
            case 19: return(Omega_12_);
            case 20: return(Omega_20_);
//
            case 21: return(Epsid_00_);
            case 22: return(Epsid_11_);
            case 23: return(Epsid_22_);
            case 24: return(Epsid_01_);
            case 25: return(Epsid_12_);
//
            case 26: return(Epsid_20_);
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
            case  5: E0_          = p.toDouble();  break;
            case  6: Poisson0_    = p.toDouble();  break;
            case  7: a1_          = p.toDouble();  break;
            case  8: a2_          = p.toDouble();  break;
            case  9: a3_          = p.toDouble();  break;
            case 10: a4_          = p.toDouble();  break;
            case 11: C0_          = p.toDouble();  break;
            case 12: C1_          = p.toDouble();  break;
            case 13: alpha_       = p.toDouble();  break;
            case 14: Debug_       = p.toDouble();  break;
            case 15: Omega_00_    = p.toDouble();  break;
            case 16: Omega_11_    = p.toDouble();  break;
            case 17: Omega_22_    = p.toDouble();  break;
            case 18: Omega_01_    = p.toDouble();  break;
            case 19: Omega_12_    = p.toDouble();  break;
            case 20: Omega_20_    = p.toDouble();  break;
            case 21: Epsid_00_    = p.toDouble();  break;
            case 22: Epsid_11_    = p.toDouble();  break;
            case 23: Epsid_22_    = p.toDouble();  break;
            case 24: Epsid_01_    = p.toDouble();  break;
            case 25: Epsid_12_    = p.toDouble();  break;
            case 26: Epsid_20_    = p.toDouble();  break;
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
        E0_         = em->E0_      ;
        Poisson0_   = em->Poisson0_;
        a1_         = em->a1_      ;
        a2_         = em->a2_      ;
        a3_         = em->a3_      ;
        a4_         = em->a4_      ;
        C0_         = em->C0_      ;
        C1_         = em->C1_      ;
        alpha_      = em->alpha_   ;
        Debug_      = em->Debug_   ;
        Omega_00_   = em->Omega_00_;
        Omega_11_   = em->Omega_11_;
        Omega_22_   = em->Omega_22_;
        Omega_01_   = em->Omega_01_;
        Omega_12_   = em->Omega_12_;
        Omega_20_   = em->Omega_20_;
        Epsid_00_   = em->Epsid_00_;
        Epsid_11_   = em->Epsid_11_;
        Epsid_22_   = em->Epsid_22_;
        Epsid_01_   = em->Epsid_01_;
        Epsid_12_   = em->Epsid_12_;
        Epsid_20_   = em->Epsid_20_;
        
    }

    void Modeldsid::initialize(UByte dim,State *s) {
        ConstitutiveModel::initialize(dim,s);
        // initialize mean pressure
        Omega[0] = Omega_00_;
        Omega[1] = Omega_11_;
        Omega[2] = Omega_22_;
        Omega[3] = Omega_01_;
        Omega[4] = Omega_12_;
        Omega[5] = Omega_20_;
        effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_);
        Epsid[0] = Epsid_00_;
        Epsid[1] = Epsid_11_;
        Epsid[2] = Epsid_22_;
        Epsid[3] = Epsid_01_;
        Epsid[4] = Epsid_12_;
        Epsid[5] = Epsid_20_;

        dstran[0] = s->stnE_.s11();
        dstran[1] = s->stnE_.s22();
        dstran[2] = s->stnE_.s33();
        dstran[3] = s->stnE_.s12();
        dstran[4] = s->stnE_.s23();
        dstran[5] = s->stnE_.s13();

        Stress[0] = s->stnS_.s11();
        Stress[1] = s->stnS_.s22();
        Stress[2] = s->stnS_.s33();
        Stress[3] = s->stnS_.s12();
        Stress[4] = s->stnS_.s13();
        Stress[5] = s->stnS_.s23();

        double E00 = 1./Matdom[0][0];
        double E11 = 1./Matdom[1][1];
        double E22 = 1./Matdom[2][2];
        double G01 = Matdom[3][3];
        double G12 = Matdom[4][4];
        double G20 = Matdom[5][5];
        double Emax;
        Emax = MAX(E00, MAX(E11,E22));
        Shear_ = MAX(G01, MAX(G12,G20));
        Bulk = Emax*Shear_/3./(3.*Shear_-Emax);
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
        int ntens = Omega.size();
        double zero = 0.;
        double deps=0.;
        double fd, fdt, XL;
        int iopt;
        const double tol = 1e-6, ITmax = 250, tiny = 1e-3;
        for (int i=0; i<ntens; i++) {deps+= dstran[i]*dstran[i];}
        deps=sqrt(deps);
        for (int i=0; i<ntens; i++) {
            dSig[i]=0.;
            for (int j=0; j<ntens; j++) {
                if (j<3) {
                    dSig += Matdom[i][j]*dstran[j];
                } else {
                    dSig += 2.* Matdom[i][j]*dstran[j];
            }
        }

        for (int i=0; i<ntens; i++) {
            Stress[i] += dSig[i];
        }
     
        damageFunction(fd, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_);
        if (fd > tol) {
            int Inc = 0;
            double fdt = fd;
            while(fdt>0 && ((fdt/fd)>tol && Inc < ITmax)) {
                iopt = 0;
                if (Inc == 0) iopt =0;
                cuttingPlaneMethod(Omega, Sigma, Matdom, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                                   C0_, C1_, alpha_, H0, Hp, dG_dY, df_dSig, temp, iopt);
                XL = fdt/(H0-HP);
                if (dstrain>tiny) XL/=1.5;
                for (int i=0; i<ntens; i++) {
                    Omega[i] += XL*dG_dY[i];
                    Epsid[i] += XL*df_dSig[i];
                    Stress[i] -= XL*temp[i];
                }
                ioptfd = 2;
                damageFunction(fdt, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
                Inc += 1;
            }
            if (Inc>=ITmax) throwout("DSID: no convergence");
        }
        ioptfd = 3;
        damageFunction(fd2, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
        
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


void effectiveStiffness(r2Tensor<double> &Matdom, const r1Tensor<double> &Omega, const double & E0_,
                        const double & Poisson0_, const double & a1_, const double & a2_,
                        const double & a3_, const double & a4_, const double & C0_,
                        const double & C1_, const double & alpha_){
     int ntens = Matdom.dim1();
     double zero = 0;
     double trOmega, b1, b2, coe1, coe2,
     r2Tensor<double> MatS(ntens,ntens,zero);
     Matdom = MatS;

       trOmega = Omega[0] + Omega[1] + Omega[2];
       b1 = (1. + Poisson0_)/E0_/2.;
       b2 = Poisson0_/E0_;

       coe1 = 2.;
       coe2 = 4.;

       MatS[0][0]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[0]+2.*a4_*trOmega;
       MatS[0][1]=-b2+2.*a1_*trOmega+a3_*(Omega[0]+Omega[1]);
       MatS[0][2]=-b2+2.*a1_*trOmega+a3_*(Omega[2]+Omega[0]);
       MatS[0][3]=coe1*(a2_*Omega[3]+a3_*Omega[3]);

       MatS[1][0]=MatS[0][1];
       MatS[1][1]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[1]+2.*a4_*trOmega;
       MatS[1][2]=-b2+2.*a1_*trOmega+a3_*(Omega[2]+Omega[1]);
       MatS[1][3]=coe1*(a2_*Omega[3]+a3_*Omega[3]);

       MatS[2][0]=MatS[0][2];
       MatS[2][1]=MatS[1][2];
       MatS[2][2]=2.*b1-b2+2.*a1_*trOmega+2.*(a2_+a3_)*Omega[2]+2.*a4_*trOmega;
       MatS[2][3]=coe1*(a3_*Omega[3]);

       MatS[3][0]=MatS[0][3];
       MatS[3][1]=MatS[1][3];      
       MatS[3][2]=MatS[2][3];       
       MatS[3][3]=coe2*(b1+0.5*a2_*(Omega[0]+Omega[1])+a4_*trOmega);
 
       
//       IF (NTENS.EQ.6] THEN
         MatS[0][4]=coe1*(a3_*Omega[4]);
         MatS[0][5]=coe1*(a2_*Omega[5]+a3_*Omega[5]);

         MatS[1][4]=coe1*(a2_*Omega[4]+a3_*Omega[4]);
         MatS[1][5]=coe1*(a3_*Omega[5]);

         MatS[2][4]=coe1*(a2_*Omega[4]+a3_*Omega[4]);
         MatS[2][5]=coe1*(a2_*Omega[5]+a3_*Omega[5]);

         MatS[3][4]=coe2*0.5*a2_*Omega[5];
         MatS[3][5]=coe2*0.5*a2_*Omega[4];

         MatS[4][0]=MatS[0][4];
         MatS[4][1]=MatS[1][4];      
         MatS[4][2]=MatS[2][4];       
         MatS[4][3]=MatS[3][4];
         MatS[4][4]=coe2*(b1+0.5*a2_*(Omega[2]+Omega[1])+a4_*trOmega);
         MatS[4][5]=coe2*0.5*a2_*Omega[3];

         MatS[5][0]=MatS[0][5];
         MatS[5][1]=MatS[1][5];       
         MatS[5][2]=MatS[2][5];       
         MatS[5][3]=MatS[3][5];
         MatS[5][4]=MatS[4][5];
         MatS[5][5]=coe2*(b1+0.5*a2_*(Omega[2]+Omega[0])+a4_*trOmega);


//       ENDIF

        gaussj(MatS,Matdom);
     
}

void damageFunction(double & fd, const r1Tensor<double> &Sigma, const r1Tensor<double> &Omega, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_, const int & ioptfd) {

     int ntens = Sigma.size();
     double zero = 0;
     double trSigma, trOmega, trSigSig, trY, SS;
     r1Tensor<double> SigSig(ntens,zero),P1Y(ntens,zero),sij(ntens,zero),e1(ntens,zero),yd1(ntens,zero);
     r2Tensor<double> P1(ntens,ntens,zero);
     for (int i=0; i<3; i++) e1[i]=1;
     
     trSigma = Sigma[1] + Sigma[2] + Sigma[3];
     trOmega = Omega[1] + Omega[2] + Omega[3];

     Aik_Bkj(Sigma,Sigma,SigSig);

     trSigSig = SigSig[1] + SigSig[2] + SigSig[3];

     for (int i=0; i<ntens; i++)
          yd1[i]= a1_*trSigma**2.*e1[i]+a2_*SigSig[i]+a3_*trSigma*Sigma[i]+a4_*trSigSig*e1[i];
     

     MATP1(P1,Sigma);

     for (int i=0; i<ntens; ++i) {
            P1Y[i] = 0.;
         for (int j=0; j<ntens; j++) {
             if (j<3) {
                P1Y[i] += P1[i][j]*yd1[j];
             } else {
                P1Y[i] += 2.*P1[i][j]*yd1[j];
             }  
         }
     }

     trY = P1Y[1] + P1Y[2] + P1Y[3];

     for (int i=0; i<ntens; ++i)
          sij[i] = P1Y[i]-1./3.*trY*e1[i];

     SS=0.;
     for (int i=0; i<ntens; ++i) {
         if (i<3) {
            SS += sij[i]*sij[i];
         } else {
            SS += 2.*sij[i]*sij[i];
         }  
     }

    /*
    fd = alpha_*trY-C0_-C1_*trOmega;
    if (fd>0.) {
        throwout('DSID:  Stresses are in tension');
        cout << "IOPT = " << IOPTfd << endl;
        cout << "Sigma = " << Sigma[1] << " " << Sigma[2] << " " << Sigma[3] << " " << Sigma[4] << endl;
    } 
    */
    //THE SIGN BEFORE alpha_ IS "+" DUE TO MECHANICAL CONVENTION
    fd = sqrt(0.5*SS)+alpha_*trY-C0_-C1_*trOmega;   

}

void matP1(r2Tensor<double> &P1, const r1Tensor<double> &Sigma) {
    int ntens = Sigma.size();
    int m=3;
    r1Tensor<double> s(m), anan1(ntens), anan2(ntens), anan3(ntens);
    r2Tensor<double> sigma1(m,m), sigma0(m,m);
    const double one=1.,two=2.,tol=1e-6;
    vectorToTensor(Sigma, sigma1, one);
    sigma0 = sigma1;
    Unsymmeig h(sigma0);
    for (int i=0; i<m; i++) {
        if (abs(h.wri[i].real())<tol) h.wri[i].real()=0.;
    }
    for (int i=0; i<m; i++) {
        if (h.wri[i].real()>=0){
            s[i] = 1.;
        } else {
            s[i] = -1.;
        }
    }
    for (int i=0; i<m; i++) {
        anan1[i] = h.zz[i][0]*h.zz[i][0];
        anan2[i] = h.zz[i][1]*h.zz[i][1];
        anan3[i] = h.zz[i][2]*h.zz[i][2];
    }

    anan1[3] = h.zz[0][0]*h.zz[1][0];
    anan2[3] = h.zz[0][1]*h.zz[1][1];
    anan3[3] = h.zz[0][2]*h.zz[1][2];

    anan1[4] = h.zz[1][0]*h.zz[2][0];
    anan2[4] = h.zz[1][1]*h.zz[2][1];
    anan3[4] = h.zz[1][2]*h.zz[2][2];

    anan1[5] = h.zz[2][0]*h.zz[0][0];
    anan2[5] = h.zz[2][1]*h.zz[0][1];
    anan3[5] = h.zz[2][2]*h.zz[0][2];

    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            P1[i][j] = s[0] * anan1[i]*anan1[j] +
                       s[1] * anan2[i]*anan2[j] +
                       s[2] * anan3[i]*anan3[j];
        }
    }
}

void matP2(r2Tensor<double> &P2, const r1Tensor<double> &Sigma);
    int ntens = Sigma.size();
    int m=3;i
    double res;
    r1Tensor<double> s(m), anan1(ntens), anan2(ntens), anan3(ntens);
    r2Tensor<double> sigma1(m,m), sigma0(m,m);
    const double one=1.,two=2.,tol=1e-6;
    vectorToTensor(Sigma, sigma1, one);
    sigma0 = sigma1;
    Unsymmeig h(sigma0);
    //for (int i=0; i<m; i++) {
    //    if (abs(h.wri[i].real())<tol) h.wri[i].real()=0.;
    //}
    for (int i=0; i<m; i++) {
        res = h.wri[i].real()- MIN(h.wri[0].real(),MIN(h.wri[1].real(),h.wri[2].real()));
        if (res>0.){
            s[i] = 1.;
        } else {
            s[i] = 0.;
        }
    }
    for (int i=0; i<m; i++) {
        anan1[i] = h.zz[i][0]*h.zz[i][0];
        anan2[i] = h.zz[i][1]*h.zz[i][1];
        anan3[i] = h.zz[i][2]*h.zz[i][2];
    }

    anan1[3] = h.zz[0][0]*h.zz[1][0];
    anan2[3] = h.zz[0][1]*h.zz[1][1];
    anan3[3] = h.zz[0][2]*h.zz[1][2];

    anan1[4] = h.zz[1][0]*h.zz[2][0];
    anan2[4] = h.zz[1][1]*h.zz[2][1];
    anan3[4] = h.zz[1][2]*h.zz[2][2];

    anan1[5] = h.zz[2][0]*h.zz[0][0];
    anan2[5] = h.zz[2][1]*h.zz[0][1];
    anan3[5] = h.zz[2][2]*h.zz[0][2];

    for (int i=0; i<m; i++) {
        for (int j=0; j<m; j++) {
            P2[i][j] = s[0] * anan1[i]*anan1[j] +
                       s[1] * anan2[i]*anan2[j] +
                       s[2] * anan3[i]*anan3[j];
        }
    }

void dY_dSigFunction(const r1Tensor<double> &Sigma, r2Tensor<double> &dY_dSig,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_) {
    int ntens = Sigma.size();
    double trSig;
    trSig = Sigma[0]+Sigma[1]+Sigma[2];
    dY_dSig[0][0]=2.*a1_*trSig+2.*a2_*Sigma[0]+a3_*(Sigma[0]+trSig)+2.*a4_*Sigma[0];
    dY_dSig[0][1]=2.*a1_*trSig+a3_*Sigma[0]+2.*a4_*Sigma[1];
    dY_dSig[0][2]=2.*a1_*trSig+a3_*Sigma[0]+2.*a4_*Sigma[2];
    dY_dSig[0][3]=a2_*Sigma[3]+2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[1][0]=2.*a1_*trSig+a3_*Sigma[1]+2.*a4_*Sigma[0];
    dY_dSig[1][1]=2.*a1_*trSig+2.*a2_*Sigma[1]+a3_*(Sigma[1]+trSig)+2.*a4_*Sigma[1];
    dY_dSig[1][2]=2.*a1_*trSig+a3_*Sigma[1]+2.*a4_*Sigma[2];
    dY_dSig[1][3]=a2_*Sigma[3]+2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[2][0]=2.*a1_*trSig+a3_*Sigma[2]+2.*a4_*Sigma[0];
    dY_dSig[2][1]=2.*a1_*trSig+a3_*Sigma[2]+2.*a4_*Sigma[1];
    dY_dSig[2][2]=2.*a1_*trSig+2.*a2_*Sigma[2]+a3_*(Sigma[2]+trSig)+2.*a4_*Sigma[2];
    dY_dSig[2][3]=2.*a4_*Sigma[3];
                                                                                   
    dY_dSig[3][0]=a2_*Sigma[3]+a3_*Sigma[3];
    dY_dSig[3][1]=a2_*Sigma[3]+a3_*Sigma[3];
    dY_dSig[3][2]=a3_*Sigma[3];
    dY_dSig[3][3]=0.5*a2_*(Sigma[0]+Sigma[1])+0.5*a3_*trSig;
                                                                                   
    dY_dSig[0][4]=2.*a4_*Sigma[4];
    dY_dSig[0][5]=a2_*Sigma[5]+2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[1][4]=a2_*Sigma[4]+2.*a4_*Sigma[4];
    dY_dSig[1][5]=2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[2][4]=a2_*Sigma[4]+2.*a4_*Sigma[4];
    dY_dSig[2][5]=a2_*Sigma[5]+2.*a4_*Sigma[5];
                                                                                   
    dY_dSig[3][4]=0.5*a2_*Sigma[5];
    dY_dSig[3][5]=0.5*a2_*Sigma[4];
                                                                                   
    dY_dSig[4][0]=a3_*Sigma[4];
    dY_dSig[4][1]=a2_*Sigma[4]+a3_*Sigma[4];
    dY_dSig[4][2]=a2_*Sigma[4]+a3_*Sigma[4];
    dY_dSig[4][3]=dY_dSig[3][4];
    dY_dSig[4][4]=0.5*a2_*(Sigma[1]+Sigma[2])+0.5*a3_*trSig;
    dY_dSig[4][5]=0.5*a2_*Sigma[3];
                                                                                   
    dY_dSig[5][0]=a2_*Sigma[5]+a3_*Sigma[5];
    dY_dSig[5][1]=a3_*Sigma[5];
    dY_dSig[5][2]=a2_*Sigma[5]+a3_*Sigma[5];
    dY_dSig[5][3]=dY_dSig[3][5];
    dY_dSig[5][4]=dY_dSig[4][5];
    dY_dSig[5][5]=0.5*a2_*(Sigma[2]+Sigma[0])+0.5*a3_*trSig;
}

// SUB FD_LAM
void cuttingPlaneMethod(const r1Tensor<double> &Omega, const r1Tensor<double> &Sigma, r2tensor<double> &Matdom,const double & E0_,
             const double & Poisson0_, const double & a1_, const double & a2_,
             const double & a3_, const double & a4_, const double & C0_,
             const double & C1_, const double & alpha_, double & H0,
             double & Hp, r1Tensor<double> &dG_dY, r1Tensor<double> &df_dSig,
             r1Tensor<double> &temp, const int & iopt){
    int ntens = Omega.size();
    double zero=0.;
    double trSigma, trOmega, trSigSig, P1yd1e, f1f1, f2f2;
    r1Tensor<double> e(ntens), yd1(ntens), SigSig(ntens), f1ij(ntens),f2ij(ntens);
    r1Tensor<double> P1yd1(ntens), f2p2(ntens), ep1(ntens),df_dY(ntens), df_dOmega(ntens);
    r1Tensor<double> f1p3(ntens),zeros(ntens,zero),temp1(ntens); 
    r2Tensor<double> P1(ntens,ntens), dY_dSig(ntens,ntens),eep1(ntens,ntens);
    r2Tensor<double> P3(ntens,ntens), P2(ntens,ntens), temp2(ntens,ntens);
    r3Tensor<double> dS_dOmega; 

    //DATA ONE,TWO,HALF / 1.0D0,2.0D0,0.5D0 /

    for (int i=0; i<ntens; i++) {
        if (i<3) {
            e[i]=1.;
        } else {
            e[i]=0.;
        }
    }

    trSigma = Sigma[0]+Sigma[1]+Sigma[2];
    trOmega = Omega[0]+Omega[1]+Omega[2];

    Aik_Bkj(Sigma,Sigma,SigSig);

    trSigSig =SigSig[0]+SigSig[1]+SigSig[2];

    for (int i=0; i<ntens; i++) {
        yd1[i]= a1_*trSigma**2.*e[i]+a2_*SigSig[i]+
              a3_*trSigma*Sigma[i]+a4_*trSigSig*e[i];
    }

    dY_dSigFunction(Sigma,dY_dSig,
            E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);
    matP1(P1,Sigma);

    matP2(P2,Sigma);
    for (int i=0;i<ntens;i++){
        P1yd1[i]=0.;
        for (int j=0;j<ntens;j++){
            if (j<3) {
                P1yd1[i]+=P1[i][j]*yd1[j];
            } else {
                P1yd1[i]+=2.*P1[i][j]*yd1[j];
            }
         }   
    }
    P1yd1e=P1yd1[0]+P1yd1[1]+P1yd1[2];

    for (int i=0; i<ntens; i++) {
        f2ij[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                f2ij[i]+=P2[i][j]*yd1[j];
            } else {
                f2ij[i]+=2.*P2[i][j]*yd1[j];
            }
         }
     }

    for (int i=0; i<ntens; i++){
        f2p2[i]=0.;
        for (int j=0;j<ntens;j++){
            if (j<3) {   
                f2p2[i]+=P2[j][i]*f2ij[j];
            } else {
                f2p2[i]+=2.*P2[j][i]*f2ij[j];
            }  
         }
     }
      
     f2f2=0.;
     for (int i=0; i<ntens; i++) {
         if (i<3) {
             f2f2+=f2ij[i]*f2ij[i];
         } else {
             f2f2+=2.*f2ij[i]*f2ij[i];
         }  
     }

    for (int i=0; i<ntens; i++) {
        ep1[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                ep1[i]+=P1[j][i]*e[j];
            } else {
                ep1[i]+=2.*P1[j][i]*e[j];
            }  
        }
    }
     
    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            eep1[i][j]=e[i]*ep1[j];
        }
    }

    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            P3[i][j]=P1[i][j]-1./3.*eep1[i][j];
        }
    }

    for (int i=0; i<ntens; i++) {
        f1ij[i]=P1yd1[i]-1./3.*P1yd1e*e[i];
        df_dOmega[i]=-C1_*e[i];
        if (f2f2 == 0.) {
            dG_dY[i]=0.;
        } else {  
            dG_dY[i]=f2p2[i]/sqrt(2.0*f2f2);
        }
    }

      
    for (int i=0; i<ntens; i++) {
        f1p3[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                f1p3[i]+=f1ij[j]*P3[j][i];
            } else {
                f1p3[i]+=2.*f1ij[j]*P3[j][i];
            }  
        }
    }
      
    f1f1=0.
    for (int i=0; i<ntens; i++) {
        if (i<3) {
            f1f1+=f1ij[i]*f1ij[i];
        } else {
            f1f1+=2.*f1ij[i]*f1ij[i];
        }  
    }

    for (int i=0; i<ntens; i++) {
        df_dY[i]=f1p3[i]/sqrt(2.*f1f1)+alpha_*ep1[i];
    }
     
    for (int i=0; i<ntens; i++) {
        df_dSig[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3) {
                df_dSig[i]+=df_dY[j]*dY_dSig[j][i];
            } else {
                df_dSig[i]+=2.*df_dY[j]*dY_dSig[j][i];
            }  
        }
    }

      
    HP=0.;
    for (int i=0; i<ntens; i++) {
        if(i<3){
            HP+=df_dOmega[i]*dG_dY[i];
        } else {
            HP+=2.*df_dOmega[i]*dG_dY[i];
        } 
    }

    H0=0.;

    Mat_dS_dOmega(dS_dOmega,E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);

    for (int i=0; i<ntens; i++) {
        for (int j=0; j<ntens; j++) {
            temp2[i][j]=0.;
            for (int k=0; k<ntens; k++) {
                if (k<3) {
                    temp2[i][j]=temp2[i][j]+Sigma[k]*dS_dOmega[k][i][j];
                } else {
                    temp2[i][j]=temp2[i][j]+2.*Sigma[k]*dS_dOmega[k][i][j];
                }
            }
        }
    }

    for (int i=0; i<ntens; i++) {
        temp1[i]=0.;
        for (int j=0; j<ntens; j++) {
            if (j<3){
                temp1[i]+=temp2[i][j]*dG_dY[j];
            } else {
                temp1[i]+=2.0*temp2[i][j]*dG_dY[j];
            }
        }
    }

    for (int i=0; i<ntens; i++) {
         temp1[i]+=df_dSig[i];
    }
     
    if(iopt.eq.1) {
      effectiveStiffness(Matdom,Omega,
               E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_,alpha_);
    }
    for (int i=0; i<ntens; i++) {
        temp[i]=0.;
        for (int j=0; j<ntens; j++) {
            if(j<3){
                temp[i]+=Matdom[i][j]*temp1[j];
            } else {
                temp[i]+=2.*Matdom[i][j]*temp1[j];
            } 
        }
    }
    for (int i=0; i<ntens; i++) {
        if(i<3){
            H0+=df_dSig[i]*temp[i];
        } else {
            H0+=2.*df_dSig[i]*temp[i];
        } 
    }
}

void Mat_dS_dOmega(r3Tensor<double> &dS_dO, const double & E0_,
                   const double & Poisson0_, const double & a1_, const double & a2_,
                   const double & a3_, const double & a4_, const double & C0_,
                   const double & C1_, const double & alpha_){

// ------- page 1 -------
      dS_dO[0][0][0]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;  
      dS_dO[0][1][0]=2.*a1_+a3_                 ;
      dS_dO[0][2][0]=2.*a1_+a3_                 ;
      dS_dO[0][3][0]=0.                         ;
                                                  
      dS_dO[1][0][0]=dS_dO[0][1][0]             ;
      dS_dO[1][1][0]=2.*a1_+2.*a4_              ;
      dS_dO[1][2][0]=2.*a1_                     ;
      dS_dO[1][3][0]=0.                         ;
                                                  
      dS_dO[2][0][0]=dS_dO[0][2][0]             ;
      dS_dO[2][1][0]=dS_dO[1][2][0]             ;
      dS_dO[2][2][0]=2.*a1_+2.*a4_              ;
      dS_dO[2][3][0]=0.                         ;
                                                  
      dS_dO[3][0][0]=dS_dO[0][3][0]             ;
      dS_dO[3][1][0]=dS_dO[1][3][0]             ;
      dS_dO[3][2][0]=dS_dO[2][3][0]             ;
      dS_dO[3][3][0]=0.5*a2_+a4_                ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][0]=0.                         ;
      dS_dO[0][5][0]=0.                         ;
      dS_dO[1][4][0]=0.                         ;
      dS_dO[1][5][0]=0.                         ;
      dS_dO[2][4][0]=0.                         ;
      dS_dO[2][5][0]=0.                         ;
      dS_dO[3][4][0]=0.                         ;
      dS_dO[3][5][0]=0.                         ;
      dS_dO[4][0][0]=dS_dO[0][4][0]             ;
      dS_dO[4][1][0]=dS_dO[1][4][0]             ;
      dS_dO[4][2][0]=dS_dO[2][4][0]             ;
      dS_dO[4][3][0]=dS_dO[3][4][0]             ;
      dS_dO[4][4][0]=a4_                        ;
      dS_dO[4][5][0]=0.                         ;
                                                  
      dS_dO[5][0][0]=dS_dO[0][5][0]             ;
      dS_dO[5][1][0]=dS_dO[1][5][0]             ;
      dS_dO[5][2][0]=dS_dO[2][5][0]             ;
      dS_dO[5][3][0]=dS_dO[3][5][0]             ;
      dS_dO[5][4][0]=dS_dO[4][5][0]             ;
      dS_dO[5][5][0]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 2 -------                       --
      dS_dO[0][0][1]=2.*a1_+2.*a4_              ;
      dS_dO[0][1][1]=2.*a1_+a3_                 ;
      dS_dO[0][2][1]=2.*a1_                     ;
      dS_dO[0][3][1]=0.                         ;
                                                  
      dS_dO[1][0][1]=dS_dO[0][1][1]             ;
      dS_dO[1][1][1]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;
      dS_dO[1][2][1]=2.*a1_+a3_                 ;
      dS_dO[1][3][1]=0.                         ;
                                                  
      dS_dO[2][0][1]=dS_dO[0][2][1]             ;
      dS_dO[2][1][1]=dS_dO[1][2][1]             ;
      dS_dO[2][2][1]=2.*a1_+2.*a4_              ;
      dS_dO[2][3][1]=0.                         ;
                                                  
      dS_dO[3][0][1]=dS_dO[0][3][1]             ;
      dS_dO[3][1][1]=dS_dO[1][3][1]             ;
      dS_dO[3][2][1]=dS_dO[2][3][1]             ;
      dS_dO[3][3][1]=0.5*a2_+a4_                ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][1]=0.                         ;
      dS_dO[0][5][1]=0.                         ;
      dS_dO[1][4][1]=0.                         ;
      dS_dO[1][5][1]=0.                         ;
      dS_dO[2][4][1]=0.                         ;
      dS_dO[2][5][1]=0.                         ;
      dS_dO[3][4][1]=0.                         ;
      dS_dO[3][5][1]=0.                         ;
      dS_dO[4][0][1]=dS_dO[0][4][1]             ;
      dS_dO[4][1][1]=dS_dO[1][4][1]             ;
      dS_dO[4][2][1]=dS_dO[2][4][1]             ;
      dS_dO[4][3][1]=dS_dO[3][4][1]             ;
      dS_dO[4][4][1]=0.5*a2_+a4_                ;
      dS_dO[4][5][1]=0.                         ;
                                                  
      dS_dO[5][0][1]=dS_dO[0][5][1]             ;
      dS_dO[5][1][1]=dS_dO[1][5][1]             ;
      dS_dO[5][2][1]=dS_dO[2][5][1]             ;
      dS_dO[5][3][1]=dS_dO[3][5][1]             ;
      dS_dO[5][4][1]=dS_dO[4][5][1]             ;
      dS_dO[5][5][1]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 3 -------                       --
      dS_dO[0][0][2]=2.*a1_+2.*a4_              ;
      dS_dO[0][1][2]=2.*a1_                     ;
      dS_dO[0][2][2]=2.*a1_+a3_                 ;
      dS_dO[0][3][2]=0.                         ;
                                                  
      dS_dO[1][0][2]=dS_dO[0][1][2]             ;
      dS_dO[1][1][2]=2.*a1_+2.*a4_              ;
      dS_dO[1][2][2]=2.*a1_+a3_                 ;
      dS_dO[1][3][2]=0.                         ;
                                                  
      dS_dO[2][0][2]=dS_dO[0][2][2]             ;
      dS_dO[2][1][2]=dS_dO[1][2][2]             ;
      dS_dO[2][2][2]=2.*a1_+2.*a2_+2.*a3_+2.*a4_;
      dS_dO[2][3][2]=0.                         ;
                                                  
      dS_dO[3][0][2]=dS_dO[0][3][2]             ;
      dS_dO[3][1][2]=dS_dO[1][3][2]             ;
      dS_dO[3][2][2]=dS_dO[2][3][2]             ;
      dS_dO[3][3][2]=a4_                        ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][2]=0.                         ;
      dS_dO[0][5][2]=0.                         ;
      dS_dO[1][4][2]=0.                         ;
      dS_dO[1][5][2]=0.                         ;
      dS_dO[2][4][2]=0.                         ;
      dS_dO[2][5][2]=0.                         ;
      dS_dO[3][4][2]=0.                         ;
      dS_dO[3][5][2]=0.                         ;
      dS_dO[4][0][2]=dS_dO[0][4][2]             ;
      dS_dO[4][1][2]=dS_dO[1][4][2]             ;
      dS_dO[4][2][2]=dS_dO[2][4][2]             ;
      dS_dO[4][3][2]=dS_dO[3][4][2]             ;
      dS_dO[4][4][2]=0.5*a2_+a4_                ;
      dS_dO[4][5][2]=0.                         ;
                                                  
      dS_dO[5][0][2]=dS_dO[0][5][2]             ;
      dS_dO[5][1][2]=dS_dO[1][5][2]             ;
      dS_dO[5][2][2]=dS_dO[2][5][2]             ;
      dS_dO[5][3][2]=dS_dO[3][5][2]             ;
      dS_dO[5][4][2]=dS_dO[4][5][2]             ;
      dS_dO[5][5][2]=0.5*a2_+a4_                ;
//      ENDIF                                     
                                                  
// ------- page 4 -------                       --
      dS_dO[0][0][3]=0.                         ;
      dS_dO[0][1][3]=0.                         ;
      dS_dO[0][2][3]=0.                         ;
      dS_dO[0][3][3]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[1][0][3]=dS_dO[0][1][3]             ;
      dS_dO[1][1][3]=0.                         ;
      dS_dO[1][2][3]=0.                         ;
      dS_dO[1][3][3]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[2][0][3]=dS_dO[0][2][3]             ;
      dS_dO[2][1][3]=dS_dO[1][2][3]             ;
      dS_dO[2][2][3]=0.                         ;
      dS_dO[2][3][3]=0.5*a3_                    ;
                                                  
      dS_dO[3][0][3]=dS_dO[0][3][3]             ;
      dS_dO[3][1][3]=dS_dO[1][3][3]             ;
      dS_dO[3][2][3]=dS_dO[2][3][3]             ;
      dS_dO[3][3][3]=0.                         ;
                                                  
//      IF[NTENS.EQ.6]THEN                        
      dS_dO[0][4][3]=0.                         ;
      dS_dO[0][5][3]=0.                         ;
      dS_dO[1][4][3]=0.                         ;
      dS_dO[1][5][3]=0.                         ;
      dS_dO[2][4][3]=0.                         ;
      dS_dO[2][5][3]=0.                         ;
      dS_dO[3][4][3]=0.                         ;
      dS_dO[3][5][3]=0.                         ;
      dS_dO[4][0][3]=dS_dO[0][4][3]             ;
      dS_dO[4][1][3]=dS_dO[1][4][3]             ;
      dS_dO[4][2][3]=dS_dO[2][4][3]             ;
      dS_dO[4][3][3]=dS_dO[3][4][3]             ;
      dS_dO[4][4][3]=0.                         ;
      dS_dO[4][5][3]=0.25*a2_                   ;
                                                  
      dS_dO[5][0][3]=dS_dO[0][5][3]             ;
      dS_dO[5][1][3]=dS_dO[1][5][3]             ;
      dS_dO[5][2][3]=dS_dO[2][5][3]             ;
      dS_dO[5][3][3]=dS_dO[3][5][3]             ;
      dS_dO[5][4][3]=dS_dO[4][5][3]             ;
      dS_dO[5][5][3]=0.                         ;
//      ENDIF                                     
                                                  
//      IF[NTENS.EQ.6]THEN                        
// ------- page 5 -------                       --
      dS_dO[0][0][4]=0.                         ;
      dS_dO[0][1][4]=0.                         ;
      dS_dO[0][2][4]=0.                         ;
      dS_dO[0][3][4]=0.                         ;
      dS_dO[0][4][4]=0.5*a3_                    ;
      dS_dO[0][5][4]=0.                         ;
                                                  
      dS_dO[1][0][4]=dS_dO[0][1][4]             ;
      dS_dO[1][1][4]=0.                         ;
      dS_dO[1][2][4]=0.                         ;
      dS_dO[1][3][4]=0.                         ;
      dS_dO[1][4][4]=0.5*a2_+0.5*a3_            ;
      dS_dO[1][5][4]=0.                         ;
                                                  
      dS_dO[2][0][4]=dS_dO[0][2][4]             ;
      dS_dO[2][1][4]=dS_dO[1][2][4]             ;
      dS_dO[2][2][4]=0.                         ;
      dS_dO[2][3][4]=0.                         ;
      dS_dO[2][4][4]=0.5*a2_+0.5*a3_            ;
      dS_dO[2][5][4]=0.                         ;
                                                  
      dS_dO[3][0][4]=dS_dO[0][3][4]             ;
      dS_dO[3][1][4]=dS_dO[1][3][4]             ;
      dS_dO[3][2][4]=dS_dO[2][3][4]             ;
      dS_dO[3][3][4]=0.                         ;
      dS_dO[3][4][4]=0.                         ;
      dS_dO[3][5][4]=0.25*a2_                   ;
                                                  
      dS_dO[4][0][4]=dS_dO[0][4][4]             ;
      dS_dO[4][1][4]=dS_dO[1][4][4]             ;
      dS_dO[4][2][4]=dS_dO[2][4][4]             ;
      dS_dO[4][3][4]=dS_dO[3][4][4]             ;
      dS_dO[4][4][4]=0.                         ;
      dS_dO[4][5][4]=0.                         ;
                                                  
      dS_dO[5][0][4]=dS_dO[0][5][4]             ;
      dS_dO[5][1][4]=dS_dO[1][5][4]             ;
      dS_dO[5][2][4]=dS_dO[2][5][4]             ;
      dS_dO[5][3][4]=dS_dO[3][5][4]             ;
      dS_dO[5][4][4]=dS_dO[4][5][4]             ;
      dS_dO[5][5][4]=0.                         ;
                                                  
// ------- page 6 -------                       --
      dS_dO[0][0][5]=0.                         ;
      dS_dO[0][1][5]=0.                         ;
      dS_dO[0][2][5]=0.                         ;
      dS_dO[0][3][5]=0.                         ;
      dS_dO[0][4][5]=0.                         ;
      dS_dO[0][5][5]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[1][0][5]=dS_dO[0][1][5]             ;
      dS_dO[1][1][5]=0.                         ;
      dS_dO[1][2][5]=0.                         ;
      dS_dO[1][3][5]=0.                         ;
      dS_dO[1][4][5]=0.                         ;
      dS_dO[1][5][5]=0.5*a3_                    ;
                                                  
      dS_dO[2][0][5]=dS_dO[0][2][5]             ;
      dS_dO[2][1][5]=dS_dO[1][2][5]             ;
      dS_dO[2][2][5]=0.                         ;
      dS_dO[2][3][5]=0.                         ;
      dS_dO[2][4][5]=0.                         ;
      dS_dO[2][5][5]=0.5*a2_+0.5*a3_            ;
                                                  
      dS_dO[3][0][5]=dS_dO[0][3][5]             ;
      dS_dO[3][1][5]=dS_dO[1][3][5]             ;
      dS_dO[3][2][5]=dS_dO[2][3][5]             ;
      dS_dO[3][3][5]=0.                         ;
      dS_dO[3][4][5]=0.25*a2_                   ;
      dS_dO[3][5][5]=0.                         ;
                                                  
      dS_dO[4][0][5]=dS_dO[0][4][5]             ;
      dS_dO[4][1][5]=dS_dO[1][4][5]             ;
      dS_dO[4][2][5]=dS_dO[2][4][5]             ;
      dS_dO[4][3][5]=dS_dO[3][4][5]             ;
      dS_dO[4][4][5]=0.                         ;
      dS_dO[4][5][5]=0.                         ;
                                                  
      dS_dO[5][0][5]=dS_dO[0][5][5]             ;
      dS_dO[5][1][5]=dS_dO[1][5][5]             ;
      dS_dO[5][2][5]=dS_dO[2][5][5]             ;
      dS_dO[5][3][5]=dS_dO[3][5][5]             ;
      dS_dO[5][4][5]=dS_dO[4][5][5]             ;
      dS_dO[5][5][5]=0.                         ;
//      ENDIF

}


// EOF

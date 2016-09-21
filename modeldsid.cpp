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

    Modeldsid::Modeldsid(): Bulk_(0.0),Shear_(0.0),E0_(0.0),Poisson0_(0.0),a1_(0.0),
               a2_(0.0),a3_(0.0),a4_(0.0),C0_(0.0),C1_(0.0),
               alpha_(0.0),Debug_(0.0),Omega_00_(0.0),Omega_11_(0.0),Omega_22_(0.0),
               Omega_01_(0.0),Omega_12_(0.0),Omega_20_(0.0),Epsid_00_(0.0),Epsid_11_(0.0),
               Epsid_22_(0.0),Epsid_01_(0.0),Epsid_12_(0.0),Epsid_20_(0.0) {   }

    String Modeldsid::getProperties(void) const {
        return L"bulk,shear,E0,poisson0,a1,"
               L"a2,a3,a4,C0,C1,"
               L"alpha,debug,omega00,omega11,omega22,"
               L"omega01,omega12,omega20,epsid00,epsid11,"
               L"epsid22,epsid01,epsid12,epsid20";

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
            case 3:  return(E0_)      ;
            case 4:  return(Poisson0_);
            case 5:  return(a1_)      ;
            case 6:  return(a2_)      ;
            case 7:  return(a3_)      ;
            case 8:  return(a4_)      ;
            case 9:  return(C0_)      ;
            case 10: return(C1_)      ;
            case 11: return(alpha_)   ;
            case 12: return(Debug_)   ;
            case 13: return(Omega_00_);
            case 14: return(Omega_11_);
            case 15: return(Omega_22_);
            case 16: return(Omega_01_);
            case 17: return(Omega_12_);
            case 18: return(Omega_20_);
            case 19: return(Epsid_00_);
            case 20: return(Epsid_11_);
            case 21: return(Epsid_22_);
            case 22: return(Epsid_01_);
            case 23: return(Epsid_12_);
            case 24: return(Epsid_20_);
        }
        return(0.0);
    }

    void Modeldsid::setProperty(UInt ul,const Variant &p,UInt restoreVersion) {
        ConstitutiveModel::setProperty(ul,p,restoreVersion);
        switch (ul) {
            case  1: Bulk_        = p.toDouble();  break;
            case  2: Shear_       = p.toDouble();  break;
            case  3: E0_          = p.toDouble();  break;
            case  4: Poisson0_    = p.toDouble();  break;
            case  5: a1_          = p.toDouble();  break;
            case  6: a2_          = p.toDouble();  break;
            case  7: a3_          = p.toDouble();  break;
            case  8: a4_          = p.toDouble();  break;
            case  9: C0_          = p.toDouble();  break;
            case 10: C1_          = p.toDouble();  break;
            case 11: alpha_       = p.toDouble();  break;
            case 12: Debug_       = p.toDouble();  break;
            case 13: Omega_00_    = p.toDouble();  break;
            case 14: Omega_11_    = p.toDouble();  break;
            case 15: Omega_22_    = p.toDouble();  break;
            case 16: Omega_01_    = p.toDouble();  break;
            case 17: Omega_12_    = p.toDouble();  break;
            case 18: Omega_20_    = p.toDouble();  break;
            case 19: Epsid_00_    = p.toDouble();  break;
            case 20: Epsid_11_    = p.toDouble();  break;
            case 21: Epsid_22_    = p.toDouble();  break;
            case 22: Epsid_01_    = p.toDouble();  break;
            case 23: Epsid_12_    = p.toDouble();  break;
            case 24: Epsid_20_    = p.toDouble();  break;
        }
    }

    void Modeldsid::copy(const ConstitutiveModel *m) {
        ConstitutiveModel::copy(m);
        const Modeldsid *em = dynamic_cast<const Modeldsid *>(m);
        if (!em) throw std::runtime_error("Internal error: constitutive model dynamic cast failed.");
        Bulk_       = em->Bulk_    ;
        Shear_      = em->Shear_   ;
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
        Omega[0] = Omega_00_;
        Omega[1] = Omega_11_;
        Omega[2] = Omega_22_;
        Omega[3] = Omega_01_;
        Omega[4] = Omega_12_;
        Omega[5] = Omega_20_;
        effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_);

        double E00 = 1./Matdom[0][0];
        double E11 = 1./Matdom[1][1];
        double E22 = 1./Matdom[2][2];
        double G01 = Matdom[3][3];
        double G12 = Matdom[4][4];
        double G20 = Matdom[5][5];
        double Emax;
        Emax = MAX(E00, MAX(E11,E22));
        Shear_ = MAX(G01, MAX(G12,G20));
        Bulk_ = Emax*Shear_/3./(3.*Shear_-Emax);
    }

    static const UInt dOmega_00 =  0;
    static const UInt dOmega_11 =  1;
    static const UInt dOmega_22 =  2;
    static const UInt dOmega_01 =  3;
    static const UInt dOmega_12 =  4;
    static const UInt dOmega_20 =  5;
    static const UInt dEpsid_00 =  6;
    static const UInt dEpsid_11 =  7;
    static const UInt dEpsid_22 =  8;
    static const UInt dEpsid_01 =  9;
    static const UInt dEpsid_12 = 10;
    static const UInt dEpsid_20 = 11;
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
            s->working_[dOmega_00] = 0.0;
            s->working_[dOmega_11] = 0.0;
            s->working_[dOmega_22] = 0.0;
            s->working_[dOmega_01] = 0.0;
            s->working_[dOmega_12] = 0.0;
            s->working_[dOmega_20] = 0.0;
            s->working_[dEpsid_00] = 0.0;
            s->working_[dEpsid_11] = 0.0;
            s->working_[dEpsid_22] = 0.0;
            s->working_[dEpsid_01] = 0.0;
            s->working_[dEpsid_12] = 0.0;
            s->working_[dEpsid_20] = 0.0;
        }

        /* --- trial elastic stresses --- */
        int ntens = Omega.size();
        double zero = 0.;
        double deps = 0.;
        double fd, fdt, fd2, XL;
        int iopt, ioptfd;
        const double tol = 1e-6, ITmax = 250, tiny = 1e-3;
        r1Tensor<double> Omega0(ntens), Epsid0(ntens);
        Omega[0] = Omega_00_;
        Omega[1] = Omega_11_;
        Omega[2] = Omega_22_;
        Omega[3] = Omega_01_;
        Omega[4] = Omega_12_;
        Omega[5] = Omega_20_;
        Omega0 = Omega;
        effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_);
        Epsid[0] = Epsid_00_;
        Epsid[1] = Epsid_11_;
        Epsid[2] = Epsid_22_;
        Epsid[3] = Epsid_01_;
        Epsid[4] = Epsid_12_;
        Epsid[5] = Epsid_20_;
        Epsid0 = Epsid;

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
        Stress[4] = s->stnS_.s23();
        Stress[5] = s->stnS_.s13();
        for (int i=0; i<ntens; i++) {deps += dstran[i]*dstran[i];}
        deps = sqrt(deps);
        for (int i=0; i<ntens; i++) {
            dSig[i] = 0.;
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

        s->stnS_.rs11() = Stress[0];
        s->stnS_.rs22() = Stress[1];
        s->stnS_.rs33() = Stress[2];
        s->stnS_.rs12() = Stress[3];
        s->stnS_.rs23() = Stress[4];
        s->stnS_.rs13() = Stress[5];
        // default settings, altered below if model is found to be failing
        s->viscous_ = true;  // allow viscous strains

        if (!canFail()) return;

        ioptfd = 1;
        damageFunction(fd, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
        if (fd > tol) {
            s->iworking_[Pind] = 1;
            s->state_ |= mShearNow;
            int Inc = 0;
            double fdt = fd;
            while(fdt>0 && ((fdt/fd)>tol && Inc < ITmax)) {
                iopt = 0;
                if (Inc == 0) iopt =0;
                cuttingPlaneMethod(Omega, Stress, Matdom, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                                   C0_, C1_, alpha_, H0, Hp, dG_dY, df_dSig, temp, iopt);
                XL = fdt/(H0-HP);
                if (deps>tiny) XL/=1.5;
                for (int i=0; i<ntens; i++) {
                    Omega[i] += XL*dG_dY[i];
                    Epsid[i] += XL*df_dSig[i];
                    Stress[i] -= XL*temp[i];
                }
                ioptfd = 2;
                damageFunction(fdt, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                       C0_, C1_, alpha_, ioptfd);
                ++Inc;
            }
            if (Inc>=ITmax) throw std::runtime_error("DSID: no convergence");
            ioptfd = 3;
            damageFunction(fd2, Stress, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                           C0_, C1_, alpha_, ioptfd);
            fd = fd2; 
            s->stnS_.rs11() = Stress[0];
            s->stnS_.rs22() = Stress[1];
            s->stnS_.rs33() = Stress[2];
            s->stnS_.rs12() = Stress[3];
            s->stnS_.rs23() = Stress[4];
            s->stnS_.rs13() = Stress[5];
        }



        /* update zone stresses and strains --- */
        Double dVol = s->getSubZoneVolume();
        //std::cout << " dPNew = " << dPNew <<std::endl;
        s->working_[dOmega_00] += (Omega[0] - Omega0[0]) * dVol;    
        s->working_[dOmega_11] += (Omega[1] - Omega0[1]) * dVol;
        s->working_[dOmega_22] += (Omega[2] - Omega0[2]) * dVol;
        s->working_[dOmega_01] += (Omega[3] - Omega0[3]) * dVol;
        s->working_[dOmega_12] += (Omega[4] - Omega0[4]) * dVol;
        s->working_[dOmega_20] += (Omega[5] - Omega0[5]) * dVol;  
        s->working_[dEpsid_00] += (Epsid[0] - Epsid0[0]) * dVol;    
        s->working_[dEpsid_11] += (Epsid[1] - Epsid0[1]) * dVol;
        s->working_[dEpsid_22] += (Epsid[2] - Epsid0[2]) * dVol;
        s->working_[dEpsid_01] += (Epsid[3] - Epsid0[3]) * dVol;    
        s->working_[dEpsid_12] += (Epsid[4] - Epsid0[4]) * dVol;
        s->working_[dEpsid_20] += (Epsid[5] - Epsid0[5]) * dVol;

        /* --- the last zone has been processed, update parameters: --- */
        if (s->sub_zone_==s->total_sub_zones_-1) {
            dVol = 1.0 / s->getZoneVolume();

            if (s->overlay_==2) dVol *= 0.5;
            s->working_[dOmega_00] *= dVol;
            s->working_[dOmega_11] *= dVol;
            s->working_[dOmega_22] *= dVol;
            s->working_[dOmega_01] *= dVol;
            s->working_[dOmega_12] *= dVol;
            s->working_[dOmega_20] *= dVol;
            s->working_[dEpsid_00] *= dVol;    
            s->working_[dEpsid_11] *= dVol;
            s->working_[dEpsid_22] *= dVol;
            s->working_[dEpsid_01] *= dVol;    
            s->working_[dEpsid_12] *= dVol;
            s->working_[dEpsid_20] *= dVol;

            Omega_00_ += s->working_[dOmega_00];
            Omega_11_ += s->working_[dOmega_11];
            Omega_22_ += s->working_[dOmega_22];
            Omega_01_ += s->working_[dOmega_01];
            Omega_12_ += s->working_[dOmega_12];
            Omega_20_ += s->working_[dOmega_20];
            Epsid_00_ += s->working_[dEpsid_00];
            Epsid_11_ += s->working_[dEpsid_11];
            Epsid_22_ += s->working_[dEpsid_22];
            Epsid_01_ += s->working_[dEpsid_01];
            Epsid_12_ += s->working_[dEpsid_12];
            Epsid_20_ += s->working_[dEpsid_20];

            Omega[0] = Omega_00_;
            Omega[1] = Omega_11_;
            Omega[2] = Omega_22_;
            Omega[3] = Omega_01_;
            Omega[4] = Omega_12_;
            Omega[5] = Omega_20_;
            effectiveStiffness(Matdom, Omega, E0_, Poisson0_, a1_, a2_, a3_, a4_,
                               C0_, C1_, alpha_);
            
            double E00 = 1./Matdom[0][0];
            double E11 = 1./Matdom[1][1];
            double E22 = 1./Matdom[2][2];
            double G01 = Matdom[3][3];
            double G12 = Matdom[4][4];
            double G20 = Matdom[5][5];
            double Emax;
            Emax = MAX(E00, MAX(E11,E22));
            Shear_ = MAX(G01, MAX(G12,G20));
            Bulk_ = Emax*Shear_/3./(3.*Shear_-Emax);
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

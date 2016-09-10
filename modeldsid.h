// modeldsid.cpp
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/07  Hao Xu   First release 

#ifndef MODELDSID_H_
#define MODELDSID_H_

#pragma once

#include "C:\Program Files\Itasca\Flac3d500\pluginfiles\models\src\conmodel.h"

namespace models {
  class Modeldsid : public ConstitutiveModel  {
    public:
      // Creators
      Modeldsid();
      virtual String        getName() const { return L"dsid"; }
      virtual String        getFullName() const { return L"dsid"; }
      virtual UInt          getMinorVersion() const;
      virtual String        getProperties() const;
      virtual String        getStates() const;
      virtual Variant       getProperty(UInt index) const;
      virtual void          setProperty(UInt index, const Variant &p,UInt restoreVersion=0);

      //virtual Modeldsid *clone() const { return(NEW("Modeldsid") Modeldsid()); }
      virtual Modeldsid    *clone() const { return new Modeldsid(); }
      virtual Double        getConfinedModulus() const { return(BulkB_ + Shear_*4.0/3.0); }
      virtual Double        getShearModulus() const { return(Shear_); }
      virtual Double        getBulkModulus() const { return(Bulk_); }
      virtual void          copy(const ConstitutiveModel *mod);
      virtual void          initialize(UByte dim,State *s);
      virtual void          run(UByte dim,State *s); 

      // Optional
      virtual bool          supportsHystereticDamping() const { return(false); }
    private:
      Double Bulk_,BulkB_,Shear_,Poisson_,E0_,Poisson0_,a1_,a2_,a3_,a4_;
      Double C0_,C1_,alpha_,Debug_,Omega_11_;
      Double Omega_22_,Omega_33_,Omega_12_,Omega_23_,Omega_31_;
      Double Epsid_11_,Epsid_22_,Epsid_33_,Epsid_12_,Epsid_23_;
      Double Epsid_31_;
      Double Matdom[6][6],Omega[6],Epsid[6];
  };
}
//EOF

#endif 

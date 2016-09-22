// modeldsid.cpp
//   A Continuum Damage Model used in FLAC3D.
// History:
// 2016/09/07  Hao Xu   First release 

#ifndef MODELDSID_H_
#define MODELDSID_H_

#pragma once

#include "C:\Program Files\Itasca\Flac3d500\pluginfiles\models\src\conmodel.h"
#include "..\mathLib\arithmetic.h"
#include "..\mathLib\r1Tensor.h"
#include "..\mathLib\r2Tensor.h"
#include "..\mathLib\r3Tensor.h"
#include "..\mathLib\errInfo.h"
#include "..\mathLib\gaussj.h"
#include "..\mathLib\eigen.h"

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
      virtual Double        getConfinedModulus() const { return(Bulk_ + Shear_*4.0/3.0); }
      virtual Double        getShearModulus() const { return(Shear_); }
      virtual Double        getBulkModulus() const { return(Bulk_); }
      virtual void          copy(const ConstitutiveModel *mod);
      virtual void          initialize(UByte dim,State *s);
      virtual void          run(UByte dim,State *s); 

      // Optional
      virtual bool          supportsHystereticDamping() const { return(false); }
    private:
      Double Bulk_,Shear_,E0_,Poisson0_,a1_,a2_,a3_,a4_,C0_,C1_;
      Double alpha_,Debug_,Omega_00_,Omega_11_,Omega_22_;
      Double Omega_01_,Omega_12_,Omega_20_,Epsid_00_,Epsid_11_;
      Double Epsid_22_,Epsid_01_,Epsid_12_,Epsid_20_;

      r2Tensor<double> Matdom;
      r1Tensor<double> Omega,Epsid,dstran,Stress,dSig;
  };
}
//EOF

#endif 

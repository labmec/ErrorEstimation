//
//  TPZMRCMLagrangeMultiplier.cpp
//  PZ
//
//  Created by Philippe Devloo on 5/2/14.
//
//

#include "TPZMRCMLagrangeMultiplier.h"
#include "pzaxestools.h"



/** @brief Unique identifier for serialization purposes */
int TPZMRCMLagrangeMultiplier::ClassId() const{
    return Hash("TPZMRCMLagrangeMultiplier") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

/** @brief Saves the element data to a stream */
void TPZMRCMLagrangeMultiplier::Write(TPZStream &buf, int withclassid) const
{
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
    buf.Write(&fNStateVariables);
    buf.Write(&fDimension);
    buf.Write(&fMultiplier);
    buf.Write(&fBeta);
}

/** @brief Reads the element data from a stream */
void TPZMRCMLagrangeMultiplier::Read(TPZStream &buf, void *context)
{
    TPZDiscontinuousGalerkin::Read(buf, context);
    buf.Read(&fNStateVariables);
    buf.Write(&fDimension);
    buf.Write(&fMultiplier);
    buf.Write(&fBeta);

}

//Contribution of skeletal elements.
void TPZMRCMLagrangeMultiplier::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int nmesh = datavec.size();
    if (nmesh!=2) DebugStop();

    TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phiP = datavec[1].phi;
    int phrq = phiQ.Rows();
    int phrp = phiP.Rows();
    
//------- Block of matrix B ------
    int iq, jp;
	for(iq = 0; iq<phrq; iq++) {
		for(jp=0; jp<phrp; jp++) {
            ek(iq, phrq+jp) += fMultiplier*weight*phiQ(iq,0)*phiP(jp,0);
		}
	}
    
    
//------- Block of matrix B^T ------
    int ip, jq;
	for(ip=0; ip<phrp; ip++) {
		for(jq=0; jq<phrq; jq++) {
			ek(ip + phrq,jq) += fMultiplier*weight*phiP(ip,0)*phiQ(jq,0);
		}
	}
}

/**
 * @brief Computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since June 5, 2012
 */
void TPZMRCMLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(dataleft[0].phi.Rows() == 0 || dataright[0].phi.Rows() == 0
       || dataleft[1].phi.Rows() == 0)
    {
        DebugStop();
    }
    
    TPZFMatrix<REAL> &phiFluxL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiFluxR = dataright[0].phi;
    TPZFMatrix<REAL> &phiPressL = dataleft[1].phi;
    REAL H = data.HSize;
    
    int nfluxl = phiFluxL.Rows();
    int npressl = phiPressL.Rows();
    int nfluxr = phiFluxR.Rows();
    static int count  = 0;

    if((nfluxl+npressl+nfluxr)*fNStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nfluxl <<
        " nrowl " << npressl << " may give wrong result " << std::endl;
        count++;
    }

    int firstblock = phiFluxL.Rows()*fNStateVariables;
    int secondblock = (phiFluxL.Rows()+phiPressL.Rows())*fNStateVariables;
    
    int fluxl_block = phiFluxL.Rows()*fNStateVariables;
    // 3) phi_press_left, phi_flux_right
    if(fBeta < 1.e7)
    {
        for(int il=0; il<npressl; il++) {
            for(int jr=0; jr<nfluxr; jr++) {
                for (int ist=0; ist<fNStateVariables; ist++) {
                    ek(firstblock+fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight *  (phiPressL(il) * phiFluxR(jr));
                    ek(fNStateVariables*jr+ist+secondblock,firstblock+fNStateVariables*il+ist) += weight *  (phiPressL(il) * phiFluxR(jr));
                }
            }
        }
    }
    else
    {
        for(int il=0; il<npressl; il++) {
            for(int jl=0; jl<npressl; jl++)
            {
                ek(firstblock+il,firstblock+jl) += weight *  (phiPressL(il) * phiPressL(jl));

            }
        }
    }
    
    if(fNStateVariables != 1) DebugStop();
    // penalty contribution
    STATE diag = fBeta;
    if(fBeta == 0.)
    {
        diag = 1.;
    }
    for(int il=0; il<nfluxl; il++) {
        for(int jl=0; jl<nfluxl; jl++) {
            ek(il,jl) += weight * diag * H * phiFluxL(il) * phiFluxL(jl);
        }
        for(int ir=0; ir<nfluxr; ir++) {
            ek(ir+secondblock,il) -= weight * fBeta * H * fMultiplier * phiFluxL(il) * phiFluxR(ir);
            ek(il,ir+secondblock) -= weight * fBeta * H * fMultiplier * phiFluxL(il) * phiFluxR(ir);
        }
    }
    for(int ir=0; ir<nfluxr; ir++) {
        for(int jr=0; jr<nfluxr; jr++) {
            ek(ir+secondblock,jr+secondblock) += weight *fBeta * H * phiFluxR(jr) * phiFluxR(ir);
        }
    }
}


/**
 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ek [out] is the stiffness matrix
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZMRCMLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
//	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
//	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	
//	TPZFNMatrix<660> dphiL, dphiR;
//	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
//	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	

	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
#ifdef PZDEBUG
    if(phiL.Rows()*fNStateVariables+phiR.Rows()*fNStateVariables != ek.Rows())
    {
        DebugStop();
    }
#endif
    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
	int il,jl,ir,jr;
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		for(jr=0; jr<nrowr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight * fMultiplier * (phiL(il) * phiR(jr));
            }
		}
	}
	
    //	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		for(jl=0; jl<nrowl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * fMultiplier * (phiR(ir) * phiL(jl));
            }
		}
	}
    
}

/**
 * @brief It computes a contribution to residual vector at one integration point
 * @param data [in]
 * @param dataleft [in]
 * @param dataright [in]
 * @param weight [in]
 * @param ef [out] is the load vector
 * @since April 16, 2007
 */
void TPZMRCMLagrangeMultiplier::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}


// print the data in human readable form
void TPZMRCMLagrangeMultiplier::Print(std::ostream &out)
{
    out << __PRETTY_FUNCTION__ << std::endl;
    TPZMaterial::Print(out);
    out << "NStateVariables " << this->fNStateVariables << std::endl;
    out << "fDimension " << this->fDimension << std::endl;
    out << "fMultiplier " << this->fMultiplier << std::endl;
    out << "fBeta " << this->fBeta << std::endl;
}


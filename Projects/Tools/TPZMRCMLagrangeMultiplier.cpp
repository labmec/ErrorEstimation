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
       || dataright[1].phi.Rows() == 0)
    {
        DebugStop();
    }
    
    TPZFMatrix<REAL> &phiFluxL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiFluxR = dataright[0].phi;
    TPZFMatrix<REAL> &phiPressR = dataright[1].phi;
    TPZFMatrix<REAL> &phiL = phiFluxL;
    TPZFMatrix<REAL> &phiR = phiPressR;
    REAL H = data.HSize;
    
    int nfluxl = phiL.Rows();
    int npressr = phiR.Rows();
    int nfluxr = phiFluxR.Rows();
    static int count  = 0;

    if((nfluxl+npressr+nfluxr)*fNStateVariables != ek.Rows() && count < 20)
    {
        std::cout<<"ek.Rows() "<< ek.Rows()<<
        " nrowl " << nfluxl <<
        " nrowr " << npressr << " may give wrong result " << std::endl;
        count++;
    }

    int secondblock = ek.Rows()-phiR.Rows()*fNStateVariables;
    
    // 3) phi_I_left, phi_J_right
    for(int il=0; il<nfluxl; il++) {
        for(int jr=0; jr<npressr; jr++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(fNStateVariables*il+ist,fNStateVariables*jr+ist+secondblock) += weight *  (phiFluxL(il) * phiR(jr));
            }
        }
    }
    
    //	// 4) phi_I_right, phi_J_left
    for(int ir=0; ir<npressr; ir++) {
        for(int jl=0; jl<nfluxl; jl++) {
            for (int ist=0; ist<fNStateVariables; ist++) {
                ek(ir*fNStateVariables+ist+secondblock,jl*fNStateVariables+ist) += weight * (phiR(ir) * phiL(jl));
            }
        }
    }

    if(fNStateVariables != 1) DebugStop();
    // penalty contribution
    for(int il=0; il<nfluxl; il++) {
        for(int jl=0; jl<nfluxl; jl++) {
            ek(il,jl) += weight * fBeta * H * phiFluxL(il) * phiFluxR(jl);
        }
        for(int ir=0; ir<nfluxr; ir++) {
            ek(ir+nfluxl,il) -= weight * fBeta * H * fMultiplier * phiFluxL(il) * phiFluxR(ir);
            ek(il,ir+nfluxl) -= weight * fBeta * H * fMultiplier * phiFluxL(il) * phiFluxR(ir);
        }
    }
    for(int ir=0; ir<nfluxr; ir++) {
        for(int jr=0; jr<nfluxr; jr++) {
            ek(ir+nfluxl,jr+nfluxl) += weight *fBeta * H * phiFluxR(jr) * phiFluxR(ir);
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



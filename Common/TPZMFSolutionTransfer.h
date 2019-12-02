//  TPZMFSolutionTransfer.hpp
//  LinearTracer
//
//  Created by Jose on 11/21/19.
//

#ifndef TPZMFSolutionTransfer_hpp
#define TPZMFSolutionTransfer_hpp

#include <stdio.h>
#include <iostream>
#include "TPZMultiphysicsCompMesh.h"
#include "pzsubcmesh.h"
#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


class TPZMFSolutionTransfer
{
   
    /**
     * @brief Structure that relates the origin blocks indexes and target indexes to transfer the solution.
     */
    struct Match{
        
        /**
         * @brief Origin block index.
         */
        int64_t fblocknumber;
        /**
         * @brief Pair that relates the origin block with the target block index.
         */
        std::pair<TPZBlock<STATE> *, int64_t> fblockTarget;// block target
        
        /**
         * @brief Empty Constructor
         */
        Match() : fblocknumber(-1), fblockTarget(0,-1)
        {
            
        }
        /**
         * @brief Copy Constructor
         */
        Match(const Match & copy){
            fblocknumber=copy.fblocknumber;
            fblockTarget = copy.fblockTarget;
        }
        
        /**
         * @brief operator =
         */
        
        Match &operator=(const Match &other){
            fblocknumber=other.fblocknumber;
            fblockTarget = other.fblockTarget;
            return *this;
        }
        
        /**
         * @brief Default Destructor
         */
        ~Match(){
        }
        
        /**
         * @brief Transfer the solution from the multiphysic mesh to the atomic meshes.
         * @param mfmesh target mesh.
         */
        
        void TransferFromMultiphysics(TPZCompMesh * mfmesh);
        /**
         * @brief Transfer the solution from the atomic meshes to the multiphysic mesh.
         * @param cmesh is the multiphysics mesh.
         * @note Internaly are taken the corresponding blocks.
         */
        void TransferToMultiphysics(TPZCompMesh * cmesh);
        
        /// Print the corresponding blocks
        void Print(std::ostream &out, TPZCompMesh *cmesh);
    };
    
    /**
     * @brief Structure that stock the matchs.
     */
    struct MeshTransferData{
        /**
         * @brief Origin mesh.
         */
        TPZCompMesh * fcmesh_ori;
        /**
         * @brief Matchs Vector.
         */
        TPZStack<Match> fconnecttransfer;
        
        /**
         * @brief Empty constructor
         */
        MeshTransferData(){
            
        }
        /**
         * @brief Copy constructor
         */
        MeshTransferData(const MeshTransferData & copy){
            fcmesh_ori=copy.fcmesh_ori;
            fconnecttransfer=copy.fconnecttransfer;
           
        }
        /**
         * @brief operator =
         */
        MeshTransferData &operator=(const MeshTransferData &other){
            fcmesh_ori=other.fcmesh_ori;
            fconnecttransfer=other.fconnecttransfer;
            return *this;
        }
        /**
         * @brief Default destructor
         */
        ~MeshTransferData(){
            
        }
        /**
         * @brief Build all teh matches between the multiphysic and the atomic meshes.
           @param cmesh Origin mesh
         */
        void BuildTransferData(TPZCompMesh* cmesh);
        
        /**
         * @brief Transfer the solution from the multiphysic mesh to the atomic meshes for every match stored in fconnecttransfer
         */
        void TransferFromMultiphysics();
        
        /**
         * @brief Transfer the solution from the atomic meshes to the multiphysic mesh for every match stored in fconnecttransfer
         */
        void TransferToMultiphysics();
       
        /**
         * @brief print the data structure
         */
        void Print(std::ostream &out);
    };
    public:
        /**
         * @brief Objects vector MeshTransferData, the transference should be done for the multiphysics mesh and their substructures.
         */
        TPZStack<MeshTransferData> fmeshTransfers;
        /**
         * @brief Build all the MeshTransferData for the multiphysic and their substructure.
         @param mfcmesh Origin mesh.
         */
        void BuildTransferData(TPZCompMesh* mfcmesh);
        /**
         * @brief Transfer the solution from the multiphysics mesh to the atomica meshes for every MeshTransferData stored in fmeshTransfers.
         */
        void TransferFromMultiphysics();
        /**
         * @brief Transfer teh solution from the atomic meshes to the multiphysic mesh for every MeshTransferData stored in fmeshTransfers.
         */
        void TransferToMultiphysics();
        /**
                Print the data structure
         */
    void Print(std::ostream &out = std::cout);
    
        
};

#endif /* TPZMFSolutionTransfer_hpp */

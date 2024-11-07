#pragma once

#include "PNPNodeData.h"
#include "PNPBaseEquation.h"
#include "PNPInputData.h"


class PNPManufacturedSoln : public PNPBaseEquation<PNPNodeData> {
 public:
      explicit PNPManufacturedSoln(int nsd, SurfaceIntegrationBehavior surface_Integration);

      void fillEssBC() override;
};


/// Constructor
PNPManufacturedSoln::PNPManufacturedSoln(int nsd, SurfaceIntegrationBehavior surface_Integration = ENABLE_SURFACE_INTEGRATION)
        : PNPBaseEquation<PNPNodeData>(nsd, surface_Integration) {
}

 void PNPManufacturedSoln::fillEssBC() {
     //const int nsd = p_grid_->nsd();  // TODO will be wrong when loading from Gmsh meshes
     const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
     const double t = this->t_;
     double dt = this->dt_;

     this->initEssBC();

     for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {

         ZEROPTV p = p_data_->p_grid_->GetNode(nodeID)->location();

         // Does the node lie on the left boundary?
         if (p_grid_->BoNode(nodeID, LEFT)){
             double phi = -1.0;
             p_data_->GetNodeData(nodeID).u[PNPNodeData::PHI_IDX] = phi;
             specifyValue(nodeID, PNPNodeData::PHI_IDX, 0.0);
         }

         // Does the node lie on the right boundary?
         if (p_grid_->BoNode(nodeID, RIGHT)){
             double phi = 1.0;
             p_data_->GetNodeData(nodeID).u[PNPNodeData::PHI_IDX] = phi;
             specifyValue(nodeID, PNPNodeData::PHI_IDX, 0.0);
         }

     }
 }












#pragma once

#include "PNPNodeData.h"
#include "PNPBaseEquation.h"
#include "PNPInputData.h"


class PNPManufacturedSoln : public PNPBaseEquation<PNPNodeData> {
 public:
      explicit PNPManufacturedSoln(int nsd);

      ZEROPTV velocityField(const FEMElm &fe, const double t) override;

      double NPForcing(const FEMElm &fe, const int species_id, double t) override;
      double PoissonForcing(const FEMElm &fe, double t) override;

      double calc_C_at(const ZEROPTV &location, const double &t, int species_idx);
      double calc_Phi_at(const ZEROPTV &location, const double &t);

      ZEROPTV calc_grad_C_at(const ZEROPTV &location, const double &t, int species_idx);
      ZEROPTV calc_grad_Phi_at(const ZEROPTV &location, const double &t);

      void fillEssBC() override;
      void calcbe_weak(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be) override;
};


/// Constructor
PNPManufacturedSoln::PNPManufacturedSoln(int nsd)
        : PNPBaseEquation<PNPNodeData>(nsd) {
}

 ZEROPTV PNPManufacturedSoln::velocityField(const FEMElm &fe, const double t) {
     return ZEROPTV(0.0, 0.0, 0.0);
 }

 double PNPManufacturedSoln::NPForcing(const FEMElm &fe, const int species_id, double t) {
     return 0.0;
 }

 double PNPManufacturedSoln::PoissonForcing(const FEMElm &fe, double t) {
     return 0.0;
 }

 void PNPManufacturedSoln::fillEssBC() {
     //const int nsd = p_grid_->nsd();  // TODO will be wrong when loading from Gmsh meshes
     const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
     const double t = this->t_;
     double dt = this->dt_;

     this->initEssBC();

     for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {

         //PrintStatus("Checking node: ", nodeID, "for Dirichlet");

         ZEROPTV p = p_data_->p_grid_->GetNode(nodeID)->location();

         // Does the node lie on the left boundary?
         if (p_grid_->BoNode(nodeID, LEFT)){
             double phi = -1.0;
             p_data_->GetNodeData(nodeID).u[PNPNodeData::PHI_IDX] = phi;
         }

         // Does the node lie on the right boundary?
         if (p_grid_->BoNode(nodeID, RIGHT)){
             double phi = 1.0;
             p_data_->GetNodeData(nodeID).u[PNPNodeData::PHI_IDX] = phi;
         }

     }
 }

 double PNPManufacturedSoln::calc_C_at(const ZEROPTV &location, const double &t, int species_idx) {
    return 0.0;
 }

 double PNPManufacturedSoln::calc_Phi_at(const ZEROPTV &location, const double &t) {
    return 0.0;
 }


 ZEROPTV PNPManufacturedSoln::calc_grad_C_at(const ZEROPTV &location, const double &t, int species_idx) {
    return ZEROPTV(0.0,0.0,0.0);
 }

 ZEROPTV PNPManufacturedSoln::calc_grad_Phi_at(const ZEROPTV &location, const double &t) {
     return ZEROPTV(0.0, 0.0, 0.0);
 }


 void PNPManufacturedSoln::calcbe_weak(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be) {

 }












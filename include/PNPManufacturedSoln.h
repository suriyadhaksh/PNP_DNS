 #pragma once
//#include <talyfem/equations/NSBaseEquation.h>
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

     const ZEROPTV p = fe.position();

     double u = cos(2 * M_PI * t)*sin(2*M_PI*p.x())*cos(2*M_PI*p.y());
     double v = - cos(2 * M_PI * t)*cos(2*M_PI*p.x())*sin(2*M_PI*p.y());
     double w = 0.0;

     return ZEROPTV(u, v, w);
 }

 double PNPManufacturedSoln::NPForcing(const FEMElm &fe, const int species_id, double t) {
     const ZEROPTV p = fe.position();
     double x = p.x();
     double y = p.y();

     double force = 0.0;
     switch (species_id) {
         case 0:
             force = 4*pow(M_PI,2)*z_[0]*pow(cos(2*M_PI*t),2)*pow(cos(2*M_PI*x),2)*pow(cos(2*M_PI*y),2)
                     +8*pow(M_PI,2)*cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y)
                     -2*M_PI*pow(cos(2*M_PI*t),2)*pow(cos(2*M_PI*x),2)*cos(2*M_PI*y)*sin(2*M_PI*y)
                     -2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*t)*sin(2*M_PI*y)
                     -2*M_PI*pow(cos(2*M_PI*t),2)*cos(2*M_PI*y)*pow(sin(2*M_PI*x),2)*sin(2*M_PI*y)
                     -8*pow(M_PI,2)*z_[0]*pow(cos(2*M_PI*t),2)*pow(cos(2*M_PI*x),2)*pow(sin(2*M_PI*y),2)
                     +4*pow(M_PI,2)*z_[0]*pow(cos(2*M_PI*t),2)*pow(sin(2*M_PI*x),2)*pow(sin(2*M_PI*y),2);

             break;

         case 1:
             force = 8*pow(M_PI,2)*cos(2*M_PI*t)*cos(2*M_PI*y)*sin(2*M_PI*x)
                     +2*M_PI*pow(cos(2*M_PI*t),2)*cos(2*M_PI*x)*pow(cos(2*M_PI*y),2)*sin(2*M_PI*x)
                     -2*M_PI*cos(2*M_PI*y)*sin(2*M_PI*t)*sin(2*M_PI*x)
                     -16*pow(M_PI,2)*z_[1]*pow(cos(2*M_PI*t),2)*cos(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*x)*sin(2*M_PI*y)
                     +2*M_PI*pow(cos(2*M_PI*t),2)*cos(2*M_PI*x)*sin(2*M_PI*x)*pow(sin(2*M_PI*y),2);

             break;

         default:
             PrintError("Forcing term: Species ID exceeds the number of Species");

     }

     return force;
 }

 double PNPManufacturedSoln::PoissonForcing(const FEMElm &fe, double t) {

     const ZEROPTV p = fe.position();
     double x = p.x();
     double y = p.y();

     double force = z_[1]*cos(2*M_PI*t)*cos(2*M_PI*y)*sin(2*M_PI*x)
             +16*pow(lambda_,2)*pow(M_PI,2)*cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y)
             +z_[0]*cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y);


     return force;
 }

 void PNPManufacturedSoln::fillEssBC() {
     //const int nsd = p_grid_->nsd();  // TODO will be wrong when loading from Gmsh meshes
     const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
     const int phi_idx = PNPNodeData::PHI_IDX;
     const double t = this->t_;
     double dt = this->dt_;

     this->initEssBC();

     //PrintStatus("Fill Essential BC initialised");

     for (int nodeID = 0; nodeID < this->p_grid_->n_nodes(); nodeID++) {

         //PrintStatus("Checking node: ", nodeID, "for Dirichlet");

         ZEROPTV p = p_data_->p_grid_->GetNode(nodeID)->location();

         for (int boundary = LEFT; boundary <= FRONT; boundary++) {
             //Does the node lie on the specific boundary?
             if (p_grid_->BoNode(nodeID, boundary)){

                 //Loop over all species
                 for (int species_idx = 0; species_idx < noOfSpecies; species_idx++) {
                     //Is this Dirichlet imposition?
                     if(BoundaryConditionArray_(boundary, species_idx) == DIRICHLET) {
                         double c = calc_C_at(p,  t+dt, species_idx);
                         p_data_->GetNodeData(nodeID).u[PNPNodeData::C_IDX + species_idx] = c;
                     }
                 }

                 //Check on Phi
                 if(BoundaryConditionArray_(boundary, phi_idx) == DIRICHLET) {
                     double phi = calc_Phi_at(p, t+dt);
                     p_data_->GetNodeData(nodeID).u[PNPNodeData::PHI_IDX] = phi;
                 }

             }

         }

     }
 }

 double PNPManufacturedSoln::calc_C_at(const ZEROPTV &location, const double &t, int species_idx) {

     double value = 0.0;
     double x = location.x();
     double y = location.y();

     switch (species_idx){
         case 0:
             value = cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y);
             break;

         case 1:
             value = cos(2*M_PI*t)*cos(2*M_PI*y)*sin(2*M_PI*x);
             break;

         default:
             PrintError("Calculate C at: Species ID exceeds the number of Species");
             return 0;
     }
     return value;
 }

 double PNPManufacturedSoln::calc_Phi_at(const ZEROPTV &location, const double &t) {

     double x = location.x();
     double y = location.y();
     return -1 * cos(2*M_PI*t)*cos(2*M_PI*x)*sin(2*M_PI*y);
 }


 ZEROPTV PNPManufacturedSoln::calc_grad_C_at(const ZEROPTV &location, const double &t, int species_idx) {
        double dC_dx = 0.0;
        double dC_dy = 0.0;
        double dC_dz = 0.0;

        double x = location.x();
        double y = location.y();

        switch (species_idx){
            case 0:
                dC_dx = -2*M_PI*cos(2*M_PI*t)*sin(2*M_PI*x)*sin(2*M_PI*y);
                dC_dy = 2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*cos(2*M_PI*y);

                break;

            case 1:
                dC_dx = 2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*cos(2*M_PI*y);
                dC_dy = -2*M_PI*cos(2*M_PI*t)*sin(2*M_PI*x)*sin(2*M_PI*y);

                break;

            default:
                PrintError("Calculate GradC at: Species ID exceeds the number of Species");
                return ZEROPTV(0.0, 0.0, 0.0);

        }

     return ZEROPTV(dC_dx, dC_dy, dC_dz);
 }

 ZEROPTV PNPManufacturedSoln::calc_grad_Phi_at(const ZEROPTV &location, const double &t) {

     double x = location.x();
     double y = location.y();

     double dPhi_dx = 2*M_PI*cos(2*M_PI*t)*sin(2*M_PI*x)*sin(2*M_PI*y);
     double dPhi_dy = -2*M_PI*cos(2*M_PI*t)*cos(2*M_PI*x)*cos(2*M_PI*y);
     double dPhi_dz = 0;


     return ZEROPTV(dPhi_dx, dPhi_dy, dPhi_dz);
 }


 void PNPManufacturedSoln::calcbe_weak(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be) {

     const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;

     const int nbf = fe.nbf();
     const double detSideJxW = fe.detJxW();

     const ZEROPTV &p = fe.position();
     const ZEROPTV &normal = fe.surface()->normal();

     const int C1Index = 0;
     if (BoundaryConditionArray_(sideInd, C1Index) == NEUMANN) {
         for (int a = 0; a < nbf; a++) {
             double gradCxNormal = calc_grad_C_at(p, t_ + dt_, C1Index).innerProduct(normal);
             be((noOfSpecies + 1) * a + C1Index) += - fe.N(a) * gradCxNormal * detSideJxW;
         }
     }

     const int C2Index = 1;
     if (BoundaryConditionArray_(sideInd, C2Index) == NEUMANN) {
         for (int a = 0; a < nbf; a++) {
             double gradCxNormal = calc_grad_C_at(p, t_ + dt_, C2Index).innerProduct(normal);
             be((noOfSpecies + 1) * a + C2Index) += - fe.N(a) * gradCxNormal * detSideJxW;
         }
     }

     const int phi_idx = PNPNodeData::PHI_IDX;
     if(BoundaryConditionArray_(sideInd, phi_idx) == NEUMANN) {
         const double lambda = lambda_;
         double gradPhixNormal = calc_grad_Phi_at(p,t_ + dt_).innerProduct(normal);
         for (int a = 0; a < nbf; a++){

             for (int i = 0; i < noOfSpecies; i++) {
                 double c = calc_C_at(p, t_ + dt_, i);
                 be((noOfSpecies + 1)*a + i) +=  - fe.N(a) * z_[i] * c * gradPhixNormal * detSideJxW;
             }

             be((noOfSpecies + 1)*a + phi_idx) += - 2 * lambda * lambda * fe.N(a) * gradPhixNormal * detSideJxW;

         }
     }


 }












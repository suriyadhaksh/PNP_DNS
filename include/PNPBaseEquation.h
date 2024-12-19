/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //

#pragma once

#include <limits>
#include <talyfem/fem/nonlinear_equation.h>
#include <talyfem/utils/macros.h>
#include <libconfig.h++>
#include "PNPInputData.h"


// used by PNPBaseEquation to detect missing fields in NodeData

DEFINE_HAS_MEMBER(has_PNP_C_IDX, C_IDX);
DEFINE_HAS_MEMBER(has_PNP_C_PREV_IDX, C_PREV_IDX);
DEFINE_HAS_MEMBER(has_PNP_C_PREV_2_IDX, C_PREV_2_IDX);
DEFINE_HAS_MEMBER(has_PNP_PHI_IDX, PHI_IDX);

enum BoundaryCondition {
    DIRICHLET = 0,
    NEUMANN = 1
};

enum BoundaryIndices {
    LEFT = 1,
    RIGHT = 2,
    BOTTOM = 3,
    TOP = 4,
    BACK = 5,
    FRONT = 6
};

template<typename NodeData>
class PNPBaseEquation : public NonlinearEquation<NodeData> {
    // give some more user-friendly compiler errors than "PNP_C_IDX not found"
    static_assert(has_PNP_C_IDX<NodeData>::value,
                  "NodeData is missing PNP_C_IDX. "
                  "To fix this, add something like enum { PNP_C_IDX = 0 }; inside your NodeData class.");
    static_assert(has_PNP_C_PREV_IDX<NodeData>::value,
                  "NodeData is missing PNP_C_PREV_IDX. "
                  "To fix this, add something like enum { PNP_C_PREV_IDX = 3 }; inside your NodeData class.");
    static_assert(has_PNP_C_PREV_2_IDX<NodeData>::value,
                  "NodeData is missing PNP_C_PREV_2_IDX. "
                  "To fix this, add something like enum { PNP_C_PREV_2_IDX = 5 }; inside your NodeData class.");
    static_assert(has_PNP_PHI_IDX<NodeData>::value,
                  "Your NodeData is missing PNP_PHI_IDX. "
                  "To fix this, add something like enum { PNP_PHI_IDX = 2 }; inside your NodeData class.");


    //Implementation of PNP. Currently only supports SuPG.


public:
    /**
     * @param nsd nsd of equation (e.g. input_data.nsd), used to know how many degrees of freedom are needed
     * @param surface_integration Set to "ENABLE_SURFACE_INTEGRATION" to enable weak BC. You must also override the enable_weak_bc method.
     * integrands4side is set by default for weak bc, if you want to use integrands4side for neumann bc, override integrands4side.
     */

    // indices to timing arrays. These are locations in the timers_ array that
    // correspond to specific stages of the code that we wish to time.
    static const int kTimerSolve = 0;  ///< index to timer for solve process
    static const int kTimerAssemble = 1;  ///< index to timer for assemble
    static const int kTimerUpdate = 2;  ///< index to timer for update process

    explicit PNPBaseEquation(int nsd, SurfaceIntegrationBehavior surface_integration = ENABLE_SURFACE_INTEGRATION)
            : NonlinearEquation<NodeData>(surface_integration) {

        for (int species_idx = 0; species_idx < PNPNodeData::NO_OF_SPECIES; species_idx++) {
            this->addDof(species_idx, PNPNodeData::C_IDX + species_idx);
        }
        this->addDof(PNPNodeData::PHI_IDX, PNPNodeData::PHI_IDX);

        // calculate second derivative (needed for SuPG)
        this->add_basis_flag(BASIS_SECOND_DERIVATIVE);

        timers_[kTimerSolve].set_label("Solve");
        timers_[kTimerAssemble].set_label("Assemble");
        timers_[kTimerUpdate].set_label("Update");
    }

    void Integrands(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) override {
        // PrintStatus("We are inside Integrands");
        if (use_stab_){
            IntegrandsNew(fe, Ae,  be);
        }
        else{
            IntegrandsGeneric(fe, Ae, be);
        }

        // IntegrandsGeneric(fe, Ae, be);
    }


    void Integrands4side(const FEMElm &fe, const int sideInd, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) override {
        // calcAe_weak(fe, sideInd, Ae);
        // calcbe_weak(fe, sideInd, be);
    }

    void copyBoundaryConditions(const ZeroMatrix<int> &BoundaryConditionArray) {
        BoundaryConditionArray_ = BoundaryConditionArray;
    }

    void setParams(const double lambda, double gridSpacing, double * z, bool use_stab) {
        lambda_ = lambda;
        h_ = gridSpacing;
        z_ = z;
        use_stab_ = use_stab;
    }

    void IntegrandsSUPG(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) {

        const int nsd = fe.nsd();
        const int n_basis_functions = fe.nbf();
        const double detJxW = fe.detJxW();
        const double t = this->t_;
        double dt = this->dt_;
        const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
        const int phi_idx = PNPNodeData::PHI_IDX;

        const double h = h_;

        ZEROPTV velocity = velocityField(fe, t + dt);

        std::vector<double> c (noOfSpecies);
        std::vector<double> c_prev (noOfSpecies);
        std::vector<double> c_prev_2 (noOfSpecies);
        std::vector<double> forcingC(noOfSpecies);

        ZeroMatrix<double> dc;
        ZeroMatrix<double> netVelocity;
        dc.redim(noOfSpecies, nsd);
        netVelocity.redim(noOfSpecies, nsd);

        ZEROPTV dphi;
        ZEROPTV d2phi_dir_dir;
        for (int dir = 0; dir < nsd; dir++) {
            dphi(dir) = this->p_data_->valueDerivativeFEM(fe, phi_idx, dir);
            d2phi_dir_dir(dir) = this->p_data_->value2DerivativeFEM(fe, phi_idx, dir, dir);
        }

        std::vector<double> nernstPlanck(noOfSpecies);
        std::vector<double> tau(noOfSpecies);

        for (int species_idx = 0; species_idx < noOfSpecies; species_idx++){
            c[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_IDX + species_idx);
            c_prev[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_IDX + species_idx);
            c_prev_2[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_2_IDX + species_idx);
            forcingC[species_idx] = 0.0;

            for(int dir = 0; dir < nsd; dir ++) {
                dc(species_idx, dir) = this->p_data_->valueDerivativeFEM(fe,PNPNodeData::C_IDX + species_idx,dir);
            }

            //To calculate residual
            double temporal = ( 3 * c[species_idx] - 4 * c_prev[species_idx] + c_prev_2[species_idx] ) / (2 * dt);
            double advection = 0.0;
            double diffusion = 0.0;
            double electromobility = 0.0;

            //SUPG
            double netVelocitySumOfSquares = 0.0;

            for(int dir = 0; dir < nsd; dir++) {
                advection += velocity(dir) * dc(species_idx, dir);
                //diffusion += this->p_data_->value2DerivativeFEM(fe, PNPNodeData::C_IDX + species_idx, dir, dir);
                electromobility += z_[species_idx]*dphi(dir)*dc(species_idx,dir);
                //electromobility += z_[species_idx]*c[species_idx]*d2phi_dir_dir(dir)
                //+ z_[species_idx]*dphi(dir)*dc(species_idx,dir);

                netVelocity(species_idx,dir) = velocity(dir) - z_[species_idx] * dphi(dir);
                netVelocitySumOfSquares += netVelocity(species_idx,dir) * netVelocity(species_idx,dir);

            }

            //nernstPlanck[species_idx] = temporal + advection - diffusion - electromobility - forcingC[species_idx];
            nernstPlanck[species_idx] = temporal + advection - electromobility;

            //SUPG calculation
            double net_vel_mag = sqrt(netVelocitySumOfSquares);
            double Pe = net_vel_mag * h / (2  * 1); //1 is the diffusivity coeff

            tau[species_idx] = 1 / sqrt( 4/(dt*dt) + 4* net_vel_mag * net_vel_mag/(h*h) + 9/(Pe*Pe));
            //tau[species_idx] = 0.0;


        }

        //double forcePhi = PoissonForcing(fe, t + dt);
        double forcePhi = 0.0;


        //Stiffness Matrix
        for (int a = 0; a < n_basis_functions; a++) {

            for (int b = 0; b < n_basis_functions; b++) {

                double F = fe.N(a) * fe.N(b) * detJxW;

                double G = 0.0;
                double H = 0.0;
                double I = 0.0;


                for (int dir = 0; dir < nsd; dir++) {

                    G += fe.N(a) * velocity(dir) * fe.dN(b, dir) * detJxW;
                    H += fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    I += fe.dN(a, dir) * dphi(dir) * fe.N(b) * detJxW;
                }

                for (int i = 0; i < noOfSpecies; i++) {

                    double J = 0.0;
                    double K = 0.0;
                    double L = 0.0;
                    for (int dir = 0; dir < nsd; dir++){
                        J += fe.dN(a,dir) * netVelocity(i,dir);
                        K += fe.dN(b,dir) * netVelocity(i,dir) * detJxW;
                        L += fe.dN(b,dir) * z_[i] * dc(i,dir)  * detJxW;
                    }

                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + i) += 3 * F / (2 * dt)
                                                                            + G
                                                                            + H
                                                                            + I * z_[i]
                                                                            + 3 * J * tau[i] * fe.N(b) * detJxW / (2 * dt);
                    + J * K * tau[i] ;

                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + phi_idx) += H * z_[i] * c[i]
                                                                                  - J * L * tau[i]
                                                                                  + H * tau[i] * z_[i] * nernstPlanck[i];


                    Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + i) += - F * z_[i];

                    //Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + i) += F * z_[i]
                    //+ J * tau[i] * fe.N(b) * z_[i] * detJxW;
                }

                Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + phi_idx) += 2 * lambda_ * lambda_ * H;

            }

            //RHS vector

            double M = 0.0;
            double U = 0.0;
            for (int i = 0; i < noOfSpecies; i++) {
                double N = fe.N(a) * ( 3 * c[i] - 4 * c_prev[i] + c_prev_2[i] ) / 2 * detJxW;
                double O = 0.0;
                double P = 0.0;
                double Q = 0.0;
                double R = 0.0;

                for (int dir = 0; dir < nsd; dir++) {
                    O += fe.N(a) * velocity(dir) * dc(i,dir) * detJxW;
                    P += fe.dN(a,dir) * dc(i,dir) * detJxW;
                    Q += fe.dN(a,dir) * netVelocity(i,dir) * tau[i] * nernstPlanck[i] * detJxW;
                    R += fe.dN(a,dir) * z_[i] * c[i] * dphi(dir) * detJxW;
                }

                double forceNP = fe.N(a)  * forcingC[i] * detJxW;

                be((noOfSpecies + 1)*a + i) += N / dt + (O + P + Q + R) - forceNP;

                //for the poisson equation
                M += z_[i] * c[i];
            }

            double S = 0.0;
            double delta = 0.0;
            for(int dir = 0; dir < nsd; dir++) {
                S += fe.dN(a,dir) * dphi(dir);
                delta += fe.dN(a, dir) * netVelocity(0,dir) * tau[0];
            }

            /*be((noOfSpecies + 1)*a + phi_idx) += ( fe.N(a) + delta ) * M * detJxW
                    - 2 * lambda_ * lambda_ * S * detJxW
                    - (fe.N(a) + delta ) * forcePhi * detJxW;*/

            /*
            be((noOfSpecies + 1)*a + phi_idx) += - fe.N(a) * M * detJxW
                                                 + 2 * lambda_ * lambda_ * S * detJxW
                                                 - (fe.N(a) ) * forcePhi * detJxW;
                                                 */
            be((noOfSpecies + 1)*a + phi_idx) += - fe.N(a) * M * detJxW
                                                 + 2 * lambda_ * lambda_ * S * detJxW;
        }

    }

    // Function created for MATH666 coursework. Supg implementation
    void IntegrandsNew(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be){
        const int nsd = fe.nsd();
        const int n_basis_functions = fe.nbf();
        const double detJxW = fe.detJxW();
        const double h = h_;
        const double t = this->t_;
        double dt = this->dt_;
        const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
        const int phi_idx = PNPNodeData::PHI_IDX;
        ZEROPTV pt = fe.position();

        // ------------ solution guess vectors & arrays ------------ //
        std::vector<double> c (noOfSpecies);
        std::vector<double> c_prev (noOfSpecies);
        std::vector<double> c_prev_2 (noOfSpecies);
        ZeroMatrix<double> dc;
        dc.redim(noOfSpecies, nsd);
        ZEROPTV dphi;

        // ------------ fill all solution guess vectors & arrays ------------ //
        for (int dir = 0; dir < nsd; dir++) {
            dphi(dir) = this->p_data_->valueDerivativeFEM(fe, phi_idx, dir);
        }

        for (int species_idx = 0; species_idx < noOfSpecies; species_idx++){

            c[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_IDX + species_idx);
            c_prev[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_IDX + species_idx);
            c_prev_2[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_2_IDX + species_idx);

            for(int dir = 0; dir < nsd; dir ++) {
                dc(species_idx, dir) = this->p_data_->valueDerivativeFEM(fe,PNPNodeData::C_IDX + species_idx, dir);
            }
        }

        // ----------- SUPG terms ------------- //
        std::vector<double> res(noOfSpecies);
        std::vector<double> tau(noOfSpecies);

        // ------------ fill all SuPG guess vectors & arrays ------------ //
        for (int species_idx = 0; species_idx < noOfSpecies; species_idx++){
            res[species_idx] = (c[species_idx] - c_prev[species_idx])/dt;
            double normOfVelocity = 0.0;
            for(int dir = 0; dir < nsd; dir ++) {
                normOfVelocity += (z_[species_idx] * dphi(dir)) * (z_[species_idx] * dphi(dir));
                res[species_idx] += - z_[species_idx] * dphi(dir) * dc(species_idx, dir);
            }
            normOfVelocity = sqrt(normOfVelocity);
            double Pe = normOfVelocity * h / (2  * 1); //1 is the diffusivity coeff
            tau[species_idx] = 1 / sqrt( 4/(dt*dt) + 4* normOfVelocity * normOfVelocity/(h*h) + 9/(Pe*Pe));
        }


        // ------------ assemble the jacobian ------------ //
        for (int a = 0; a < n_basis_functions; a++) {

            for (int b = 0; b < n_basis_functions; b++) {

                // Nernst Planck terms
                //double Dtemporal_DC = 3 * fe.N(a) * fe.N(b) / (2 * dt) * detJxW;
                double Dtemporal_DC = fe.N(a) * fe.N(b) / dt * detJxW;
                double Ddiffusion_DC = 0.0;
                double Delectromigration_DC = 0.0;
                double Delectromigration_DPhi = 0.0;

                // supg terms
                std::vector<double> dres_dC(noOfSpecies);
                std::vector<double> dres_dphi(noOfSpecies);
                double Dsupg_Dres = 0.0;

                for (int i = 0; i < noOfSpecies; i++){
                    dres_dC[i] = fe.N(b)/dt;
                    dres_dphi[i] = 0.0;
                    for (int dir = 0; dir < nsd; dir++){
                        dres_dC[i] += - z_[i]*fe.dN(b,dir) * dphi(dir);
                        dres_dphi[i] += - z_[i] * dc(i,dir) * fe.dN(b, dir);
                    }
                }

                // Poisson terms
                double Dchargedensity_DC = fe.N(a) * fe.N(b) * detJxW;
                double Dphilaplacian_Dphi = 0.0;

                for (int dir = 0; dir < nsd; dir++) {
                    Ddiffusion_DC += fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    Delectromigration_DC += fe.dN(a, dir) * dphi(dir) * fe.N(b) * detJxW;
                    Delectromigration_DPhi += fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    Dphilaplacian_Dphi += 2 * lambda_ * lambda_ * fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    Dsupg_Dres += fe.dN(a,dir) * dphi(dir) * detJxW;
                }

                for (int i = 0; i < noOfSpecies; i++) {
                    // Nernst Planck
                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + i) += Dtemporal_DC + Ddiffusion_DC + Delectromigration_DC * z_[i]
                            - Delectromigration_DC * z_[i] * tau[i] * res[i]
                            - Dsupg_Dres * z_[i] * c[i] * tau[i] * dres_dC[i];
                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + phi_idx) += Delectromigration_DPhi * z_[i] * c[i]
                            - Delectromigration_DPhi * z_[i] * c[i] * tau[i] * res[i]
                            - Dsupg_Dres * z_[i] * tau[i] * dres_dphi[i];

                    // Poisson equation
                    Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + i) += - Dchargedensity_DC  * z_[i];

                }

                // Poisson equation
                Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + phi_idx) += Dphilaplacian_Dphi;

            }

            // ------------ assemble the rhs ------------ //
            double chargedensity = 0.0;
            for (int i = 0; i < noOfSpecies; i++) {
                // double temporal = fe.N(a) * ( 3 * c[i] - 4 * c_prev[i] + c_prev_2[i] ) / (2 * dt)* detJxW;
                double temporal = fe.N(a) * (c[i] - c_prev[i]) / dt * detJxW;
                double diffusion = 0.0;
                double electromigration = 0.0;
                double supg = 0.0;

                for (int dir = 0; dir < nsd; dir++) {
                    diffusion += fe.dN(a, dir) * dc(i, dir) * detJxW;
                    electromigration += fe.dN(a, dir) * z_[i] * c[i] * dphi(dir) * detJxW;
                    supg += - tau[i] * res[i] * z_[i] * fe.dN(a, dir) * dphi(dir) * detJxW;
                    double dphi_1 = dphi(dir);
                }

                // Nernst planck
                be((noOfSpecies + 1)*a + i) += temporal + diffusion + electromigration + supg;

                //for the poisson equation
                chargedensity += fe.N(a) * z_[i] * c[i] * detJxW;
            }

            double philaplacian = 0.0;
            for(int dir = 0; dir < nsd; dir++) {
                philaplacian += 2 * lambda_ * lambda_ * fe.dN(a, dir) * dphi(dir) * detJxW;
            }

            // poisson
            be((noOfSpecies + 1)*a + phi_idx) += - chargedensity + philaplacian;
        }
    }

    void IntegrandsGeneric(const FEMElm &fe, ZeroMatrix<double> &Ae, ZEROARRAY<double> &be) {

        const int nsd = fe.nsd();
        const int n_basis_functions = fe.nbf();
        const double detJxW = fe.detJxW();
        const double t = this->t_;
        double dt = this->dt_;
        const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
        const int phi_idx = PNPNodeData::PHI_IDX;
        ZEROPTV pt = fe.position();

        // ------------ solution guess vectors & arrays ------------ //
        std::vector<double> c (noOfSpecies);
        std::vector<double> c_prev (noOfSpecies);
        std::vector<double> c_prev_2 (noOfSpecies);
        ZeroMatrix<double> dc;
        dc.redim(noOfSpecies, nsd);
        ZEROPTV dphi;

        // ------------ fill all solution guess vectors & arrays ------------ //
        for (int dir = 0; dir < nsd; dir++) {
            dphi(dir) = this->p_data_->valueDerivativeFEM(fe, phi_idx, dir);
        }

        for (int species_idx = 0; species_idx < noOfSpecies; species_idx++){

            c[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_IDX + species_idx);
            c_prev[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_IDX + species_idx);
            c_prev_2[species_idx] = this->p_data_->valueFEM(fe, PNPNodeData::C_PREV_2_IDX + species_idx);

            for(int dir = 0; dir < nsd; dir ++) {
                dc(species_idx, dir) = this->p_data_->valueDerivativeFEM(fe,PNPNodeData::C_IDX + species_idx, dir);
            }
        }


        // ------------ assemble the jacobian ------------ //
        for (int a = 0; a < n_basis_functions; a++) {

            for (int b = 0; b < n_basis_functions; b++) {

                // Nernst Planck terms
                //double Dtemporal_DC = 3 * fe.N(a) * fe.N(b) / (2 * dt) * detJxW;
                double Dtemporal_DC = fe.N(a) * fe.N(b) / dt * detJxW;
                double Ddiffusion_DC = 0.0;
                double Delectromigration_DC = 0.0;
                double Delectromigration_DPhi = 0.0;

                // Poisson terms
                double Dchargedensity_DC = fe.N(a) * fe.N(b) * detJxW;
                double Dphilaplacian_Dphi = 0.0;

                for (int dir = 0; dir < nsd; dir++) {
                    Ddiffusion_DC += fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    Delectromigration_DC += fe.dN(a, dir) * dphi(dir) * fe.N(b) * detJxW;
                    Delectromigration_DPhi += fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                    Dphilaplacian_Dphi += 2 * lambda_ * lambda_ * fe.dN(a, dir) * fe.dN(b, dir) * detJxW;
                }

                for (int i = 0; i < noOfSpecies; i++) {
                    // Nernst Planck
                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + i) += Dtemporal_DC + Ddiffusion_DC + Delectromigration_DC * z_[i];
                    Ae((noOfSpecies + 1)*a + i, (noOfSpecies + 1)*b + phi_idx) += Delectromigration_DPhi * z_[i] * c[i];

                    // Poisson equation
                    Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + i) += - Dchargedensity_DC  * z_[i];

                }

                // Poisson equation
                Ae((noOfSpecies + 1)*a + phi_idx, (noOfSpecies + 1)*b + phi_idx) += Dphilaplacian_Dphi;

            }

            // ------------ assemble the rhs ------------ //
            double chargedensity = 0.0;
            for (int i = 0; i < noOfSpecies; i++) {
                // double temporal = fe.N(a) * ( 3 * c[i] - 4 * c_prev[i] + c_prev_2[i] ) / (2 * dt)* detJxW;
                double temporal = fe.N(a) * (c[i] - c_prev[i]) / dt * detJxW;
                double diffusion = 0.0;
                double electromigration = 0.0;

                for (int dir = 0; dir < nsd; dir++) {
                    diffusion += fe.dN(a, dir) * dc(i, dir) * detJxW;
                    electromigration += fe.dN(a, dir) * z_[i] * c[i] * dphi(dir) * detJxW;
                    double dphi_1 = dphi(dir);
                }

                // Nernst planck
                be((noOfSpecies + 1)*a + i) += temporal + diffusion + electromigration;

                //for the poisson equation
                chargedensity += fe.N(a) * z_[i] * c[i] * detJxW;
            }

            double philaplacian = 0.0;
            for(int dir = 0; dir < nsd; dir++) {
                philaplacian += 2 * lambda_ * lambda_ * fe.dN(a, dir) * dphi(dir) * detJxW;
            }

            // poisson
            be((noOfSpecies + 1)*a + phi_idx) += - chargedensity + philaplacian;
        }
    }


private:

    double h_;
    MPITimer timers_[3];


protected:

    double * z_;
    double lambda_;
    bool use_stab_;
    ZeroMatrix<int> BoundaryConditionArray_ ;

    virtual ZEROPTV velocityField(const FEMElm &fe, const double t) {
        return ZEROPTV(0.0, 0.0, 0.0);
    }

    virtual void calcAe_weak(const FEMElm &fe, int sideInd, ZeroMatrix<double> &Ae) {
        //design this function in the inherited manufacturedsol class
    }

    virtual void calcbe_weak(const FEMElm &fe, int sideInd, ZEROARRAY<double> &be) {
        //design this function in the inherited manufacturedsol class
    }

};
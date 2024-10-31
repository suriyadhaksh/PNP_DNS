/*
  Copyright 2014-2017 Baskar Ganapathysubramanian

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

#include "PNPNodeData.h"
#include "PNPInputData.h"
#include "PNPManufacturedSoln.h"
#include <random>

class PNPGridField : public GridField<PNPNodeData> {
 public:
  PNPGridField(PNPInputData *input_data_in) : input_data_(input_data_in) {};

  void SetIC(int nsd) {

      const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
      const int phi_idx = PNPNodeData::PHI_IDX;
      const int c_mms_idx = PNPNodeData::C_MMS_IDX;
      const int phi_mms_idx = PNPNodeData::PHI_MMS_IDX;

      for (int nodeID = 0; nodeID < p_grid_->n_nodes(); nodeID++) {
          PNPNodeData* pData = &(GetNodeData(nodeID));
          ZEROPTV pt = p_grid_->GetNode(nodeID)->location();

          for (int i = 0; i < noOfSpecies; i++) {
              //double c_init = pnpeq_->calc_C_at(pt, 0, i);
              double c_init = 1.0; // for leading order case
              pData->u[i] = c_init;
              pData->u_prev[i] = c_init;
              pData->u_prev2[i] = c_init;
              pData->u[c_mms_idx + i] = c_init;
          }

          // double phi_init = pnpeq_->calc_Phi_at(pt,0);
          double phi_init = 0.0;
          pData->u[phi_idx] = phi_init;
          pData->u[phi_mms_idx] = phi_init;
      }
  }

  void setMMSGridField(double t) {
      const int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
      const int c_mms_idx = PNPNodeData::C_MMS_IDX;
      const int phi_mms_idx = PNPNodeData::PHI_MMS_IDX;

      for (int node_id = 0; node_id < p_grid_->n_nodes(); node_id++) {
          PNPNodeData* pData = &(GetNodeData(node_id));
          ZEROPTV p = p_grid_->GetNode(node_id)->location();
          for (int i = 0; i < noOfSpecies; i++) {
              double c_mms = pnpeq_->calc_C_at(p, t, i);
              pData->value(c_mms_idx+i) = c_mms;
          }
          double phi_mms = pnpeq_->calc_Phi_at(p,t);
          pData->value(phi_mms_idx) = phi_mms;
      }
  }

  void setMMS(PNPManufacturedSoln *pnpeq) {
      pnpeq_ = pnpeq;
  }

  void PrintTotalConcentration(double t){
      std::vector<double> TotalConcentration;
      std::fstream filePtErrorManufacSol;

      int noOfSpecies = PNPNodeData::NO_OF_SPECIES;

      std::string fileName = "Total_C_data.plt";

      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (!rank) {
          TotalConcentration = CalcTotalConcentration(this->input_data_, t, this->input_data_->ifDD);

          filePtErrorManufacSol.open(fileName, std::ios::app);

          filePtErrorManufacSol << t;
          for (int i = 0; i < noOfSpecies; i++) {
              filePtErrorManufacSol << "\t" << TotalConcentration[i];
          }
          filePtErrorManufacSol <<  std::endl;

          for (int i = 0; i < noOfSpecies; i++) {
              std::string str = "Total_C" + std::to_string(i) + " = ";
              PrintInfo(str, TotalConcentration[i]);
          }
          filePtErrorManufacSol.close();
      }

  }

  /// Compare with manufactured solution
  void PrintError(double t) {
    std::vector<double> L2Error;
    std::fstream filePtErrorManufacSol;

    int noOfSpecies = PNPNodeData::NO_OF_SPECIES;

    std::string fileName = "Error_data_" + std::to_string(noOfSpecies) + "_Species.plt";

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank) {
      L2Error = CalcL2Error(this->input_data_, t, this->input_data_->ifDD);

      filePtErrorManufacSol.open(fileName, std::ios::app);

      filePtErrorManufacSol << t;
      for (int i = 0; i < noOfSpecies; i++) {
          filePtErrorManufacSol << "\t" << L2Error[i];
      }
      filePtErrorManufacSol << "\t" << L2Error[PNPNodeData::PHI_IDX] <<  std::endl;

      for (int i = 0; i < noOfSpecies; i++) {
          std::string str = "L2Error_C" + std::to_string(i) + " = ";
          PrintInfo(str, L2Error[i]);
      }

      PrintInfo("L2Error_Phi = ", L2Error[PNPNodeData::PHI_IDX]);
      filePtErrorManufacSol.close();
    }
  }

  std::vector<double> CalcTotalConcentration(const InputData *input_data, double t,
                                               bool ifDD) const {

      FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);

      int noOfSpecies = PNPNodeData::NO_OF_SPECIES;

      std::vector<double> C_total(noOfSpecies, 0.0);

      const double n_elements = p_grid_->n_elements();
      for (int elm_id = 0; elm_id < n_elements; elm_id++) {
          fe.refill(elm_id, input_data->basisFunction, 0);
          while (fe.next_itg_pt()) {
              const double detJxW = fe.detJxW();

              for (int i = 0; i < noOfSpecies; i++) {
                  double calculatedCSol = valueFEM(fe, i);
                  C_total[i] +=  calculatedCSol * detJxW;
              }

          }
      }

      return C_total;
  }


/** Calculates the Relative Error when the solution for the test case is
* compared with manufactured solutions
* Returned vector is always in 3D format, <u, v, w, p>, in 2D, w is zero.
* @param input_data
* @param t
* @return returns a std::vector of errors
*/
  std::vector<double> CalcL2Error(const InputData *input_data, double t,
                                  bool ifDD) const {
    FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);

    int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
    int phi_idx = PNPNodeData::PHI_IDX;
    int errorVectorSize = PNPNodeData::NUM_VARS;

    std::vector<double> l2_error(errorVectorSize, 0.0);

    const double n_elements = p_grid_->n_elements();
    for (int elm_id = 0; elm_id < n_elements; elm_id++) {
      fe.refill(elm_id, input_data->basisFunction, 0);
      while (fe.next_itg_pt()) {
        const double detJxW = fe.detJxW();


        for (int i = 0; i < noOfSpecies; i++) {
            double manufacturedCSol = pnpeq_->calc_C_at(fe.position(),t,i);
            double calculatedCSol = valueFEM(fe, i);
            l2_error[i] += (manufacturedCSol - calculatedCSol) * (manufacturedCSol - calculatedCSol) * detJxW;
        }

        double manufacturedPhiSol = pnpeq_->calc_Phi_at(fe.position(),t);
        double calculatedPhiSol = valueFEM(fe, phi_idx);
        l2_error[phi_idx] += (manufacturedPhiSol - calculatedPhiSol) * (manufacturedPhiSol - calculatedPhiSol) * detJxW;

      }
    }

    for (int i = 0; i < noOfSpecies; i++) {
        l2_error[i] = sqrt(l2_error[i]);
    }
    l2_error[phi_idx] = sqrt(l2_error[phi_idx]);

    return l2_error;

  }


  std::vector<double> CalcRelL2Error(const InputData *input_data, double t,
                                    bool ifDD) const {
      FEMElm fe(p_grid_, BASIS_FIRST_DERIVATIVE | BASIS_POSITION);

      int noOfSpecies = PNPNodeData::NO_OF_SPECIES;
      int phi_idx = PNPNodeData::PHI_IDX;
      int errorVectorSize = PNPNodeData::NUM_VARS;

      std::vector<double> l2_error(errorVectorSize);
      std::vector<double> absManufacSol(errorVectorSize);
      std::vector<double> relError(errorVectorSize);

      const double n_elements = p_grid_->n_elements();
      for (int elm_id = 0; elm_id < n_elements; elm_id++) {
          fe.refill(elm_id, input_data->basisFunction, 0);
          while (fe.next_itg_pt()) {
              const double detJxW = fe.detJxW();

              std::vector<double> manufacSol(errorVectorSize);
              std::vector<double> calcSol(errorVectorSize);

              for (int i = 0; i < noOfSpecies; i++) {
                  PrintStatus("Calculating error for species ", i);
                  double localL2error = 0.0;
                  double manufacturedCSol = pnpeq_->calc_C_at(fe.position(),t,i);
                  double calculatedCSol = valueFEM(fe, i);
                  absManufacSol[i] += manufacturedCSol * manufacturedCSol * detJxW;
                  l2_error[i] += (manufacturedCSol - calculatedCSol) * (manufacturedCSol - calculatedCSol) * detJxW;
              }

              double manufacturedPhiSol = pnpeq_->calc_Phi_at(fe.position(),t);
              double calculatedPhiSol = valueFEM(fe, phi_idx);
              absManufacSol[phi_idx] += manufacturedPhiSol * manufacturedPhiSol * detJxW;
              l2_error[phi_idx] += (manufacturedPhiSol - calculatedPhiSol) * (manufacturedPhiSol - calculatedPhiSol) * detJxW;

          }
      }

      double eps = 1e-12;

      for (int i = 0; i < errorVectorSize; i++) {
          l2_error[i] = sqrt(l2_error[i]);

          if (sqrt(absManufacSol[i] >= eps)) {
              relError[i] = l2_error[i] / sqrt(absManufacSol[i]);
          }

          else {
              PrintInfo("The absolute sum of dof no.", i, " is close to zero! This is normal in 2D");
              relError[i] = l2_error[i];
          }
      }

      return relError;
  }


 private:
  PNPInputData *input_data_;
  PNPManufacturedSoln* pnpeq_;
};

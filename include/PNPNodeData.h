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

class PNPNodeData {
 public:
  enum ValueIndex {

    NO_OF_SPECIES = 2, //No of species

    C_IDX = 0, //Species Concentration
    PHI_IDX = NO_OF_SPECIES, //Electrical Potential

    C_MMS_IDX = PHI_IDX + 1, //Species Concentration - MMS
    PHI_MMS_IDX = C_MMS_IDX + NO_OF_SPECIES, //Electric Potential - MMS

    NUM_VARS = (NO_OF_SPECIES + 1) * 2, //No of variables (Includes variables and MMS)

    C_PREV_IDX = C_IDX + NUM_VARS,

    C_PREV_2_IDX = C_PREV_IDX + NO_OF_SPECIES,

  };

  double u[NUM_VARS];
  double u_prev[NO_OF_SPECIES]; //U^{n-1} Array
  double u_prev2[NO_OF_SPECIES]; //U^{n-2} Array

  inline double& value(int index) {
    if (index >= 0 && index < C_PREV_IDX) // 0, 1 = C, 2 = PHI, 3, 4 = C_MMS, 5 = PHI_MMS
      return u[index];
    else if (index >= C_PREV_IDX && index < C_PREV_2_IDX)  //6, 7 = C_PREV
      return u_prev[index - C_PREV_IDX];
    else if (index >= C_PREV_2_IDX && index < C_PREV_2_IDX + NO_OF_SPECIES)  //8, 9 = C_PREV_PREV
      return u_prev2[index - C_PREV_2_IDX];
    else
      throw TALYException() << "Invalid NodeData index";
  }

  inline const double& value(int index) const {
    return const_cast<PNPNodeData*>(this)->value(index);
  }

  static const char* name(int index) {
    switch (index) {
        case C_IDX:
            return "C1";
            break;

        case C_IDX + 1:
            return "C2";
            break;

        case PHI_IDX:
            return "Phi";
            break;

        case C_MMS_IDX:
            return "C1_mms";
            break;

        case C_MMS_IDX + 1:
            return "C2_mms";
            break;

        case PHI_MMS_IDX:
            return "Phi_mms";
            break;

        default:
            throw TALYException() << "Invalid NodeData index";
    }

  }

  static int valueno() {
    return NUM_VARS;
  }

  void UpdateDataStructures() {
      for (int idx = 0; idx < NO_OF_SPECIES; idx++) {
          u_prev2[idx] = u_prev[idx];
          u_prev[idx] = u[idx];
      }
  }
};

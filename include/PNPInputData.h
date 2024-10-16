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

#include <talyfem/input_data/input_data.h>
#include <string>
#include <PNPNodeData.h>
#include <PNPManufacturedSoln.h>


struct PNPInputData : public InputData {
  bool ifPrintPltFiles;  ///< whether to print .plt files at start and end
  bool shouldFail;  ///< whether or not this run.sh should fail (for unit testing)
  bool use_bdf2_; //?<whether to use BDF2 or not. This code has no implementation for non-BDF2 cases
  double dt;
  double totalT;
  double lambda; //Normalised Debye length
  double z[PNPNodeData::NO_OF_SPECIES]; //Valency array of species

  PNPInputData()
      : InputData(),
        ifPrintPltFiles(true),
        shouldFail(false),
        totalT(1.0),
        dt(1e-2),
        use_bdf2_(true){

      for (int idx = 0; idx < PNPNodeData::NO_OF_SPECIES; idx++) {
          z[idx] = 1.0;
      }
  }



  bool CheckInputData() const {
    if (((typeOfIC == 0) && (inputFilenameGridField == ""))
        || ((typeOfIC != 0) && (inputFilenameGridField != ""))) {
      PrintWarning("IC not set properly check!", typeOfIC, " ",
                   inputFilenameGridField);
      return false;
    }

    return InputData::CheckInputData();
  }

  ZeroMatrix<int> BoundaryConditionArray;


  bool ReadFromFile(const std::string& filename = std::string("config.txt")) {
    // Read config file and initialize basic fields
    InputData::ReadFromFile(filename);    // read the input file

    ReadValue("ifPrintPltFiles", ifPrintPltFiles);
    ReadValue("shouldFail", shouldFail);
    ReadValue("dt", dt);
    ReadValue("use_bdf2", use_bdf2_);
    ReadValue("totalT", totalT);
    ReadValue("dbLength", lambda); // read normalised Debye length

    for (int IDX = 0; IDX < PNPNodeData::NO_OF_SPECIES; IDX++) {
        std::string str = "val_" + std::to_string(IDX); // read valency
        ReadValue(str, z[IDX]);
    }

    setBoundary();

    return true;
  }

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


  void setBoundary () {

      const int noOfBoundaries = 6;
      BoundaryConditionArray.redim(noOfBoundaries + 1, PNPNodeData::NUM_VARS);

      //Boundary indices start from 1. Hence the first row is left unused.

      const int c1Index = 0;
      BoundaryConditionArray(LEFT, c1Index) = NEUMANN;
      BoundaryConditionArray(RIGHT, c1Index) = NEUMANN;
      BoundaryConditionArray(BOTTOM, c1Index) = NEUMANN;
      BoundaryConditionArray(TOP, c1Index) = NEUMANN;
      BoundaryConditionArray(BACK, c1Index) = NEUMANN;
      BoundaryConditionArray(FRONT, c1Index) = NEUMANN;

      //1 - C2 Boundary Imposition
      const int c2Index = 1;
      BoundaryConditionArray(LEFT, c2Index) = NEUMANN;
      BoundaryConditionArray(RIGHT, c2Index) = NEUMANN;
      BoundaryConditionArray(BOTTOM, c2Index) = NEUMANN;
      BoundaryConditionArray(TOP, c2Index) = NEUMANN;
      BoundaryConditionArray(BACK, c2Index) = NEUMANN;
      BoundaryConditionArray(FRONT, c2Index) = NEUMANN;

      //2 - Phi Boundary Imposition
      const int phiIndex = PNPNodeData::PHI_IDX;
      BoundaryConditionArray(LEFT, phiIndex) = DIRICHLET;
      BoundaryConditionArray(RIGHT, phiIndex) = DIRICHLET;
      BoundaryConditionArray(BOTTOM, phiIndex) = DIRICHLET;
      BoundaryConditionArray(TOP, phiIndex) = DIRICHLET;
      BoundaryConditionArray(BACK, phiIndex) = DIRICHLET;
      BoundaryConditionArray(FRONT, phiIndex) = DIRICHLET;

  }
};

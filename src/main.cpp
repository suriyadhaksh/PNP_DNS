
#include <talyfem/talyfem.h>
#include <math.h>
#include <fstream>
#include <string>

#include <PNPInputData.h>
#include <PNPNodeData.h>
#include <PNPGridField.h>
//#include <derived_example.h>
#include <PNPManufacturedSoln.h>
#include <UtilFunctions.h>
using namespace TALYFEMLIB;
static char help[] = "Solves the Poisson-Nernst-Planck equations!";


inline bool SetIC(PNPGridField& data, PNPInputData& idata) {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  switch (idata.typeOfIC) {
    case 1:
      data.SetIC(idata.nsd);
      return true;
    case 0:
      try {
        load_gf(&data, &idata);
        data.UpdateDataStructures();
      } catch(const TALYException& e) {
        e.print();
        PrintWarning("Failed to load GridField data!");
        return false;
      }
      return true;
    default:
      if (rank == 0) std::cerr << "IC not set up " << std::endl;
      return false;
  }
}


int main(int argc, char **args) {
  PetscInitialize(&argc, &args, NULL, help);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  PNPInputData input_data;
  GRID* p_grid = NULL;
  {
    PNPGridField data(&input_data);

    // Read input data from file "config.txt"
    if (!input_data.ReadFromFile()) {
      throw TALYException() << "Error reading config.txt";
    }

    // Loging inputdata to std::clog or file
    if (rank == 0) {
      std::clog << input_data;

      std::fstream filePtErrorManufacSol;
      std::string fileName = "Error_data_" + std::to_string(PNPNodeData::NO_OF_SPECIES) + "_Species.plt";
      filePtErrorManufacSol.open(fileName, std::ios::app);
      filePtErrorManufacSol << "VARIABLES = \"t\"";

      for (int i = 0; i < PNPNodeData::NO_OF_SPECIES; i++) {
          filePtErrorManufacSol << " \"C_" << i << "\"";
      }

      filePtErrorManufacSol << " \"Phi\"" << std::endl;
      filePtErrorManufacSol.close();
    }

    //  Check if inputdata is complete
    if (!input_data.CheckInputData()) {
      throw TALYException() << "[ERR] Problem with input data, check the config file!";
    }

    // Based on inputdata create Grid
    CreateGrid(p_grid, &input_data);

    // check gaussian quadrature
    /*FEMElm fe(p_grid, BASIS_ALL);
    fe.refill(0, input_data.basisRelativeOrder);
    std::cout << fe.n_itg_pts() << "\n";
    exit(0);*/

    // Construct gridfield based on Grid
    data.redimGrid(p_grid);
    data.redimNodeData();

    // Set Initial Conditions
    if (!SetIC(data, input_data)) {
      PrintResults("Failed to set initial conditions", false,
                   input_data.shouldFail);
      delete p_grid;
      throw TALYException() << "Problem with IC, not loaded";
    }
    PrintStatus("IC set");

    // Set Solver parameters
    const int nOfDofPerNode = PNPNodeData::NO_OF_SPECIES + 1;  // number of degree of freedom per node
    PNPManufacturedSoln PNPEq2(input_data.nsd);

    PNPEq2.setParams(input_data.lambda , input_data.L[0]/input_data.Nelem[0], input_data.z);
    PNPEq2.copyBoundaryConditions(input_data.BoundaryConditionArray);


    PNPEq2.redimSolver(p_grid, nOfDofPerNode, false, input_data.basisRelativeOrder);

    PNPEq2.setData(&data);
    PrintStatus("solver initialized!");

    PNPEq2.fillEssBC();
    PrintStatus("Dirichlet BC applied");
    data.setMMS(&PNPEq2);

    save_gf(&data, &input_data, "data_guess.plt", 0.0);

    int n_time_steps = int( input_data.totalT / input_data.dt );
    int no_of_frames = 70;
    int time_skip = n_time_steps / no_of_frames;

    double t = 0.0;
    int timeStepCounter = 0;
    double dt = input_data.dt;
    while (t < input_data.totalT) {

      PrintStatus("Solver time step:", timeStepCounter);

      PNPEq2.Solve(dt, t);

      // copy u into u_prev
      data.UpdateDataStructures();

      PNPEq2.fillEssBC();

      data.setMMSGridField(t+dt);

      t += dt;
      timeStepCounter++;

      if (timeStepCounter%time_skip == 0) {
          std::string name = Suffix("data.plt", timeStepCounter);
          save_gf(&data, &input_data, (char *) name.c_str(), t);
      }

      PrintStatus("Time Step: ", timeStepCounter, " Time: ", t);
      data.PrintError(t);


      //std::string name = Suffix("data.plt", t / dt);
      //save_gf(&data, &input_data, (char *) name.c_str(), t);
      //PrintInfo("time t = ", t);
//        data.PrintError(t);

    }

    PrintStatus("Finished Solve!");

    std::string filename = "data_final.plt";
    save_gf(&data, &input_data, filename.c_str(), t);
    PrintStatus("Final data sent to file.");


/*    if (input_data.ifPrintPltFiles) {
      save_gf(&data, &input_data, "data.plt", 0.0);
    }*/
  }

  // clean up
  DestroyGrid(p_grid);

  // PrintResults("Completed successfully.", true, input_data.shouldFail);

  PetscFinalize();
  return 0;
}

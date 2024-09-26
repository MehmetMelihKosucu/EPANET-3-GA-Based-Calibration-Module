/* EPANET 3 Genetic Algorithm Based Calibration Module
 *
 * Copyright (c) 2016 Open Water Analytics
 * Distributed under the MIT License (see the LICENSE file for details).
 *
 */

#include <string>
#include "../epanet3.h"
//#include "..\\Calibration\\GeneticAlgorithm.h" // For Windows
#include "../Calibration/GeneticAlgorithm.h" // For Linux
//#include "..\\Calibration\\Chromosome.h" // For Windows
#include "../Calibration/Chromosome.h" // For Linux
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <random>
#include <bitset>
#include <vector>
#include <cmath>
#include <numeric>
#include <time.h>
#include <iomanip>
#include <time.h>
#include <sstream>
#include "omp.h"

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cout << "\nCorrect syntax is: epanet3 inpFile rptFile (outFile)\n";
        return 1;
    }
    int err = 0;

    clock_t start_t = clock();

    //... check number of command line arguments
    if (argc < 3)
    {
        std::cout << "\nCorrect syntax is: epanet3 inpFile rptFile (outFile)\n";
        return 0;
    }

    //... retrieve file names from command line
    const char* inpFile = argv[1];
    const char* rptFile = argv[2];
    const char* outFile = "";
    if (argc > 3) outFile = argv[3];

    EN_Project project = EN_createProject();

    
    int reportFileStatus = EN_openReportFile(rptFile, project);
    if (reportFileStatus != 0) {
        std::cerr << "Failed to open the EPANET report file.\n";
        EN_deleteProject(project);
        return reportFileStatus;
    }

    int loadStatus = EN_loadProject(inpFile, project);
    if (loadStatus != 0) {
        std::cerr << "Failed to load the EPANET project file.\n";
        EN_deleteProject(project);
        return loadStatus;
    }

    int openOutputFileStatus = EN_openOutputFile(outFile, project);
    if (openOutputFileStatus != 0) {
        std::cerr << "Failed to open the EPANET output file.\n";
        EN_deleteProject(project);
        return openOutputFileStatus;
    }

    //int runStatus = EN_runProject(project, "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Simulated.dat"); // For Windows
    EN_runProject(project, "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated.dat"); // For Linux
    
    double mutationRate = 0.1;
    double crossoverRate = 0.9;
    int populationSize = 100;
    int generations = 250;

    GA::GeneticAlgorithm ga(mutationRate, crossoverRate, populationSize, generations, project);

     // Initialize population
    int initStatus = ga.initializePopulation();
    if (initStatus != 0) {
        std::cerr << "Failed to initialize population with status: " << initStatus << std::endl;
        EN_deleteProject(project);
        return initStatus;
    }

    ga.run();

    // Re-initialize the project for applying the best solution
    EN_deleteProject(project);
    project = EN_createProject();

    // Reload the original project files
    if (EN_openReportFile(rptFile, project) != 0 ||
        EN_loadProject(inpFile, project) != 0 ||
        EN_openOutputFile(outFile, project) != 0) {
        std::cerr << "Failed to re-initialize the EPANET project for best solution.\n";
        EN_deleteProject(project);
        return 1;
    }

    Chromosome bestSolution = ga.getBestSolution();
    bestSolution.applyChromosomeToProject(project);

    /* EN_Project bestProject = EN_createProject();
    EN_openReportFile(rptFile, bestProject);
    EN_loadProject(inpFile, bestProject);
    EN_openOutputFile(outFile, bestProject); // */

    //int runStatus = EN_runProject(project, "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Simulated_best.dat"); // For Windows
    if (EN_runProject(project, "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated_best.dat") != 0) { // For Linux
        std::cerr << "Failed to run EPANET simulation.\n";
        EN_deleteProject(project);
        return 1;
    }
    // */
    if (!err) {
        clock_t end_t = clock();
        double cpu_t = ((double)(end_t - start_t)) / CLOCKS_PER_SEC;
        std::stringstream ss;
        ss << "\n  Simulation completed in ";
        ss.str("");
        if (cpu_t < 0.001) ss << "< 0.001 sec.";
        else ss << std::setprecision(3) << cpu_t << " sec.";

        std::cout << "\n Writing report ...                           ";
        if (EN_writeReport(project) != 0 || (outFile && EN_saveOutput(project) != 0)) {
            std::cerr << "Failed to save reports or output files.\n";
            EN_deleteProject(project);
            return 1;
        }
        std::cout << "\n    Simulation completed.                         \n";
        std::cout << "\n... EPANET completed in " << ss.str() << "\n";
    }

    EN_deleteProject(project);
    //EN_deleteProject(bestProject);
    return 0;
}   
    

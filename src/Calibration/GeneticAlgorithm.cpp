/* EPANET 3 Genetic Algorithm Based Calibration Module
 *
 *
 */
#include "GeneticAlgorithm.h"
#include "Chromosome.h"
#include <algorithm>
#include <random>
#include <iostream>
#include <numeric>
#include <functional>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <bitset>
#include <time.h>
#include <unordered_map>
#include "epanet3.h"
#include <map>
#include <fstream>
#include <filesystem>
#include <mutex>
#include <ctime>
#include "omp.h"
#include "utilities.h"
#include <cstdio>
#include <deque>



namespace GA {

    // Constructor
    GeneticAlgorithm::GeneticAlgorithm(double mutationRate, double crossoverRate, size_t populationSize, size_t generations, EN_Project& project)
        : mutationRate(mutationRate), crossoverRate(crossoverRate), populationSize(populationSize), generations(generations), project(project),
        bestOverall(populationSize), bestOverallFitness(std::numeric_limits<double>::max()) {
        population.reserve(populationSize);
    }

    int GeneticAlgorithm::initializePopulation() {
        int numPipes;
        EN_Project baseProject = EN_createProject();

        int returnStatus = 0;

        //int reportFileStatus = EN_openReportFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.rpt", baseProject); // For Windows
        int reportFileStatus = EN_openReportFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.rpt", baseProject); // For Linux
        if (reportFileStatus != 0) {
            std::cerr << "Failed to open the EPANET report file.\n";
            EN_deleteProject(baseProject);
            return reportFileStatus;
        }

        //int loadStatus = EN_loadProject("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.inp", baseProject); // For Windows
        if (EN_loadProject("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.inp", baseProject) != 0) { // For Linux
            std::cerr << "Failed to load the EPANET project file." << std::endl;
            EN_deleteProject(baseProject);
            return -1;
        }

        //int openOutputFileStatus = EN_openOutputFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.out", baseProject); // For Windows
        int openOutputFileStatus = EN_openOutputFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.out", project); // For Linux
        if (openOutputFileStatus != 0) {
            std::cerr << "Failed to open the EPANET output file.\n";
            EN_deleteProject(baseProject);
            return openOutputFileStatus;
        }

        if (EN_getCount(EN_PIPECOUNT, &numPipes, baseProject) != 0) {
            std::cerr << "Error getting pipe count from EPANET." << std::endl;
            EN_deleteProject(baseProject);
            return -1;
        }

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(5, 140);
        
        population.clear();
        population.reserve(populationSize);

        for (int i = 0; i < populationSize; ++i) {

            Chromosome individual(numPipes);
            for (int j = 0; j < numPipes; ++j) {
                double roughness = dis(gen);
                std::string binaryRoughness = Chromosome::doubleToBinaryString(roughness, 5, 140, 16);
                individual.setGene(j, binaryRoughness);
            }

            EN_Project clonedProject = EN_createProject();

            // Copy the Simulated.dat file to a new file for each individual
            //std::string simFileCopy = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Simulated_Ind_" + std::to_string(i) + ".dat"; // For Windows
            std::string simFileCopy = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated_Ind_" + std::to_string(i) + ".dat"; // For Linux

            //int reportFileStatus = EN_openReportFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.rpt", clonedProject); // For Windows
            int reportFileStatus = EN_openReportFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.rpt", clonedProject); // For Linux
            if (reportFileStatus != 0) {

                std::cerr << "Failed to open the EPANET report file.\n";
                EN_deleteProject(clonedProject);
                returnStatus = reportFileStatus;
            }
            //int loadStatus = EN_loadProject("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.inp", clonedProject); // For Windows
            if (EN_loadProject("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.inp", clonedProject) != 0) { // For Linux
       
                std::cerr << "Failed to load the EPANET project file." << std::endl;
                EN_deleteProject(clonedProject);
                returnStatus = -1;
            }
            //int openOutputFileStatus = EN_openOutputFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.out", clonedProject); // For Windows
            int openOutputFileStatus = EN_openOutputFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.out", clonedProject); // For Linux
            if (openOutputFileStatus != 0) {
 
                std::cerr << "Failed to open the EPANET output file.\n";
                EN_deleteProject(clonedProject);
                returnStatus = openOutputFileStatus;
            }

            individual.applyChromosomeToProject(clonedProject);
            if (EN_runProject(clonedProject, simFileCopy.c_str()) != 0) {
 
                std::cerr << "Failed to run EPANET simulation for individual " << i << std::endl;
                EN_deleteProject(project);
                continue;
            }

             // Retrieve simulated values from the EPANET project
            auto simulatedValues = Chromosome::readSimulatedValues(simFileCopy.c_str());

            // Compute fitness after the EPANET simulation
            computeFitnessAfterRun(clonedProject, individual, i, 0, simulatedValues);

            std::cout << "Fitness of individual " << i << ": " << individual.getFitness() << "\n";
            population.push_back(individual);
            EN_deleteProject(clonedProject);  // Clean up the project after processing each individual
        }

        EN_deleteProject(baseProject);

        std::cout << "Population initialized with " << population.size() << " individuals." << std::endl;

        return returnStatus;
    }

    void GeneticAlgorithm::run() {
        //int generations = 100; // Example value
        if (population.empty()) {
            throw std::runtime_error("Population not initialized.");
        }
        //std::ofstream outFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\fitness_values.dat"); // For Windows
        std::ofstream outFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/fitness_values.dat"); // For Linux
        if (!outFile.is_open()) {
            throw std::runtime_error("Could not open file for writing.");
        }

       for (size_t gen = 0; gen < static_cast<size_t>(generations); ++gen) { // Updated to size_t

            std::cout << "Generation " << gen << "\n";

            evolve();

            Chromosome best = getBestSolution();
            std::cout << "Best Fitness in Generation " << gen << ": " << best.getFitness() << "\n";
            outFile << "Generation " << gen << ": " << best.getFitness() << "\n";

            /* if (!roughnessFile.is_open()) {
                throw std::runtime_error("Could not open file for writing.");
            } // */

            for (size_t i = 0; i < population.size(); ++i) {
                const auto& individual = population[i];
                //roughnessFile << "Individual " << i << ": Fitness: " << individual.getFitness() << ", Roughness Values: ";
                for (const auto& gene : individual.getGenes()) {
                    double realValue = population[i].decodeBinaryToReal(gene, 5, 140);
                }
            }
           
            //generateBestSolutionSimulatedFile("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\best_simulated", gen); // For Windows
            generateBestSolutionSimulatedFile("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/best_simulated", gen); // For Linux

            if (best.getFitness() < bestOverallFitness) {
                bestOverallFitness = best.getFitness();
                bestOverall = best;
            }

            std::cout << "Best Overall Fitness: " << bestOverall.getFitness() << "\n";
            outFile << "Best Overall Fitness: " << bestOverall.getFitness() << "\n";
        }

        outFile.close();
    }

    void GeneticAlgorithm::generateBestSolutionSimulatedFile(const std::string& baseFilename, int generationNumber) {
        Chromosome bestSolution = getBestSolution();
        EN_Project bestProject = EN_createProject();

        //const char* inpFile = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.inp"; // For Windows
        const char* inpFile = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.inp"; // For Linux
        //const char* rptFile = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.rpt"; // For Windows
        const char* rptFile = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.rpt"; // For Linux
        //const char* outFile = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\c-town_true_network-Uncalibrated-EPA3.out"; // For Windows
        const char* outFile = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.out"; // For Linux
        //std::string baseFilename = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated";

        int reportFileStatus = EN_openReportFile(rptFile, bestProject);
        if (reportFileStatus != 0) {
            std::cerr << "Failed to open the EPANET report file.\n";
            EN_deleteProject(bestProject);
            return;
        }

        int loadStatus = EN_loadProject(inpFile, bestProject);
        if (loadStatus != 0) {
            std::cerr << "Failed to load the EPANET project file.\n";
            EN_deleteProject(bestProject);
            return;
        }

        int openOutputFileStatus = EN_openOutputFile(outFile, bestProject);
        if (openOutputFileStatus != 0) {
            std::cerr << "Failed to open the EPANET output file.\n";
            EN_deleteProject(bestProject);
            return;
        }

        bestSolution.applyChromosomeToProject(bestProject);

        if (EN_runProject(bestProject, baseFilename.c_str()) != 0) {
            std::cerr << "Failed to run EPANET simulation.\n";
            EN_deleteProject(bestProject);
            return;
        }

        // Construct the new filename with the generation number
        std::string newFilename = baseFilename + "_gen_" + std::to_string(generationNumber) + ".dat";

        // Copy the Simulated.dat file to the new file using standard file I/O
        //std::string sourceFile = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Simulated.dat"; // For Windows
        std::string sourceFile = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated.dat"; // For Linux
        std::ifstream src(sourceFile, std::ios::binary);
        std::ofstream dst(newFilename, std::ios::binary);

        if (!src.is_open()) {
            std::cerr << "Failed to open source file: " << sourceFile << "\n";
            EN_deleteProject(bestProject);
            return;
        }

        if (!dst.is_open()) {
            std::cerr << "Failed to open destination file: " << newFilename << "\n";
            EN_deleteProject(bestProject);
            return;
        }

        dst << src.rdbuf();

        src.close();
        dst.close();

        EN_deleteProject(bestProject);
    }

    std::vector<double> GeneticAlgorithm::readRoughnessValuesFromInputFile(const std::string& filePath) {
        std::ifstream inputFile(filePath);
        if (!inputFile) {
            throw std::runtime_error("Failed to open the input file: " + filePath);
        }

        std::vector<double> roughnessValues;
        std::string line;
        bool inPipesSection = false;

        while (std::getline(inputFile, line)) {
            if (line.find("[PIPES]") != std::string::npos) {
                inPipesSection = true;
                continue;
            }
            if (inPipesSection && (line.empty() || line[0] == ';' || line.find("[") != std::string::npos)) {
                if (line.find("[") != std::string::npos && line.find("[PIPES]") == std::string::npos) {
                    inPipesSection = false;
                }
                continue;
            }
            if (inPipesSection && !line.empty() && line[0] != ';') {
                std::istringstream lineStream(line);
                std::string id, node1, node2, length, diameter, roughness, minorLoss, status;
                if (lineStream >> id >> node1 >> node2 >> length >> diameter >> roughness >> minorLoss >> status) {
                    roughnessValues.push_back(std::stod(roughness));
                }
            }
        }
        inputFile.close();
        return roughnessValues;
    }

    Chromosome GeneticAlgorithm::selectParent() {
        // Calculate the total "fitness" for selection
        double totalFitness = std::accumulate(population.begin(), population.end(), 0.0,
            [](double sum, const Chromosome& chrom) {
                return sum + 1.0 / (chrom.getFitness() + 1.0);
            });

        // Randomly select a point in the total "fitness"
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, totalFitness);
        double value = dis(gen);

        // Select the parent based on the accumulated "fitness"
        double runningTotal = 0;
        for (const Chromosome& chrom : population) {
            runningTotal += 1.0 / (chrom.getFitness() + 1.0);
            if (runningTotal >= value) {
                return chrom;
            }
        }
        // If for some reason we do not return inside the loop, return the last chromosome
        return population.back();
    }

    void GeneticAlgorithm::evolve() {
        size_t eliteCount = 2;
        size_t mu = population.size();
        size_t lambda = mu * 2;

        // Sort the population by fitness
        std::stable_sort(population.begin(), population.end(), [](const Chromosome &a, const Chromosome &b) {
            return a.getFitness() < b.getFitness();
        });

        Chromosome bestInCurrentGeneration = population[0];

        if (bestInCurrentGeneration.getFitness() < bestOverallFitness) {
            bestOverallFitness = bestInCurrentGeneration.getFitness();
            bestOverall = bestInCurrentGeneration;
        }

        std::vector<Chromosome> selectedPopulation;
        for (size_t i = 0; i < population.size(); ++i) {
            selectedPopulation.push_back(rouletteWheelSelection());
        }

        std::vector<Chromosome> newPopulation;
        newPopulation.reserve(lambda);

        // Preserve elites
        for (size_t i = 0; i < eliteCount; ++i) {
            Chromosome elite(selectedPopulation[i].getGenes().size());
            elite.copySolution(selectedPopulation[i]);
            newPopulation.push_back(elite);
        }

        // Ensure the overall best individual is included
        newPopulation.push_back(bestOverall);

        double generationFactor = 1.0 - (static_cast<double>(generationCount) / generations);

        while (newPopulation.size() < lambda) {
            Chromosome parent1 = rouletteWheelSelection();
            Chromosome parent2 = rouletteWheelSelection();
            Chromosome parent3 = rouletteWheelSelection();
            Chromosome parent4 = rouletteWheelSelection();

            Chromosome offspring = differentialEvolution(parent1, parent2, parent3, parent4);
            offspring.mutate(mutationRate, generationFactor);

            newPopulation.push_back(std::move(offspring));
        }

        evaluateFitnessInParallel(newPopulation);

        double diversity = calculateDiversity(population);
        double diversityThreshold = 3;

        if (diversity < diversityThreshold) {
            // Aggressively increase mutation rate to restore diversity
            mutationRate = std::min(mutationRate * 3, maxMutationRate);
            std::cout << "Increasing mutation rate due to low diversity after selection.\n";

            // Optionally apply mutation to the entire population to introduce diversity
            for (auto& individual : population) {
                individual.mutate(mutationRate, generationFactor);
            }
        } else {
            // Gradually adjust mutation rate downwards if diversity is sufficient
            mutationRate = std::max(mutationRate * 0.9, minMutationRate);
        }

        // Ensure mutation rate stays within defined bounds
        mutationRate = std::max(minMutationRate, std::min(mutationRate, maxMutationRate));

        population = selectFittestMuLambda(population, newPopulation);

        // Debugging output to track diversity
        std::cout << "Population diversity: " << diversity << std::endl;

        // Write roughness and fitness values to a file for the current generation
        writeRoughnessAndFitnessValues(newPopulation, generationCount);

        // Increment the generation count
        generationCount++;
    }

    void GeneticAlgorithm::randomRestart() {
        // Sort the population by fitness in descending order (worst-fit individuals at the end)
        std::sort(population.begin(), population.end(), [](const Chromosome &a, const Chromosome &b) {
            return a.getFitness() > b.getFitness(); // Sort by descending fitness
        });

        // Determine how much of the population to reinitialize
        size_t restartCount = population.size() / 10; // Restart 10% of the population

        // Get the correct gene size from a good individual (i.e., the first individual in the sorted population)
        size_t correctGeneSize = population.front().getGenes().size();
        size_t correctGeneLength = population.front().getGenes().front().size();

        std::cout << "Correct gene size for reinitialization: " << correctGeneSize << "\n";
        std::cout << "Correct gene length for reinitialization: " << correctGeneLength << "\n";

        // Reinitialize the worst-fit individuals
        for (size_t i = 0; i < restartCount; ++i) {
            size_t worstIndex = population.size() - 1 - i; // Indices of the worst individuals

            // Initialize new Chromosome with the correct gene size and length
            Chromosome newChromosome(correctGeneSize);
            for (size_t j = 0; j < correctGeneSize; ++j) {
                newChromosome.setGene(j, Chromosome::generateBinaryString(correctGeneLength));
            }

            population[worstIndex] = std::move(newChromosome);

            // Debugging output to ensure proper reinitialization
            std::cout << "Added new individual with genes of size " << correctGeneSize << " and gene length " << correctGeneLength << ".\n";
        }

        std::cout << "Random restart applied by replacing the worst " << restartCount << " individuals with new ones.\n";
    }

    bool GeneticAlgorithm::noImprovementForNGenerations(int n) {
        static int generationsWithoutImprovement = 0;
        static double bestFitness = std::numeric_limits<double>::max();

        double currentBestFitness = population[0].getFitness();
        if (currentBestFitness < bestFitness) {
            bestFitness = currentBestFitness;
            generationsWithoutImprovement = 0;
        } else {
            generationsWithoutImprovement++;
        }

        return generationsWithoutImprovement >= n;
    }

    double GeneticAlgorithm::calculateTotalFitness() const {
        double totalFitness = 0.0;
        for (const auto& chrom : population) {
            totalFitness += 1.0 / (chrom.getFitness() + 1.0); // Normalize fitness
        }
        return totalFitness;
    }

    std::vector<double> GeneticAlgorithm::computeSelectionProbabilities(double totalFitness) const {
        std::vector<double> selectionProbabilities;
        for (const auto& chrom : population) {
            selectionProbabilities.push_back((1.0 / (chrom.getFitness() + 1.0)) / totalFitness); // Normalize fitness
        }
        return selectionProbabilities;
    }

    std::vector<double> GeneticAlgorithm::computeCumulativeProbabilities(const std::vector<double>& selectionProbabilities) const {
        std::vector<double> cumulativeProbabilities;
        double cumulativeSum = 0.0;
        for (double probability : selectionProbabilities) {
            cumulativeSum += probability;
            cumulativeProbabilities.push_back(cumulativeSum);
        }
        return cumulativeProbabilities;
    }

    Chromosome GeneticAlgorithm::rouletteWheelSelection() const {
        double totalFitness = calculateTotalFitness();
        std::vector<double> selectionProbabilities = computeSelectionProbabilities(totalFitness);
        std::vector<double> cumulativeProbabilities = computeCumulativeProbabilities(selectionProbabilities);

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        double randomValue = dis(gen);

        for (size_t i = 0; i < cumulativeProbabilities.size(); ++i) {
            if (randomValue <= cumulativeProbabilities[i]) {
                return population[i];
            }
        }
        return population.back();
    }

    void GeneticAlgorithm::writeRoughnessAndFitnessValues(const std::vector<Chromosome>& population, int generation) {
        //std::string filename = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\roughness_fitness_gen_" + std::to_string(generation) + ".dat"; // For Windows
        std::string filename = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/roughness_fitness_gen_" + std::to_string(generation) + ".dat"; // For Linux
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Failed to open " << filename << " for writing.\n";
            return;
        }

        for (size_t i = 0; i < population.size(); ++i) {
            const auto& chrom = population[i];
            outFile << "Individual " << i << " Fitness: " << chrom.getFitness() << "\n";
            outFile << "Roughness values: ";
            for (const auto& gene : chrom.getGenes()) {
                double roughness = chrom.decodeBinaryToReal(gene, 5.0, 140.0);
                outFile << roughness << " ";
            }
            outFile << "\n";
        }

        outFile.close();
    }

    // Function to calculate diversity of the population
    double GeneticAlgorithm::calculateDiversity(const std::vector<Chromosome>& population) const {
        double diversity = 0.0;
        for (size_t i = 0; i < population.size(); ++i) {
            for (size_t j = i + 1; j < population.size(); ++j) {
                diversity += population[i].distance(population[j]);
            }
        }
        diversity /= (population.size() * (population.size() - 1)) / 2.0;
        return diversity;
    }

    // Function for Binary Tournament Selection
    Chromosome GeneticAlgorithm::binaryTournamentSelection() {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, population.size() - 1);

        // Select two random candidates
        Chromosome candidate1 = population[dis(gen)];
        Chromosome candidate2 = population[dis(gen)];

        // Introduce a higher selection pressure
        double candidate1Fitness = candidate1.getFitness();
        double candidate2Fitness = candidate2.getFitness();

        // If fitness is equal, use dynamic selection based on the generation count
        if (candidate1Fitness == candidate2Fitness) {
            std::uniform_real_distribution<> disReal(0.0, 1.0);
            double dynamicProbability = 0.7 + 0.3 * (generationCount / static_cast<double>(generations)); // Increasing bias over time
            return (disReal(gen) <= dynamicProbability) ? candidate1 : candidate2;
        }

        // Otherwise, return the better candidate with a higher probability
        return (candidate1Fitness < candidate2Fitness) ? candidate1 : candidate2;
    }

    // Function to select the fittest individuals (Mu-Lambda Strategy)
    std::vector<Chromosome> GeneticAlgorithm::selectFittestMuLambda(const std::vector<Chromosome>& parents, const std::vector<Chromosome>& offspring) {
        // Combine parent and offspring populations
        std::vector<Chromosome> combinedPopulation;
        combinedPopulation.reserve(parents.size() + offspring.size());
        combinedPopulation.insert(combinedPopulation.end(), parents.begin(), parents.end());
        combinedPopulation.insert(combinedPopulation.end(), offspring.begin(), offspring.end());

        // Sort combined population in ascending order of fitness
        std::sort(combinedPopulation.begin(), combinedPopulation.end(), [](const Chromosome& a, const Chromosome& b) {
            return a.getFitness() < b.getFitness();
            });

        // Ensure we do not exceed the size of the combined population
        size_t newSize = std::min(parents.size(), combinedPopulation.size());
        std::vector<Chromosome> fittestPopulation(combinedPopulation.begin(), combinedPopulation.begin() + newSize);

        return fittestPopulation;
    }

    Chromosome GeneticAlgorithm::getBestSolution() const {
        return *std::min_element(population.begin(), population.end(),
            [](const Chromosome& a, const Chromosome& b) {
                //return a.getFitness() > b.getFitness(); fitnessın büyük olması daha iyi olur
				return a.getFitness() < b.getFitness(); // fitnessın küçük olması daha iyi olur
            });
    }

    void GeneticAlgorithm::evaluateFitnessInParallel(std::vector<Chromosome>& population) {
        int numThreads = 20;
        omp_set_num_threads(numThreads);
        double bestFitness = std::numeric_limits<double>::max();
        std::string bestInputFile;  // Track the path to the best individual's input file

        //std::string base_path = "E:\\EPANET-GA-Parallel\\Networks\\C-Town\\"; // For Windows
        std::string base_path = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/"; // For Linux

        #pragma omp parallel
        {
            double localBestFitness = std::numeric_limits<double>::max();
            std::string localBestInputFile;
            std::string threadInputFile, threadOutputFile, threadReportFile, threadSimFile;

            #pragma omp for schedule(static, 200)
            for (size_t i = 0; i < population.size(); ++i) {
                int thread_id = omp_get_thread_num();
                bool simulationSuccess = true;

                // Generate unique file names for each thread
                threadInputFile = base_path + "c-town_true_network-Uncalibrated-EPA3_thread_" + std::to_string(thread_id) + "_" + std::to_string(i) + ".inp";
                threadOutputFile = base_path + "c-town_true_network-Uncalibrated-EPA3_thread_" + std::to_string(thread_id) + "_" + std::to_string(i) + ".dat";
                threadReportFile = base_path + "c-town_true_network-Uncalibrated-EPA3_thread_" + std::to_string(thread_id) + "_" + std::to_string(i) + ".rpt";
                threadSimFile = base_path + "Simulated_thread_" + std::to_string(thread_id) + "_" + std::to_string(i) + ".dat";

                // Copy base input file to the thread-specific input file
                #pragma omp critical
                {
                    if (!Utilities::copyFile(base_path + "c-town_true_network-Uncalibrated-EPA3.inp", threadInputFile)) {
                        std::cerr << "Error: Unable to copy input file to " << threadInputFile << std::endl;
                        simulationSuccess = false;
                    }
                }

                if (!simulationSuccess) {
                    continue; // Skip to the next individual if file copying fails
                }

                EN_Project project = EN_createProject();
                if (!project) {
                    std::cerr << "Failed to create EPANET project for thread " << thread_id << std::endl;
                    simulationSuccess = false;
                }

                if (simulationSuccess) {
                    int reportFileStatus = EN_openReportFile(threadReportFile.c_str(), project);
                    if (reportFileStatus != 0) {
                        std::cerr << "Failed to open the EPANET report file for thread " << thread_id << std::endl;
                        EN_deleteProject(project);
                        simulationSuccess = false;
                    }
                }

                if (simulationSuccess) {
                    int loadStatus;
                    #pragma omp critical
                    {
                        loadStatus = EN_loadProject(threadInputFile.c_str(), project);
                    }
                    if (loadStatus != 0) {
                        std::cerr << "Failed to load the EPANET project for thread " << thread_id << std::endl;
                        EN_deleteProject(project);
                        simulationSuccess = false;
                    }
                }

                if (simulationSuccess) {
                    int openOutputFileStatus = EN_openOutputFile(threadOutputFile.c_str(), project);
                    if (openOutputFileStatus != 0) {
                        std::cerr << "Failed to open the EPANET output file for thread " << thread_id << std::endl;
                        EN_deleteProject(project);
                        simulationSuccess = false;
                    }
                }

                if (simulationSuccess) {
                    #pragma omp critical
                    {
                        std::vector<double> appliedRoughnessValues = population[i].applyChromosomeToProject(project);
                        int runStatus = EN_runProject(project, threadSimFile.c_str());
                        if (runStatus != 0) {
                            std::cerr << "EPANET simulation failed for individual " << i << " with error code: " << runStatus << std::endl;
                            simulationSuccess = false;
                        }
                    }
                }

                if (simulationSuccess) {
                    double fitness;
                    #pragma omp critical
                    {
                        auto simulatedValues = Chromosome::readSimulatedValues(threadSimFile.c_str());
                        fitness = computeFitnessAfterRun(project, population[i], i, 0, simulatedValues);
                    }

                    if (fitness < localBestFitness) {
                        localBestFitness = fitness;
                        localBestInputFile = threadInputFile;
                    }
                }

                EN_deleteProject(project);
                
                #pragma omp critical
                {
                    if (localBestFitness < bestFitness) {
                        bestFitness = localBestFitness;
                        bestInputFile = localBestInputFile;  // Track the correct file path
                    }
                }

                // Clean up thread-specific files
                std::remove(threadOutputFile.c_str());
                std::remove(threadReportFile.c_str());
                std::remove(threadSimFile.c_str());
            }
        }

        // After the parallel block, copy the best individual's input file to the final file
        if (!bestInputFile.empty()) {
            std::string finalBestFile = base_path + "c-town_true_network-Uncalibrated-EPA3-best.inp";
            
            if (Utilities::copyFile(bestInputFile, finalBestFile)) {
                std::cout << "Best individual with fitness " << bestFitness << " copied to " << finalBestFile << std::endl;
                std::remove(bestInputFile.c_str());  // Remove the best thread-specific input file after copying
            } else {
                std::cerr << "Error: Unable to copy the best input file from " << bestInputFile << " to " << finalBestFile << std::endl;
            }
        }
    }

    double GeneticAlgorithm::computeFitnessAfterRun(EN_Project project, Chromosome& chromosome, int chromosomeIndex, int generationCount, const std::map<std::string, std::vector<double>>& simulatedValues) {
        
        // Read observed values
        //auto observedValues = Chromosome::readObservedValues("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Observed.dat"); // For Windows
        auto observedValues = Chromosome::readObservedValues("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Observed.dat"); // For Linux
        
        double sumSquaredErrors = 0.0;

        // Log for the current chromosome
        //std::ofstream detailedLog("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/detailed_fitness_log_" + std::to_string(generationCount) + "_" + std::to_string(chromosomeIndex) + ".txt");

        // Calculate max and min for normalization across all parameters
        std::vector<double> maxObserved(observedValues.begin()->second.size(), std::numeric_limits<double>::min());
        std::vector<double> minObserved(observedValues.begin()->second.size(), std::numeric_limits<double>::max());

        for (const auto& obsPair : observedValues) {
            const std::vector<double>& obsData = obsPair.second;
            for (size_t i = 0; i < obsData.size(); ++i) {
                maxObserved[i] = std::max(maxObserved[i], obsData[i]);
                minObserved[i] = std::min(minObserved[i], obsData[i]);
            }
        }

        // Calculate NSSE - Normalized Sum of Squared Errors
        for (const auto& obsPair : observedValues) {
            const std::string& id = obsPair.first;
            const std::vector<double>& obsData = obsPair.second;

            if (simulatedValues.count(id) == 0) continue; // Skip if no simulated data for this id

            const std::vector<double>& simData = simulatedValues.at(id);
            for (size_t t = 0; t < obsData.size() && t < simData.size(); ++t) {
                double observedValue = obsData[t];
                double simulatedValue = simData[t];

                // Normalize the error using max and min
                double normalizedObservedValue = (maxObserved[t] - minObserved[t] != 0) ? (observedValue - minObserved[t]) / (maxObserved[t] - minObserved[t]) : 0;
                double normalizedSimulatedValue = (maxObserved[t] - minObserved[t] != 0) ? (simulatedValue - minObserved[t]) / (maxObserved[t] - minObserved[t]) : 0;

                double error = std::abs(normalizedSimulatedValue - normalizedObservedValue);
                double squaredError = error * error;

                sumSquaredErrors += squaredError;

                // Log detailed information
                /*detailedLog << "ID: " << id << ", Time Step: " << t << ", Observed: " << observedValue << ", Simulated: " << simulatedValue
                    << ", Normalized Observed: " << normalizedObservedValue << ", Normalized Simulated: " << normalizedSimulatedValue
                    << ", Error: " << error << ", Squared Error: " << squaredError << "\n"; // */
            }
        }

        //detailedLog.close();

        // Set the fitness value (e.g., NSSE)
        double fitness = sumSquaredErrors; // Increase the magnitude by not dividing by the count
        chromosome.setFitness(fitness);

        // Return the computed fitness
        return fitness;
    }

    std::vector<double> GeneticAlgorithm::extractRoughnessValues(EN_Project project) {
        std::vector<double> roughnessValues;
        int pipeCount = 0;

        // Get the number of pipes in the EPANET project
        if (EN_getCount(EN_PIPECOUNT, &pipeCount, project) != 0) {
            throw std::runtime_error("Failed to get the number of pipes from the EPANET project.");
        }

        //std::cout << "Pipe Count: " << pipeCount << std::endl;

        // Initialize the vector with a fixed size and default values
        roughnessValues.resize(pipeCount, 0.0);

        // Loop over each pipe and get its roughness value
        for (int i = 0; i < pipeCount; ++i) { // EPANET API is 1-based
            int linkType = 0;
            double roughness = 0.0;

            // Get the type of the link
            int linkTypeStatus = EN_getLinkType(i, &linkType, project);
            //std::cout << "EN_getLinkType Status for Link " << i << ": " << linkTypeStatus << std::endl;

            if (linkTypeStatus != 0) {
                throw std::runtime_error("Failed to get the link type from the EPANET project.");
            }

            // Only process pipes
            if (linkType == EN_PIPE) {
                // Get the roughness value of the pipe
                int roughnessStatus = EN_getLinkValue(i, EN_ROUGHNESS, &roughness, project);
                //std::cout << "EN_getLinkValue Status for Link " << i << ": " << roughnessStatus << std::endl;

                if (roughnessStatus != 0) {
                    throw std::runtime_error("Failed to get the roughness value from the EPANET project.");
                }

                //std::cout << "Link " << i << " Roughness: " << roughness << std::endl;
                roughnessValues[i] = roughness; // Assign directly to the appropriate index
            }
        }

        // Ensure the size of roughnessValues is exactly pipeCount
        if (roughnessValues.size() != static_cast<size_t>(pipeCount)) {
            throw std::runtime_error("The number of extracted roughness values does not match the pipe count.");
        }

        // Debug output to verify the number of roughness values
        //std::cout << "Number of roughness values extracted: " << roughnessValues.size() << std::endl;

        return roughnessValues;
    }

    // Function to back up the EPANET input file with the given roughness values
    void GeneticAlgorithm::backupEPANETInputFile(const std::string& baseFilePath, const std::vector<double>& roughnessValues, int generation, int populationIndex, int individualIndex) {
        std::ostringstream backupFilePath;
        backupFilePath << baseFilePath.substr(0, baseFilePath.find_last_of(".")) << "_gen_" << generation << "_pop_" << populationIndex << "_ind_" << individualIndex << ".inp";

        std::ifstream inputFile(baseFilePath);
        if (!inputFile) {
            throw std::runtime_error("Failed to open the base input file: " + baseFilePath);
        }

        std::ofstream outputFile(backupFilePath.str());
        if (!outputFile) {
            throw std::runtime_error("Failed to create the backup input file: " + backupFilePath.str());
        }

        std::string line;
        bool inPipesSection = false;
        size_t pipeIndex = 0;

        while (std::getline(inputFile, line)) {
            if (line.find("[PIPES]") != std::string::npos) {
                inPipesSection = true;
            }

            if (inPipesSection && !line.empty() && line[0] != ';') {
                std::istringstream lineStream(line);
                std::string id, node1, node2, length, diameter, roughness, minorLoss, status;

                if (!(lineStream >> id >> node1 >> node2 >> length >> diameter >> roughness >> minorLoss >> status)) {
                    outputFile << line << "\n";
                    continue;
                }

                if (pipeIndex < roughnessValues.size()) {
                    roughness = std::to_string(roughnessValues[pipeIndex]);
                    pipeIndex++;
                }

                outputFile << id << "\t" << node1 << "\t" << node2 << "\t"
                    << length << "\t" << diameter << "\t" << roughness << "\t"
                    << minorLoss << "\t" << status << "\n";
            }
            else {
                outputFile << line << "\n";
            }
        }
    }

    // Helper function to read roughness values from the EPANET project
    std::vector<double> GeneticAlgorithm::readRoughnessValuesFromProject(EN_Project project) {
        std::vector<double> roughnessValues;
        int numPipes;
        EN_getCount(EN_PIPECOUNT, &numPipes, project);

        for (int i = 1; i <= numPipes; ++i) {
            double roughness;
            EN_getLinkValue(i, EN_ROUGHNESS, &roughness, project);
            roughnessValues.push_back(roughness);
        }

        return roughnessValues;
    }

    void GeneticAlgorithm::adaptiveRates() {
        double averageFitness = 0.0;
        for (const Chromosome& chrom : population) {
            averageFitness += chrom.getFitness();
        }
        averageFitness /= population.size();

        double diversity = calculateDiversity(population);

        const double minMutationRate = 0.01;
        const double maxMutationRate = 0.3;

        // Increase mutation rate if diversity is low
        if (diversity < 0.5) {  // Adjust this threshold as needed
            mutationRate = std::min(maxMutationRate, mutationRate * 1.2);
        } else {
            mutationRate = std::max(minMutationRate, mutationRate * 0.9);
        }

        const double minCrossoverRate = 0.7;
        const double maxCrossoverRate = 0.95;

        if (crossoverRate < minCrossoverRate) {
            crossoverRate = minCrossoverRate;
        } else if (crossoverRate > maxCrossoverRate) {
            crossoverRate = maxCrossoverRate;
        }
    }

    bool GeneticAlgorithm::fitnessStable() {
        static const int stabilityThreshold = 10;  // Number of generations to consider
        static std::deque<double> recentBestFitnesses;

        double currentBestFitness = population[0].getFitness();
        recentBestFitnesses.push_back(currentBestFitness);

        if (recentBestFitnesses.size() > stabilityThreshold) {
            recentBestFitnesses.pop_front();
        }

        if (recentBestFitnesses.size() < stabilityThreshold) {
            return false;
        }

        double first = recentBestFitnesses.front();
        double last = recentBestFitnesses.back();
        return std::abs(last - first) < 0.1;  // Consider fitness stable if change is small
    }

    // Differential Evolution strategy to create a new candidate solution
    Chromosome GeneticAlgorithm::differentialEvolution(const Chromosome& target, const Chromosome& donor1, const Chromosome& donor2, const Chromosome& donor3) {
        std::vector<std::string> targetGenes = target.getGenes();  // Assume getGenes() returns std::vector<std::string>
        std::vector<std::string> donor1Genes = donor1.getGenes();
        std::vector<std::string> donor2Genes = donor2.getGenes();
        std::vector<std::string> donor3Genes = donor3.getGenes();

        std::vector<std::string> newGenes(targetGenes.size());

        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        std::uniform_int_distribution<size_t> dist(0, targetGenes.size() - 1);

        size_t R = dist(gen);
        double adaptiveDifferentialWeight = 0.8 * (1.0 - static_cast<double>(generationCount) / generations);
        double adaptiveCrossoverRate = crossoverRate * (1.0 - static_cast<double>(generationCount) / generations);

        for (size_t i = 0; i < targetGenes.size(); ++i) {
            if (dis(gen) < adaptiveCrossoverRate || i == R) {
                // Decode the binary strings to real numbers
                double donor1Real = donor1.decodeBinaryToReal(donor1Genes[i], 5.0, 140.0);
                double donor2Real = donor2.decodeBinaryToReal(donor2Genes[i], 5.0, 140.0);
                double donor3Real = donor3.decodeBinaryToReal(donor3Genes[i], 5.0, 140.0);

                // Perform the arithmetic operation with adaptive weight
                double newGeneReal = donor1Real + adaptiveDifferentialWeight * (donor2Real - donor3Real);

                // Ensure the new gene value is within bounds using custom clamp
                newGeneReal = clamp(newGeneReal, 5.0, 140.0);

                // Encode the result back into a binary string
                newGenes[i] = Chromosome::doubleToBinaryString(newGeneReal, 5.0, 140.0, 16);
            } else {
                newGenes[i] = targetGenes[i];
            }
        }

        return Chromosome(newGenes);  // Use the constructor to create a new Chromosome
    }

} //namespace GA


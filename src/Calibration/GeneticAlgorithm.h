// GeneticAlgorithm.h
#ifndef GENETIC_ALGORITHM_H
#define GENETIC_ALGORITHM_H
#include "Chromosome.h"
#include "epanet3.h"
#include <vector>
#include <fstream>
#include <filesystem>

namespace GA {

    class GeneticAlgorithm {
    public:
        GeneticAlgorithm(double mutationRate, double crossoverRate, size_t populationSize, size_t generations, EN_Project& project);
        int initializePopulation();
        void run();
        double computeFitnessAfterRun(EN_Project project, Chromosome& chromosome, int chromosomeIndex, int generationCount, const std::map<std::string, std::vector<double>>& simulatedValues);
        Chromosome getBestSolution() const;
        void evolve();
        void writeRoughnessAndFitnessValues(const std::vector<Chromosome>& population, int generation);
        std::vector<double> extractRoughnessValues(EN_Project project);
        void generateBestSolutionSimulatedFile(const std::string& baseFilename, int generationNumber);
        void simulatePopulation(size_t start, size_t end, int gen);
        void randomRestart();
        bool noImprovementForNGenerations(int n);
        bool fitnessStable();
        double minMutationRate = 0.001; // Set an appropriate minimum mutation rate
        double maxMutationRate = 0.5;   // Set an appropriate maximum mutation rate
    private:
        void backupEPANETInputFile(const std::string& baseFilePath, const std::vector<double>& roughnessValues, int generation, int populationIndex, int individualIndex);
        double mutationRate;
		double crossoverRate;
        size_t populationSize;
        size_t generations;
        EN_Project& project;
        std::vector<Chromosome> population;
        Chromosome selectParent();
        Chromosome bestOverall;
        double bestOverallFitness; 
		double calculateDiversity(const std::vector<Chromosome>& population) const;
        int generationCount = 0;
        std::vector<double> readRoughnessValuesFromInputFile(const std::string& filePath);
        double calculateTotalFitness() const;
        std::vector<double> computeSelectionProbabilities(double totalFitness) const;
        std::vector<double> computeCumulativeProbabilities(const std::vector<double>& selectionProbabilities) const;
        Chromosome rouletteWheelSelection() const;
        Chromosome binaryTournamentSelection();
        std::vector<Chromosome> selectFittestMuLambda(const std::vector<Chromosome>& parents, const std::vector<Chromosome>& offspring);
        void evaluateFitnessInParallel(std::vector<Chromosome>& population);
        void adaptiveRates();
        Chromosome differentialEvolution(const Chromosome& target, const Chromosome& donor1, const Chromosome& donor2, const Chromosome& donor3);
        std::vector<double> readRoughnessValuesFromProject(EN_Project project);
        
    };

} // namespace GA


#endif // GENETIC_ALGORITHM_H



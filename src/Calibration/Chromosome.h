/* EPANET 3 Genetic Algorithm Based Calibration Module
 *
 *
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H
#pragma once
#include <vector>
#include <string>
#include "epanet3.h"
#include <unordered_map>
#include <map>


class Chromosome {
public:
    Chromosome(size_t geneSize);
    
    Chromosome() = default; // Default constructor
    Chromosome(const std::vector<std::string>& genes) : genes(genes) {} // Constructor with genes

    void computeFitness(EN_Project project, int chromosomeIndex, int generationNumber);
    
    void mutate(double mutationRate, double generationFactor);
    Chromosome crossover(const Chromosome& other, double crossoverRate) const;

    Chromosome(const Chromosome& other);  // Copy constructor
    Chromosome(Chromosome&& other) noexcept;  // Move constructor
    Chromosome& operator=(const Chromosome& other);  // Copy assignment
    Chromosome& operator=(Chromosome&& other) noexcept;  // Move assignment
    static double decodeBinaryToReal(const std::string& binary, double lowerBound, double upperBound);
    void setGene(int index, const std::string& geneValue);
    std::vector<double> applyChromosomeToProject(EN_Project project);

    std::string getGene(int index) const;
    double getFitness() const; 
    void copySolution(const Chromosome& other);
    static bool isTank(const std::string& id);
    const std::vector<std::string>& getGenes() const;
    int getIdIndex(const std::string& id, EN_Project project);
    
    static std::string doubleToBinaryString(double value, double lowerBound, double upperBound, size_t bits);
    std::vector<double> getRoughnessValues() const;
    static std::string generateBinaryString(int length);

    double distance(const Chromosome& other) const;
	void setFitness(double fitness);
    std::vector<std::string> genes;
    static std::map<std::string, std::vector<double>> readSimulatedValues(const char* filename);
    static std::map<std::string, std::vector<double>> readObservedValues(const std::string& filename);
    void reset();
private:
    static double randomDouble(double min, double max);
    
    double fitness;
};

void updateInputFileWithRoughness(const std::string& filePath, const std::vector<double>& roughnessValues);

// Custom clamp function for compatibility with C++ versions before C++17
template<typename T>
T clamp(T value, T lower, T upper) {
    return std::min(std::max(value, lower), upper);
}// */

#endif // CHROMOSOME_H

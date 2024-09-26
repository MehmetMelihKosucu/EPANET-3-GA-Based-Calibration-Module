#include "Chromosome.h"
#include <cmath>
#include <algorithm>
#include <random>
#include "epanet3.h"
#include <iostream>
#include <bitset>
#include <stdexcept>
#include <limits>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <numeric>
#include <filesystem>
#include <functional>

// Default constructor
//Chromosome::Chromosome() : fitness(0.0) {}

Chromosome::Chromosome(size_t geneSize) : genes(geneSize), fitness(0.0) {
    // Initialize genes with random binary strings
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);
    for (auto& gene : genes) {
        gene.resize(geneSize);
        for (auto& bit : gene) {
            bit = dis(gen) ? '1' : '0';
        }
    }
}

Chromosome::Chromosome(const Chromosome& other) : genes(other.genes), fitness(other.fitness) {}

Chromosome::Chromosome(Chromosome&& other) noexcept : genes(std::move(other.genes)), fitness(other.fitness) {
    other.fitness = 0.0;
}

Chromosome& Chromosome::operator=(const Chromosome& other) {
    if (this != &other) {
        genes = other.genes;
        fitness = other.fitness;
    }
    return *this;
}

Chromosome& Chromosome::operator=(Chromosome&& other) noexcept {
    if (this != &other) {
        genes = std::move(other.genes);
        fitness = other.fitness;
        other.fitness = 0.0;
    }
    return *this;
}

// Helper function to generate a random double value within a given range
double Chromosome::randomDouble(double min, double max) {
    static std::random_device rd;
    static std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(min, max);
    return dist(mt);
}

std::string Chromosome::doubleToBinaryString(double value, double lowerBound, double upperBound, size_t bits) {
    double normalizedValue = (value - lowerBound) / (upperBound - lowerBound);
    unsigned long integerRepresentation = static_cast<unsigned long>(normalizedValue * (std::pow(2, bits) - 1));
    std::bitset<64> bitsetRepresentation(integerRepresentation);
    std::string binaryString = bitsetRepresentation.to_string();
    return binaryString.substr(64 - bits);
}

// Helper function to generate a random binary string of a given length
std::string Chromosome::generateBinaryString(int length) {
    std::string binaryString;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::bernoulli_distribution distribution(0.5); // 50% chance for 0 or 1

    for (int i = 0; i < length; ++i) {
        binaryString += distribution(generator) ? '1' : '0';
    }
    return binaryString;
}
const std::vector<std::string>& Chromosome::getGenes() const {
    return genes;
}

/*// Custom clamp function for C++11 and C++14 compatibility
template<typename T>
T clamp(T value, T lower, T upper) {
    return std::min(std::max(value, lower), upper);
} */

// Function to set the value of a gene at a given index
void Chromosome::setGene(int index, const std::string& gene) {
    if (index >= 0 && static_cast<size_t>(index) < genes.size()) {
        // Validate that the value is a binary string
        if (!std::all_of(gene.begin(), gene.end(), [](char c) { return c == '0' || c == '1'; })) {
            throw std::invalid_argument("Gene value must be a binary string containing only '0' and '1'. Value provided: " + gene);
        }
        genes[index] = gene;
    }
    else {
        throw std::out_of_range("Gene index out of range");
    }
}

// Function to copy the solution from one chromosome to another
void Chromosome::copySolution(const Chromosome& other) {
    this->fitness = other.getFitness();
    this->genes = other.genes;
}

// Function to get the gene value at a given index
std::string Chromosome::getGene(int index) const {
    if (index < 0 || index >= static_cast<int>(genes.size())) {
        throw std::out_of_range("Gene index out of range");
    }
    return genes[index];
}

// Function to decode a binary string to a real value within a given range
double Chromosome::decodeBinaryToReal(const std::string &binaryString, double minValue, double maxValue) {
    // Check if binaryString is a valid binary number
    if (!std::all_of(binaryString.begin(), binaryString.end(), [](char c) { return c == '0' || c == '1'; })) {
        throw std::invalid_argument("Invalid binary string");
    }

    // Ensure binaryString is not empty and does not exceed 64 bits
    if (binaryString.empty() || binaryString.size() > 64) {
        throw std::invalid_argument("Binary string size must be between 1 and 64");
    }

    // Convert binary string to decimal value
    unsigned long long decimalValue = std::bitset<64>(binaryString).to_ullong();

    // Normalize the decimal value to a value between 0 and 1
    double normalized = static_cast<double>(decimalValue) / static_cast<double>((1ULL << binaryString.size()) - 1);

    // Scale and return the value within the specified range
    return minValue + normalized * (maxValue - minValue);
}
std::vector<double> Chromosome::applyChromosomeToProject(EN_Project project) {
    if (!project) {
        throw std::invalid_argument("Invalid EPANET project reference.");
    }

    int numPipes;
    int errorCode = EN_getCount(EN_PIPECOUNT, &numPipes, project);
    if (errorCode != 0) {
        std::cerr << "Failed to retrieve number of pipes from EPANET. Error code: " << errorCode << std::endl;
        throw std::runtime_error("Failed to retrieve number of pipes from EPANET.");
    }

    std::vector<double> roughnessValues(numPipes);

    for (int i = 0; i < numPipes; ++i) {
        try {
            roughnessValues[i] = decodeBinaryToReal(getGene(i), 5, 140);
        } catch (const std::exception& e) {
            std::cerr << "Error decoding binary to real for pipe " << i + 1 << ": " << e.what() << std::endl;
            roughnessValues[i] = -1;  // Use an invalid value to indicate failure
        }

        if (roughnessValues[i] != -1) {
            errorCode = EN_setLinkValue(i + 1, EN_ROUGHNESS, roughnessValues[i], project);
            if (errorCode != 0) {
                std::cerr << "Failed to set roughness for pipe " << i + 1 << ". Error code: " << errorCode << std::endl;
                throw std::runtime_error("Critical error setting pipe values, halting execution.");
            }
        }
    }

    // Update the input file with new roughness values
    std::string filePath = "/home/kosucu/EPANET-GA-OMP/Networks/C-Town/c-town_true_network-Uncalibrated-EPA3.inp";
    updateInputFileWithRoughness(filePath, roughnessValues);

    return roughnessValues;
}

std::vector<double> Chromosome::getRoughnessValues() const {
    size_t numGenes = genes.size();
    std::vector<double> roughnessValues(numGenes);

    for (size_t i = 0; i < numGenes; ++i) {
        try {
            roughnessValues[i] = decodeBinaryToReal(genes[i], 5, 140);
        }
        catch (const std::exception& e) {
            roughnessValues[i] = -1; // Indicate failure with an invalid value
        }
    }
    return roughnessValues;
}

void updateInputFileWithRoughness(const std::string& filePath, const std::vector<double>& roughnessValues) {
    std::ifstream inputFile(filePath);
    if (!inputFile) {
        throw std::runtime_error("Failed to open the input file: " + filePath);
    }

    std::stringstream buffer;
    buffer << inputFile.rdbuf();  // Read the entire file into a buffer
    inputFile.close();

    std::string content = buffer.str(), line;
    std::istringstream contentStream(content);
    std::ostringstream newContent;
    bool inPipesSection = false;
    size_t pipeIndex = 0;  // Keep track of pipe index

    while (getline(contentStream, line)) {
        // Check for the start of the pipes section
        if (line.find("[PIPES]") != std::string::npos) {
            inPipesSection = true;
            newContent << line << "\n";
            continue;
        }

        // Check for the end of the pipes section or an unexpected section marker
        if (inPipesSection && (line.empty() || line[0] == ';' || line.find("[") != std::string::npos)) {
            if (line.find("[") != std::string::npos && line.find("[PIPES]") == std::string::npos) {
                inPipesSection = false; // Handle case where a new section starts
            }
            newContent << line << "\n";
            continue;
        }

        if (inPipesSection && !line.empty() && line[0] != ';') {
            std::istringstream lineStream(line);
            std::string id, node1, node2, length, diameter, roughness, minorLoss, status;

            // Parse the line assuming tab or space separated values
            if (!(lineStream >> id >> node1 >> node2 >> length >> diameter >> roughness >> minorLoss >> status)) {
                newContent << line << "\n";  // Maintain original line if parsing fails
                continue;
            }

            if (pipeIndex < roughnessValues.size()) {
                roughness = std::to_string(roughnessValues[pipeIndex]);  // Update roughness
                pipeIndex++;
            }

            // Ensure formatting with fixed spacing or tabs for alignment
            newContent << id << "\t" << node1 << "\t" << node2 << "\t"
                << length << "\t" << diameter << "\t" << roughness << "\t"
                << minorLoss << "\t" << status << "\n";
        }
        else {
            newContent << line << "\n";
        }
    }

    // Write the new content back to the same file
    std::ofstream outputFile(filePath);
    if (!outputFile) {
        throw std::runtime_error("Failed to open the input file for writing: " + filePath);
    }
    outputFile << newContent.str();
}

void Chromosome::computeFitness(EN_Project project, int chromosomeIndex, int generationNumber) {
    // Read observed and simulated values from files
    //auto observedValues = readObservedValues("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Observed.dat"); // For Windows
    auto observedValues = readObservedValues("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Observed.dat"); // For Linux
    //auto simulatedValues = readSimulatedValues("E:\\EPANET-GA-Parallel\\Networks\\C-Town\\Simulated.dat"); // For Windows
    auto simulatedValues = readSimulatedValues("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/Simulated.dat"); // For Linux

    double sumSquaredErrors = 0.0;
    int count = 0;  // To keep track of the number of valid comparisons
    double tankWeight = 10.0;  // Weighting factor for tank pressure values
    double scalingFactor = 100.0; // Scaling factor to ensure errors are in a reasonable range

    std::vector<double> maxObserved(observedValues.begin()->second.size(), std::numeric_limits<double>::min());
    std::vector<double> minObserved(observedValues.begin()->second.size(), std::numeric_limits<double>::max());

    for (const auto& obsPair : observedValues) {
        const std::vector<double>& obsData = obsPair.second;
        for (size_t i = 0; i < obsData.size(); ++i) {
            maxObserved[i] = std::max(maxObserved[i], obsData[i]);
            minObserved[i] = std::min(minObserved[i], obsData[i]);
        }
    }

    //std::ofstream detailedLog("/home/kosucu/EPANET-GA-OMP/Networks/C-Town/detailed_fitness_log_" + std::to_string(generationNumber) + "_" + std::to_string(chromosomeIndex) + ".dat");

    for (const auto& obsPair : observedValues) {
        const std::string& id = obsPair.first;
        const std::vector<double>& obsData = obsPair.second;

        if (simulatedValues.count(id) == 0) continue;

        const std::vector<double>& simData = simulatedValues[id];

        for (size_t i = 0; i < obsData.size() && i < simData.size(); ++i) {
            double observedValue = obsData[i];
            double simulatedValue = simData[i];
            double error = std::abs(simulatedValue - observedValue);

            double range = maxObserved[i] - minObserved[i];
            double normalizedError = (range != 0) ? (error / range) * scalingFactor : error;
            double weightedError = isTank(id) ? normalizedError * tankWeight : normalizedError;

            sumSquaredErrors += weightedError * weightedError;
            count++;

            /*detailedLog << "Generation: " << generationNumber << ", ID: " << id << ", Observed: " << observedValue << ", Simulated: " << simulatedValue
                << ", Error: " << error << ", Normalized Error: " << normalizedError
                << ", Weighted Error: " << weightedError << "\n"; // */
        }
    }
    //detailedLog.close();

    double wsse = count > 0 ? sumSquaredErrors / count : std::numeric_limits<double>::infinity();

    fitness = wsse;
}

std::map<std::string, std::vector<double>> Chromosome::readObservedValues(const std::string& filename) {
    std::map<std::string, std::vector<double>> observedValues;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return observedValues;
    }

    std::string line;
    // Assuming the first line of the file may contain headers
    std::getline(file, line);  // Read and discard the header line

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string time;
        double value;

        // Read and discard the time column
        iss >> time;

        // Tank identifiers
        std::vector<std::string> ids = { "T1", "T2", "T3", "T4", "T5", "T6", "T7", "S1-P310", "S2-PU4", "S3-PU7", "S4-PU8", "S5-PU10" };

        for (const auto& id : ids) {
            if (!(iss >> value)) {
                std::cerr << "Error reading value for " << id << " in line: " << line << std::endl;
                continue;
            }
            observedValues[id].push_back(value);
        }
    }

    return observedValues;
}

void Chromosome::setFitness(double fitness) {
    this->fitness = fitness;
}

double Chromosome::getFitness() const {
    return fitness;
}

std::map<std::string, std::vector<double>> Chromosome::readSimulatedValues(const char* filename) {
    std::map<std::string, std::vector<double>> simulatedValues;
    std::ifstream file(filename);

    if (!file) {
        std::cerr << "Could not open file: " << filename << std::endl;
        return simulatedValues;
    }

    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read the header line from file: " << filename << std::endl;
        return simulatedValues;
    }

    std::vector<std::string> ids = { "T1", "T2", "T3", "T4", "T5", "T6", "T7", "S1-P310", "S2-PU4", "S3-PU7", "S4-PU8", "S5-PU10" };

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string time;
        double value;
        std::vector<double> values;

        if (!(iss >> time)) {
            std::cerr << "Error reading time from line: " << line << std::endl;
            continue;
        }

        bool errorReadingLine = false;
        for (const auto& id : ids) {
            if (!(iss >> value)) {
                std::cerr << "Error reading value for " << id << " in line: " << line << std::endl;
                errorReadingLine = true;
                break;
            }
            values.push_back(value);
        }

        if (!errorReadingLine) {
            for (size_t i = 0; i < ids.size(); ++i) {
                simulatedValues[ids[i]].push_back(values[i]);
            }
        }
    }

    return simulatedValues;
}

bool Chromosome::isTank(const std::string& id) {
    // Assuming tank IDs start with 'T'
    return !id.empty() && id[0] == 'T';
}

int Chromosome::getIdIndex(const std::string& id, EN_Project project) {

    // Assuming idIndexMap should be static or should be a member variable that's initialized elsewhere
    static std::map<std::string, int> idIndexMap;

    auto it = idIndexMap.find(id);
    if (it != idIndexMap.end()) {
        return it->second;  // Found the index for the id
    }

    int index = -1;
    int errorCode;

    // Create a modifiable copy of id as a C-style string
    std::vector<char> idCStr(id.begin(), id.end());
    idCStr.push_back('\0'); // Null-terminate the C-style string

    // Check the prefix of the id for determining its type
    if (id.substr(0, 1) == "T") {  // If id starts with 'T'
        errorCode = EN_getNodeIndex(idCStr.data(), &index, project);
    }
    else if (id.substr(0, 1) == "P" || id.substr(0, 2) == "PU") {  // If id starts with 'P' or 'PU'
        errorCode = EN_getLinkIndex(idCStr.data(), &index, project);
    }
    else {
        throw std::runtime_error("Unknown ID type: " + id);
    }

    if (errorCode != 0) {
        throw std::runtime_error("Failed to get index for ID: " + id + " with error code: " + std::to_string(errorCode));
    }

    idIndexMap[id] = index;
    return index;
}


void Chromosome::mutate(double mutationRate, double generationFactor) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    double adjustedMutationRate = mutationRate * generationFactor; // Adjust the mutation rate

    for (auto& gene : genes) {
        if (gene.size() < 1 || gene.size() > 64) {
            throw std::invalid_argument("Gene length must be between 1 and 64 during mutation");
        }
        for (auto& bit : gene) {
            if (dis(gen) < adjustedMutationRate) {
                bit = (bit == '0') ? '1' : '0'; // Flip the bit
            }
        }
    }
}

void Chromosome::reset() {
    // Clear the gene vector
    genes.clear(); // This clears all elements in the vector, effectively resetting it

    // Reset fitness to a default value
    fitness = 0.0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 1);

    for (auto& gene : genes) {
        gene.resize(64);  // Ensure that genes have the desired length
        for (auto& bit : gene) {
            bit = dis(gen) ? '1' : '0';
        }
    }
}

Chromosome Chromosome::crossover(const Chromosome& other, double crossoverRate) const {
    if (this->genes.size() != other.genes.size()) {
        throw std::invalid_argument("Parent chromosomes must have the same number of genes.");
    }

    Chromosome offspring(genes.size());
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dis(0, 1);

    double eta = 2.0; // SBX parameter

    // Perform SBX for each gene
    for (size_t i = 0; i < genes.size(); ++i) {
        if (dis(gen) < crossoverRate) {
            // Decode the binary strings to real values
            double parent1Value = decodeBinaryToReal(this->genes[i], 5, 140);
            double parent2Value = decodeBinaryToReal(other.genes[i], 5, 140);

            // Apply SBX
            double u = dis(gen);
            double beta;
            if (u <= 0.5) {
                beta = pow(2.0 * u, 1.0 / (eta + 1.0));
            }
            else {
                beta = pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (eta + 1.0));
            }
            double child1Value = 0.5 * ((1 + beta) * parent1Value + (1 - beta) * parent2Value);
            double child2Value = 0.5 * ((1 - beta) * parent1Value + (1 + beta) * parent2Value);

            // Introduce a small random change to encourage diversity
            std::normal_distribution<double> mutation_dist(0.0, 1.0);
            double mutation_factor = mutation_dist(gen) * 0.01; // 1% random change
            child1Value += mutation_factor;
            child2Value += mutation_factor;

            // Ensure the new gene value is within bounds
            child1Value = clamp(child1Value, 5.0, 140.0);
            child2Value = clamp(child2Value, 5.0, 140.0);

            // Encode the real values back to binary strings
            offspring.genes[i] = doubleToBinaryString(child1Value, 5, 140, genes[i].size());
        }
        else {
            // No crossover, copy the gene from one of the parents
            offspring.genes[i] = this->genes[i];
        }
    }

    return offspring;
}

double Chromosome::distance(const Chromosome& other) const {
    double dist = 0.0;

    // Check for gene size mismatch
    if (genes.size() != other.genes.size()) {
        std::cerr << "Error: Chromosome gene size mismatch in distance calculation.\n"
                  << "This gene size: " << genes.size() << "\n"
                  << "Other gene size: " << other.genes.size() << std::endl;
        throw std::runtime_error("Chromosome gene size mismatch in distance calculation.");
    }

    for (size_t i = 0; i < genes.size(); ++i) {
        if (genes[i].size() != other.genes[i].size()) {
            std::cerr << "Error: Gene size mismatch between Chromosomes at index " << i << ".\n"
                      << "This gene size: " << genes[i].size() << "\n"
                      << "Other gene size: " << other.genes[i].size() << std::endl;
            throw std::runtime_error("Gene size mismatch between Chromosomes.");
        }
    }

    // If no errors, calculate distance
    dist += std::inner_product(genes[0].begin(), genes[0].end(), other.genes[0].begin(), 0.0,
                               std::plus<double>(),
                               [](char a, char b) { return static_cast<double>(a != b); });

    return dist;
}






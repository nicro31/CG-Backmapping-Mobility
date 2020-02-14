#ifndef MEANFIELD_H_INCLUDED
#define MEANFIELD_H_INCLUDED
#include "Parse.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <fstream>
#include <time.h>
#include <omp.h>

#define outputLength 24
#define outputPrecision 8
#define dichoPrecision 1e-15

struct ConvergenceParameters
{
    doubleType m_accuracyMob;
    doubleType m_accuracyOcc;
    doubleType m_accuracyCon;
    int m_nbIterUpdate;
    int m_nbIterMax;  
};

doubleType adjustConcentration(double &pTemp, Graph &g); // Set the chemical potential to get a desired carrier concentration
doubleType computeMobility(Graph &g, bool usePrevDistrib, std::string folder, ConvergenceParameters &convParam); // Compute graph mobility, for a graph with parameters already set (field, temperature, concentration)
Vect3 findPercolationThreshold(Graph &g, std::string folder); // Find the percolation threshold (in x, y, z) for a graph
double findPercolationThresholdN(Graph &g, std::string folder); // Find the percolation threshold in nb of molecules per cluster
void prepareSystem(std::string inputGroFile, std::string inputConfigFile, std::string outputPdbFile, int generateGaussian); // Prepare pdb files and .com files
void prepareSystemGaussian(std::string inputPdbFile, std::string inputConfigFile); // Prepare .com files
int plotPercolation(Graph &g, std::string percolationFile, bool dimCluster); // Compute the percolation behavior of the graph
void computeMobilityDependance(Graph &g, std::vector<Vect3> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder, RateType rate, ConvergenceParameters &convParam); // Compute the mobility for different values of the parameters field, tempearture and concentration
Graph initGraph(std::string inputPdbFile, std::string inputConfigFile, std::string inputGraphFile); // Prepare a graph from pbd file coordFile, with parameters in configFile, and transfer integrals and energies in the graph file
void readMobilityParam(std::string configFile, std::vector<Vect3> &pField, std::vector<double> &pTemperature, std::vector<double> &pConcentration, ConvergenceParameters &pConvergence); // Read parameters related to mobility calculation from config file
void readConfigParam(std::string configFile, Graph &g); // Read parameters from configuration file
int buildGraph(std::string inputPdbFile, std::string inputConfigFile, std::string outputGraphFile, int useZindo); // Read results from ZINDO calculation, calculate energies and transfer integrals, and save all of it. Or calculate transfer with Miller-Abrahams if useZindo = 0
doubleType checkTransportConvergence(Graph &g); // Check solution of the mobility calculation

/// Export local current distribution ///
void exportLocalCurrent(Graph &g,  std::vector<Vect3> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder, RateType rate);



#endif // MEANFIELD_H_INCLUDED

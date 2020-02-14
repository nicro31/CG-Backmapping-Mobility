#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED
#include "System.h"
#include "Johnson.h"
#include "QC_Tools/src/calculateTI.hpp"

#ifdef _WIN32 || _WIN64
#include <direct.h>
#else
#include <sys/stat.h>
#endif // linux

enum HopType { intra , inter, dopantIntra, dopantInter, level };
enum RateType { maPheno, maTransfer };

struct conjSegment
{

    /// Network parameters ///
    Vect3 m_pos;
    Vect3 m_orientN;
    Vect3 m_orientP;
    std::vector<RigidFragment*> m_fragment;
    std::vector<conjSegment*> m_neighbor; // Site neighbors
    std::vector<int> m_neighborIdx; // Site neighbor index
    std::vector<std::vector<doubleType>> m_transferIntegral; // Vector of transfer integral to neighbors, for the different levels
    std::vector<Vect3> m_neighborDisplacement; // Displacement vectors to the neighbors, taking into account pbc
    std::vector<Vect3> m_neighborPbc; // What pbc (x,y and/or z, with + or -) links the neighbors ?
    std::vector<std::vector<Bound> > m_boundFragment; // Rigid fragment constituting the hopping bound
    std::vector<std::vector<int> > m_ionBridge; // Counterion forming the bound if it exists, else -1
    std::vector<std::vector<HopType> > m_hopType; // Type of hopping for each fragment bound
    std::vector<doubleType> m_hopTo; // Hopping rate to its neighbor
    std::vector<doubleType> m_hopFrom; // Hopping rate from its neighbor
    int m_clusterId; // To what cluster this conjugated segment belongs in the graph
    std::vector<bool> m_periodic; // Is it bounded by pbc
    std::vector<doubleType> m_displacementField; // Contains r.F to speed up mean field calculation

    /// Mean field parameters ///
    double m_energy; // Site energy for transport (considering for instance that a hole is always relaxing to the HOMO before hopping to another molecule)
    std::vector<double> m_energyLevels; // List of energy levels on the site (HOMO, HOMO-1, HOMO-2... or LUMO, LUMO+1, LUMO+2...)
    doubleType m_occup; // Site occupation probability
    doubleType m_occupT; // Temp site occupation probability
    doubleType m_cP; // Local chemical potential

    /// Molecule parameters ///
    int m_level; // This site corresponds to HOMO - m_level MO on the molecule or LUMO + m_level

};

struct Graph
{
    Graph() : m_transferMat(0), m_vectLength(1.5), m_vectRadius(0.2), m_vectResol(12)
    {

    }

    ~Graph()
    {
        if(m_transferMat != 0)
        {
            int n(m_site.size()/m_nbrLevels);
            for (int i = 0; i < n; i++)
            {
                delete[] m_transferMat[i];
            }
            delete[] m_transferMat;
            m_transferMat = 0;
        }

    }

    /// Input files for the graph ///
    std::string m_fileNameTop;
    std::string m_fileNameCoord;


    /// Graph structure (nodes and atoms) ///
    System m_system;
    std::vector<conjSegment> m_site;
    int m_edgeNbAll;
    int m_edgeNbIntra;
    int m_edgeNbInter;
    int m_edgeNbDopantIntra;
    int m_edgeNbDopantInter;
    int m_edgeNbAll_pbc; // Number of edges without considering edges happening because of pbc
    int m_edgeNbIntra_pbc;
    int m_edgeNbInter_pbc;
    int m_edgeNbDopantIntra_pbc;
    int m_edgeNbDopantInter_pbc;

    /// Cluster parameters ///
    int m_clusterNb;
    std::vector<Vect3> m_clusterDim;
    std::vector<int> m_clusterSiteNb;


    /// Edges threshold parameters ///
    double m_dThresholdInter;
    double m_dThresholdIntra;
    double m_dThresholdDopantInter;
    double m_dThresholdDopantIntra;
    double m_sphereInfluenceRad;

    /// External parameters for mobility computation ///
    Vect3 m_field;
    double m_temperature;
    double m_concentration; // Number of sites occupied
    std::vector<int> m_orderMeanField; // Ordering of the sites for mobility calculation

    /// Miller Abrahams parameters ///
    double m_w0_intra;
    double m_w0_inter;
    double m_w0_level;
    double m_w0_dopantIntra;
    double m_w0_dopantInter;
    double m_loc_inter;
    double m_loc_dopantInter;
    double m_loc_intra;
    double m_loc_dopantIntra;
    double m_H_intra;
    double m_H_level;

    /// Transfer Integrals Matrix ///
    std::vector<double> **m_transferMat;
    std::vector<doubleType> m_maxTransfer;
    std::vector<doubleType> m_minTransfer;
    int m_nbrLevels; // Number of considered MO on each molecule
    int m_nbrLevels2Read; // Number of levels to read from Gorilla
    std::string m_levelType; // Level type to use. HOMO for hole transport, LUMO for electron transport
    double m_dosBroadening; // DOS broadening, in eV. If not 0, all levels energies are picked from a gaussian DOS distribution with this broadening. If -1, energy levels come from ZINDO calculation.
    std::vector<double> m_levelEnergies; // Level energies from ZINDO
    int m_nbProcSharedGaussian; // Number of shared processors used for ZINDO calculations (will be written in the .com files)
    
    /// Parameters to calculate percolation curve ///
    double m_Hmin;
    double m_Hmax;
    double m_dH;


    /// Graphical parameters for VMD ///
    double m_vectLength;
    double m_vectRadius;
    double m_vectResol;


    /// Method ///
    void initSystem();
    void segmentSystem();
    void findNeighbors();
    void writeGraph(std::string fileNameGraphPdb, std::string fileNameGraphPsf, bool periodic);
    void cluster(double threshold, int clusteringMode); // Cluster the graph with a given threshold and mode. mode 0 is clustering based on linking only. mode 1 is based on hopping values
    double percolation();
    std::vector<std::vector<Edge>> createJohnsonGraph(int clusterId, Direction d);
    void computeRate(RateType rate);
    void writeCurrent(std::string fileNameOccupancyPdb, std::string fileNameCurrentPdb, std::string fileNameCurrentPsf, bool periodic);
    void gdm(double sigma, double mu, unsigned int seed, bool coarsed);
    void writeBoundGeometry(std::string fileNameGeometry);
    void writeBoundGeometryRF(std::string fileNameGeometry);
    void writeBoundGeometryRFDist(std::string fileNameGeometry);
    void writeFragmentVector(std::string fileNameFragmentVector);
    void readFragmentVector(std::string fileNameFragmentVector);
    void writeSiteVector(std::string fileNameSiteVector);
    void importTransferIntegral(std::vector<std::string> fileNameTransfer); // Import transfer integrals from gorilla build
    void writeTransferIntegral(std::string fileNameTransferPdb, std::string fileNameTransferPsf, bool periodic, int l1, int l2); // Write ln(|H|) in file to visualize with VMD ; int level : level to consider (eg transfer from HOMO to HOMO - level)
    void importEnergies(std::string fileNameEnergies); // Assign site energies based on a gorilla calculation
    void writeHistogramTransferDist(std::string fileNameHistogram, int nbBinTranfer, int nbBinDist, int level); // int level : level to consider (transfer from HOMO to HOMO - level)
    void writeRawTranfer(std::string fileNameRawTransfer);
    void writeRawRate(std::string fileNameRawRate);
    void writeRawEnergies(std::string fileNameRawEnergies);
    int writeGeometry(std::string geometryFile); // Write geometry of the links in the system
    void addMolLevels2Site(); // Add transport site corresponding to the different levels for each molecule (HOMO, HOMO-1, HOMO-2...)
    void orderSites(double latticeParam); // Arrange the site in a grid of lattice parameter latticeParam, for fastest convergence (?)
    void parseTransfer(); // Add missing transfer between sites bounded by pbc, in the case that we introduced artificially big pbc in Gorilla
    void writeSystemClustered(std::string fileNameCoord, std::string folder); // Write the system of pedot with their cluster iD (no tos)
    void writeboundDistribution(std::string boundDistrib); // write bound distribution (nb of bound linking each pair of molecule)
    void writeTransferDistanceAngleNN(std::string fileName, int l1, int l2, RateType rate); // Write geometrical info on pair of molecules to feed neural network
    void readCarrierDistribution(std::string fileNameDistrib); // read carrier distribution (i.e. occupation proba p_i) from a previous mobility calculation
    void generateGaussianCom(); // Generate Gaussian com files for the conjugated segment in the graph
    void getEnergyAndTransfer(); // Calculate transfer integrals and get energies from ZINDO calculation
    void getMillerAbrahams(); // Calculate transfer integrals with Miller-Abrahams
    int saveGraph(std::string outputGraphFile); // Save graph information (energies and transfer integrals)
    int loadGraph(std::string inputGraphFile); // Load graph information
    int writeEnergy(std::string energyFile); // Write DOS to text file
    int writeTransfer(std::string transferFile); // Write transfer integral distribution to text file
    int writeRate(std::string rateFile); // Write hopping rate distribution to text file
    int completeGraph(int seed); // Extend the graph with the different levels, attribute transfer integrals regarding pbc, eventually calculate a new DOS based on a given (e.g. gaussian distribution), calculate hopping rates
    int writeGraphVmd(std::string fileNameGraphVmd, int mode);
    void averageTransfer(); //Be sure transfer are symmetrics
};


#endif // GRAPH_H_INCLUDED

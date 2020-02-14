#include "MeanField.h"



using namespace std;


int main(int argc, char *argv[])
{
    unsigned int seed(time(NULL));
    
    ///// Available modes and corresponding list of arguments /////
    vector<string> modes{ "prep" , "prep_gaussian", "build" , "export_energy", "export_transfer", "export_rate", "export_geometry", "perco", "mobility", "export_current", "export_graph" };
    vector<vector<string>> args{
      { "input gromacs file (.gro)", "input config file (.txt)" , "output pdb file (.pdb)", "generate Gaussian files (1 for yes, 0 for no)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "output graph file (.grp)", "use ZINDO (1 for yes, 0 for no)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output energy file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output transfer file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output rate file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output geometry file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output percolation file (.txt)" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output mobility folder" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "input mobility folder" } ,
      { "input pdb file (.pdb)" , "input config file (.txt)", "input graph file (.grp)", "output graph file (.tcl)", "mode (0 for transfer, 1 for rates, 2 for current)", "input mobility folder (put dummy if not mode 2)" }
      
    };
    
    ///// Parsing arguments /////
    if(argc < 2) 
    { 
      cerr << "You should provide a mode for LOE-CTP-FRAG. Choices are:" << endl;
      for(int i(0) ; i < modes.size() ; i++)
      {
	cout << modes[i];
	if(i==modes.size()-1) {cout << endl;}
	else {cout << " / ";}
      }
      return -1;
      
    }
    
    ///// Parsing modes /////
    string mode;
    if(argc >= 2)
    {
      mode = argv[1];
      bool modeFound(false);
      for(int i(0) ; i < modes.size() ; i++)
      {
	if( mode == modes[i] )
	{
	  modeFound = true;
	  int nbArgs(argc - 2);
	  
	  if(nbArgs < args[i].size() || nbArgs > args[i].size() + 1)
	  {
	    cerr << "Wrong number of arguments for mode " << mode << ". You should provide the following arguments:" <<endl;
	    for(int j(0) ; j < args[i].size() ; j++)
	    {
	      cout << args[i][j] << " / ";
	    }
	    cout << "[optionnal] seed (unsigned integer)" << endl;
	    return -1;
	  }
	  
	  // Optionnal seed argument
	  if(nbArgs == args[i].size() + 1)
	  {
	    seed = atoi(argv[argc-1]);
	  }
	}
	
      }
      
      if(!modeFound)
      {
	cerr << "Provided mode for LOE-CTP-FRAG is incorrect. Choices are:" << endl;
	for(int i(0) ; i < modes.size() ; i++)
	{
	  cout << modes[i];
	  if(i==modes.size()-1) {cout << endl;}
	  else {cout << " / ";}
	}
	return -1;
      }
      
    }
    
    /////////////////////////////////////////////////////////
    ///// Exec program as a function of the chosen mode /////
    /////////////////////////////////////////////////////////
    
    
    
    ///// MODE = prep /////
    if( mode == "prep" )
    {
      std::string inputGroFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string outputPdbFile = argv[4];
      int generateGaussian = atoi(argv[5]);
      
      cout << "Preparing the system..." << endl;
      prepareSystem(inputGroFile, inputConfigFile, outputPdbFile, generateGaussian);
      cout << "The system has been prepared !" << endl << "Leaving LOE-CTP-FRAG" << endl;
      return 0;
        
    }
    
    ///// MODE = prep /////
    if( mode == "prep_gaussian" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      
      cout << "Preparing the system for Gaussian calculations..." << endl;
      prepareSystemGaussian(inputPdbFile, inputConfigFile);
      cout << "The system has been prepared for Gaussian calculations!" << endl << "Leaving LOE-CTP-FRAG" << endl;
      return 0;
        
    }
    
    
    ///// MODE = build /////
    if( mode == "build" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string outputGraphFile = argv[4];
      int useZindo = atoi(argv[5]);
      
      cout << "Building the graph..." << endl;
	  // Allow external program to check status via output file
	  std::ofstream buildStatus;
	  buildStatus.open("buildStatus.txt");
	  buildStatus << "running";
	  buildStatus.close();

      int building = buildGraph(inputPdbFile, inputConfigFile, outputGraphFile, useZindo);

      if(building == 0)
      {
	cout << "The graph has been built!" << endl << "Leaving LOE-CTP-FRAG" << endl;
	  buildStatus.open("buildStatus.txt");
	  buildStatus << "done";
	  buildStatus.close();
	return 0;
      }
      if(building == -1)
      {
	cerr << "Error during graph building!" << endl << "Leaving LOE-CTP-FRAG" << endl;
	  buildStatus.open("buildStatus.txt");
	  buildStatus << "crashed";
	  buildStatus.close();
	return -1;
      }
      if(building == -2)
      {
	cerr << "Error during graph saving!" << endl << "Leaving LOE-CTP-FRAG" << endl;
	  buildStatus.open("buildStatus.txt");
	  buildStatus << "crashed";
	  buildStatus.close();
	return -1;
      }
      
      return 0;
    }

  
    ///// MODE = export_energy /////
    if( mode == "export_energy" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputEnergyFile = argv[5];
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      int w = g.writeEnergy(outputEnergyFile);
      if(w != 0)
      {
	cerr << "Error during writing of the density of states in the text file " << outputEnergyFile << endl;
	return -1;
      }
      return 0;
      
    }   
    
    ///// MODE = export_transfer /////
    if( mode == "export_transfer" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputTransferFile = argv[5];
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      int w = g.writeTransfer(outputTransferFile);
      if(w != 0)
      {
	cerr << "Error during writing of the transfer integral distribution in the text file " << outputTransferFile << endl;
	return -1;
      }
      return 0;
      
    }
    
    ///// MODE = export_rate /////
    if( mode == "export_rate" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputRateFile = argv[5];
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      
      // Set up parameters to calculate rates
      std::vector<Vect3> pField;
      std::vector<double> pTemperature;
      std::vector<double> pConcentration;
      ConvergenceParameters convParam;
      readMobilityParam(inputConfigFile,pField,pTemperature,pConcentration, convParam);
      // Transform global carrier concentration in nb carrier per site
      for(int i(0) ; i < pConcentration.size() ; i++)
      {
	  pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
      }
      Vect3 field(0,0,0);
      field = pField[0] * 1e-8;
      std::cout << "Field intensity is : " << field.norm() * 1e8 << " V/cm" << std::endl;
      g.m_field = field;
      g.m_temperature = pTemperature[0];
      g.m_concentration = pConcentration[0];
      g.computeRate(maTransfer);
      
      int w = g.writeRate(outputRateFile);
      if(w != 0)
      {
	cerr << "Error during writing of the hopping rate distribution in the text file " << outputRateFile << endl;
	return -1;
      }
      return 0;
      
    }
    
    
    ///// MODE = export_geometry /////
    if( mode == "export_geometry" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputGeometryFile = argv[5];
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      int w = g.writeGeometry(outputGeometryFile);
      if(w != 0)
      {
	cerr << "Error during writing of the geometry in the text file " << outputGeometryFile << endl;
	return -1;
      }
      return 0;
      
    }
    
    ///// MODE = perco /////
    if( mode == "perco" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputPercoFile = argv[5];
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      int p = plotPercolation(g, outputPercoFile, false);
      if(p != 0)
      {
	cerr << "Error during writing of the percolation file  " << outputPercoFile << endl;
	return -1;
      }
      return 0;
      
    }

    ///// MODE = mobility /////
    if( mode == "mobility" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputMobilityFolder = argv[5];
      
      // Prepare folder to save mobility results
      outputMobilityFolder = outputMobilityFolder + "_" + std::to_string(seed);
      #ifdef _WIN32 || _WIN64
      _mkdir((std::string(_getcwd(NULL,0)) + "\\"+outputMobilityFolder).c_str());
      #else
      mkdir(outputMobilityFolder.c_str(),0777);
      #endif // linux
      outputMobilityFolder = outputMobilityFolder + "/";
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      
      std::vector<Vect3> pField;
      std::vector<double> pTemperature;
      std::vector<double> pConcentration;
      ConvergenceParameters convParam;
      readMobilityParam(inputConfigFile,pField,pTemperature,pConcentration, convParam);
      // Transform global carrier concentration in nb carrier per site
      for(int i(0) ; i < pConcentration.size() ; i++)
      {
	  pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
      }
      
      computeMobilityDependance(g, pField, pTemperature, pConcentration , outputMobilityFolder, maTransfer, convParam);
      
      return 0;
      
    }
    
    
    ///// MODE = export_current /////
    if( mode == "export_current" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string inputMobilityFolder = argv[5];
      
      inputMobilityFolder = inputMobilityFolder + "/";
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      
      std::vector<Vect3> pField;
      std::vector<double> pTemperature;
      std::vector<double> pConcentration;
      ConvergenceParameters convParam;
      readMobilityParam(inputConfigFile,pField,pTemperature,pConcentration, convParam);
      // Transform global carrier concentration in nb carrier per site
      for(int i(0) ; i < pConcentration.size() ; i++)
      {
	  pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
      }
      
      exportLocalCurrent(g, pField, pTemperature, pConcentration , inputMobilityFolder, maTransfer);
      
      return 0;
      
    }
    
    
    ///// MODE = export_graph /////
    if( mode == "export_graph" )
    {
      std::string inputPdbFile = argv[2];
      std::string inputConfigFile = argv[3];
      std::string inputGraphFile = argv[4];
      std::string outputGraphFile = argv[5];
      int mode = atoi(argv[6]);
      std::string inputMobilityFolder = argv[7];
      
      inputMobilityFolder = inputMobilityFolder + "/";
      
      Graph g = initGraph(inputPdbFile, inputConfigFile, inputGraphFile);
      g.completeGraph(seed);
      
      if(mode == 1 || mode == 2)
      {
        // Set up parameters to calculate rates
        std::vector<Vect3> pField;
        std::vector<double> pTemperature;
        std::vector<double> pConcentration;
        ConvergenceParameters convParam;
        readMobilityParam(inputConfigFile,pField,pTemperature,pConcentration, convParam);
        // Transform global carrier concentration in nb carrier per site
        for(int i(0) ; i < pConcentration.size() ; i++)
        {
  	      pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
        }
        Vect3 field(0,0,0);
        field = pField[0] * 1e-8;
        std::cout << "Field intensity is : " << field.norm() * 1e8 << " V/cm" << std::endl;
        g.m_field = field;
        g.m_temperature = pTemperature[0];
        g.m_concentration = pConcentration[0];
        g.computeRate(maTransfer); 
        
        if(mode == 2)
        {
		      std::string carrierDistributionFile(inputMobilityFolder + "distrib_F" + std::to_string(field.m_x * 1e8)  + "_" + std::to_string(field.m_y * 1e8) + "_" + std::to_string(field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt" );
		      g.readCarrierDistribution(carrierDistributionFile);
  
		      for (int n=0 ; n < g.m_site.size() ; n ++)
		      {
			      g.m_site[n].m_displacementField.resize(0);
			      for (int m(0) ; m < g.m_site[n].m_neighbor.size() ; m ++)
			      {
			        g.m_site[n].m_displacementField.push_back(g.m_site[n].m_neighborDisplacement[m]*g.m_field);
			      }
	        }
        }
        
      }     
      
      
      int w = g.writeGraphVmd(outputGraphFile, mode);
      if(w != 0)
      {
	cerr << "Error during writing of the graph in the text file " << outputGraphFile << endl;
	return -1;
      }
      return 0;
      
    }
    
    
/*    

    /// Parse the input arguments ///
    std::string mode(argv[1]);
    if( (argc < 5 && (mode != "prep" && mode != "coarse") ) || (argc < 4 && mode == "prep") || (argc < 4 && mode == "coarse") || ( argc < 7 && mode == "export_lc" ) )
    {
        cerr << "Not enough arguments !!" << '\n';
        return 0;
    }


    if( ! (mode == "prep" || mode == "perco" || mode == "mob" || mode == "export" || mode == "check" || mode == "network" || mode == "perco_plot" || mode == "system_cluster" || mode == "transfer_geom" || mode == "transfer_geomNN" || mode == "perco_plot_N" || mode == "perco_N" || mode == "coarse" || mode == "bound_distrib" || mode == "export_lc" ) )
    {
        cerr << "Mode not available !!" << '\n';
        return 0;
    }
    std::string folder;
    if( mode == "export_lc" ) { seed = atoi(argv[6]); }
    if(argc == 5 || (mode == "prep" && argc == 4) || (mode == "coarse") ) {folder = "";}
    else
    {
        folder = argv[5];
        folder = folder + "_" + std::to_string(seed);

        #ifdef _WIN32 || _WIN64
        _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
        #else
        mkdir(folder.c_str(),0777);
        #endif // linux

        folder = folder + "/";

    }

    if( mode == "prep" && argc == 5 )
    {
        folder = argv[4];
        folder = folder + "_" + std::to_string(seed);

        #ifdef _WIN32 || _WIN64
        _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
        #else
        mkdir(folder.c_str(),0777);
        #endif // linux

        folder = folder + "/";
    }

    /// Get file names
    std::string fileNameTop;
    std::string fileNameCoord;
    if(mode != "coarse")
    {
        fileNameTop = argv[2];
        fileNameCoord = argv[3];
    }
    else
    {
        fileNameCoord = argv[2];
    }
    std::string configFile;
    if( (mode != "prep") && (mode != "coarse") ) {configFile = argv[4];}
    if( mode == "coarse" ) {configFile = argv[3];}

    /// Prepare file for gorilla ///
    if(mode == "prep")
    {
        std::cout << "I will save files in " << folder << std::endl;
        prepareForGorilla(fileNameTop, fileNameCoord, folder);
    }

    /// Write bound distribution ///
    if(mode == "bound_distrib")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder, false);
        g.writeboundDistribution("boundDistrib.txt");
    }


    /// Find percolation threshold ///
    if(mode == "perco")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        Vect3 thresholdPerco;
        thresholdPerco = findPercolationThreshold(g, folder);
    }

    /// Find percolation threshold ///
    if(mode == "perco_N")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double thresholdPerco;
        thresholdPerco = findPercolationThresholdN(g, folder);
    }

    /// Export the percolation behavior curve ///
    if(mode == "perco_plot")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double threshMin(-6);
        double threshMax(2);
        double threshStep(0.1);
        plotPercolation(g, threshMin, threshMax, threshStep, folder, true);
    }

    /// Export the percolation behavior curve ///
    if(mode == "perco_plot_N")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        double threshMin(-8);
        double threshMax(3);
        double threshStep(0.01);
        plotPercolation(g, threshMin, threshMax, threshStep, folder, false);
    }


    /// Compute mobility ///
    if(mode == "mob")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        std::vector<double> pField;
        std::vector<double> pTemperature;
        std::vector<double> pConcentration;
        readMobilityParam(configFile,pField,pTemperature,pConcentration);
        // Transform global carrier concentration in nb carrier per site
        for(int i(0) ; i < pConcentration.size() ; i++)
        {
            pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
        }
        computeMobilityDependance(g, pField, pTemperature, pConcentration ,folder, maTransfer);
    }

    /// Export local current distribution ///
    if(mode == "export_lc")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        std::vector<double> pField;
        std::vector<double> pTemperature;
        std::vector<double> pConcentration;
        readMobilityParam(configFile,pField,pTemperature,pConcentration);
        // Transform global carrier concentration in nb carrier per site
        for(int i(0) ; i < pConcentration.size() ; i++)
        {
            pConcentration[i] = pConcentration[i] / (double) g.m_site.size();
        }
        exportLocalCurrent(g, pField, pTemperature, pConcentration ,folder, maTransfer);
    }


    /// Export diagonal and off-diagonal disorder ///
    if(mode == "export")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

    }


    /// Check system integrity ///
    if(mode == "check")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        cout << "Checking integrity..." << '\n';

        for(int i(0) ; i < g.m_site.size() ; i ++)
        {
            for(int j(0) ; j < g.m_site[i].m_neighborIdx.size() ; j++)
            {
                int idxNeigh(g.m_site[i].m_neighborIdx[j]);
                int id(-1);

                for(int k(0) ; k < g.m_site[idxNeigh].m_neighborIdx.size() ; k++)
                {
                    if(g.m_site[idxNeigh].m_neighborIdx[k] == i)
                    {
                        id = k;
                    }
                }

                assert(id != -1);
                if(id == -1)
                {
                    cout << "Integrity issue in the data !"  << '\n';
                    cout << "Missing neighbor"  << '\n';
                }
                else
                {
                    if( ((g.m_site[i].m_hopTo[j] - g.m_site[idxNeigh].m_hopFrom[id])/g.m_site[i].m_hopTo[j] > 1e-10) || ((g.m_site[i].m_hopFrom[j] -  g.m_site[idxNeigh].m_hopTo[id])/g.m_site[i].m_hopFrom[j] > 1e-10) )
                    {
                        cout << setprecision(20);
                        cout << "Integrity issue in the data !"  << '\n';
                        cout << "Molecules " << i << " and " << idxNeigh << '\n';
                        cout << g.m_site[i].m_hopTo[j]  << '\n';
                        cout << g.m_site[idxNeigh].m_hopFrom[id]  << '\n';
                        cout << g.m_site[i].m_hopFrom[j]  << '\n';
                        cout << g.m_site[idxNeigh].m_hopTo[id]  << '\n';
                        cout << '\n' << '\n';
                    }

                }
            }
        }

        cout << "Done." << '\n';



    }


    /// Export the transport network to visualize with VMD ///
    if(mode == "network")
    {
        string graphPdb(folder+"graph.pdb");
        string graphPsf(folder+"graph.psf");
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        double threshold(0);
        cout << "Value for clustering ? (666 is for no clustering) " << '\n';
        cin >> threshold;

        if(threshold != 666)
        {
            g.cluster(exp(threshold*log(10))/1000.0);
        }

        g.writeGraph(graphPdb,graphPsf,false);

        g.writeCurrent("occup.pdb","current.pdb","current.psf",false);
    }

    /// Export the pedot system to visualize with VMD, with molecule associated to clusterId ///
    if(mode == "system_cluster")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);

        double threshold(0);
        cout << "Value for clustering ? (666 is for no clustering) " << '\n';
        cin >> threshold;

        if(threshold != 666)
        {
            g.cluster(exp(threshold*log(10))/1000.0);
        }

        g.writeSystemClustered("system_cluster.pdb",folder);
    }

    /// Export transfer integrals and associated geometries ///
    if(mode == "transfer_geom")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        string transferGeomFile(folder+"transferGeom.txt");
	g.writeFragmentVector("orientation.tcl");
        g.writeTransferDistanceAngle(transferGeomFile,0,0, maTransfer);
        
    }

    /// Export transfer integrals and associated geometries for use with Neural Network///
    if(mode == "transfer_geomNN")
    {
        Graph g = initGraph(fileNameTop,fileNameCoord,configFile, seed, folder);
        string transferGeomFile(folder+"transferGeomNN.txt");
        g.writeTransferDistanceAngleNN(transferGeomFile,0,0, maTransfer);
        //g.writeFragmentVector("orientation.tcl");
    }

    /// Calculate morphological properties of a coarse grained MD ///
    if(mode == "coarse")
    {



    }

    */

    return 0;

}



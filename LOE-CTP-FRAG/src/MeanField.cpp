#include "MeanField.h"

using namespace std;


doubleType adjustConcentration(doubleType &pTemp, Graph &g)
{
    doubleType carrierFactor(0.0);
    if(g.m_levelType == "HOMO") {carrierFactor = -1.0;}
    if(g.m_levelType == "LUMO") {carrierFactor = 1.0;}
  
    // Adjust chemical potential to get desired carrier concentration
    doubleType cP(0); //Chemical potential (eV)
    doubleType aD,bD;
    doubleType tol(dichoPrecision); // Tolerance
    doubleType errPTmp(0);
    doubleType errP(-1e10);

    if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
    else { aD = -20, bD = 0;}
    while( (bD-aD) / 2 > tol )
    {
        pTemp = 0;
        cP = (aD + bD) / 2;

        for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / repetition ; v++)
        {
	    //CORRECT
	    //for(int cn(0) ; cn < repetition ; cn ++)
	    //{//CORRECT
	  
            for(int u(0) ; u < g.m_nbrLevels ; u++)
            {
		//CORRECT
		int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
                //int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v] + cn*(g.m_site.size() / g.m_nbrLevels / repetition));
	      
                pTemp += ( 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy - g.m_site[i].m_cP - cP) / (kB*g.m_temperature) ) ) ) ;

            //}//CORRECT
            
	    }
        }
        
        //CORRECT
        pTemp /= g.m_site.size() / repetition;
	//pTemp /= g.m_site.size();
	
        //std::cout << pTemp << "     " << g.m_concentration << '\n';
        if (pTemp < g.m_concentration ) { aD = cP;}
        else { bD = cP;}

    }

    for(int i(0) ; i < g.m_site.size() ; i++)
    {
        g.m_site[i].m_cP += cP;
        g.m_site[i].m_occup = 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy - g.m_site[i].m_cP) / (kB*g.m_temperature) ) );
        //g.m_site[i].m_occupT = g.m_site[i].m_occup;
        errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
        if(errPTmp > errP)
        {
            errP = errPTmp;

        }


    }


    return errP;
}

doubleType computeMobility(Graph &g, bool usePrevDistrib, std::string folder, ConvergenceParameters &convParam)
{
  
    /*std::cout << "Field = " << std::endl;
    std::cout << g.m_field.m_x << std::endl;
    std::cout << g.m_field.m_y << std::endl;
    std::cout << g.m_field.m_z << std::endl;*/
  
    doubleType carrierFactor(0.0);
    if(g.m_levelType == "HOMO") {carrierFactor = -1.0;}
    if(g.m_levelType == "LUMO") {carrierFactor = 1.0;}
    
    // Convergence parameters 
    doubleType accuracyMob(convParam.m_accuracyMob);
    doubleType accuracyOcc(convParam.m_accuracyOcc);
    doubleType accuracyCon(convParam.m_accuracyCon);
    int nbIterUpdate(convParam.m_nbIterUpdate);
    int nbIterMax(convParam.m_nbIterMax);
    
    doubleType relaxFactor(1.0);

    cout << "Computing mobility...";

    /// Pre-compute r.F to speed up calculation ///
    for (int i=0 ; i < g.m_site.size() ; i ++)
    {
        g.m_site[i].m_displacementField.resize(0);
        for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
        {
            g.m_site[i].m_displacementField.push_back(g.m_site[i].m_neighborDisplacement[j]*g.m_field);
	    /*double temp(g.m_site[i].m_neighborDisplacement[j]*g.m_field );
	    if(fabs(temp) != 0.0004 && fabs(temp) != 0.0008 && fabs(temp) != 0.0)
	    {
	      std::cout << "here " << temp << std::endl;
	    }*/
        }
    }
    
    /// Calculate hopping rate scaling factor ///
    std::vector<doubleType> hoppingScale;
    for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / repetition ; v++)
    {
	for(int u(0) ; u < g.m_nbrLevels ; u++)
	{
	    int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
	    doubleType sum(0);
	    
	    for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
	    {
		sum += g.m_site[i].m_hopTo[j];
	    }
	    
	    hoppingScale.push_back(sum);
	}
    }

    

    /// Init system sites occupancies (first time we compute mobility)///
    doubleType pTemp(0);
    if ( !usePrevDistrib )
    {
        //cout << '\n' << "Init system properties...";

        for (int i(0) ; i < g.m_site.size() ; i ++)
        {
            // Set occupation probability
            g.m_site[i].m_occup = 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy ) / (kB*g.m_temperature) ) );
	    
            g.m_site[i].m_occupT = g.m_site[i].m_occup;
            pTemp += g.m_site[i].m_occup / g.m_site.size();
	  
	}

        //cout << "ok" << endl;
    }

    /// Init system sites occupancies (using previous occupancy distribution as the starting point)///
    if ( usePrevDistrib )
    {
        //cout << '\n' << "Init system properties...";

        for (int i(0) ; i < g.m_site.size() ; i ++)
        {
            // Set occupation probability
            g.m_site[i].m_occup = g.m_site[i].m_occupT;
            pTemp += g.m_site[i].m_occup / g.m_site.size();
        }

        //cout << "ok" << endl;
        //cout << "pTemp = " << pTemp << endl;
    }

    /// Adjust chemical potential to get the desired carrier concentration (first time we compute mobility)///
    if( !usePrevDistrib )
    {
        //cout << "Adjust chemical potential..." << endl;
        doubleType cP(0); //Chemical potential (eV)
        doubleType aD,bD;
        doubleType tol(dichoPrecision); // Tolerance

        if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
        else { aD = -20, bD = 0;}

        while( (bD-aD) / (doubleType) 2.0 > tol )
        {
            pTemp = 0;

            cP = (aD + bD) / (doubleType) 2.0;

            for (int i(0) ; i < g.m_site.size() ; i ++)
            {


                g.m_site[i].m_occup = 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy - cP) / (kB*g.m_temperature) ) );

                g.m_site[i].m_occupT = g.m_site[i].m_occup;
                g.m_site[i].m_cP = cP;
                pTemp += g.m_site[i].m_occup / g.m_site.size();
            }

            if (pTemp < g.m_concentration ) { aD = cP;}
            else { bD = cP;}
            
        }
        
    }
    

    /// Adjust chemical potential to get the desired carrier concentration (using previous occupancy distribution as the starting point)///
    if( usePrevDistrib )
    {
        doubleType cP(0); //Chemical potential (eV)
        doubleType aD,bD;
        doubleType tol(dichoPrecision); // Tolerance

        if(fabs(g.m_concentration-pTemp)/g.m_concentration > accuracyCon)
        {
            if( pTemp < g.m_concentration ) { aD = 0; bD = 20;}
            else { aD = -20, bD = 0;}

            while( (bD-aD) / (doubleType) 2.0 > tol )
            {
                pTemp = 0;

                cP = (aD + bD) / (doubleType) 2.0;

                for (int i(0) ; i < g.m_site.size() ; i ++)
                {
                    pTemp += ( 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy - g.m_site[i].m_cP - cP) / (kB*g.m_temperature) ) ) ) / g.m_site.size();
                }

                if (pTemp < g.m_concentration ) { aD = cP;}
                else { bD = cP;}
            }

        }
        for(int i(0) ; i < g.m_site.size() ; i++)
        {
            g.m_site[i].m_cP += cP;
            g.m_site[i].m_occup = 1.0 / ( 1 + exp( (carrierFactor*g.m_site[i].m_energy - g.m_site[i].m_cP) / (kB*g.m_temperature) ) );
            g.m_site[i].m_occupT = g.m_site[i].m_occup;

        }
    }


    /// Solve equilibrium probability distribution ///
    // Temp variables
    doubleType err(1e6);
    doubleType mobilityF(0);
    doubleType errOccup(1e6);
    int nbIterTot(0);


    std::ofstream errFlux(folder + "convergence_F" + std::to_string(g.m_field.m_x * 1e8)  + "_" + std::to_string(g.m_field.m_y * 1e8) + "_" + std::to_string(g.m_field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt");
    errFlux << std::setw(outputLength) << "Number of iterations" << std::setw(outputLength) << "Error on mobility" << std::setw(outputLength)  << "Error on occupancy" << std::setw(outputLength)  << "Mobility" << std::setw(outputLength) << "MeanOccup" << std::setw(outputLength) << "Solution check" << endl;

    do
    {
        pTemp = 0;
        doubleType errPTmp(0);
        doubleType errP(-1e10);

        // Solve master equation
	int countSiteCalculated(0);
        for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / repetition ; v++)
        {
	  
	  	      	  //for(int cn(0) ; cn < repetition  ;cn++) // CORRECT
		  //{ //CORRECT
	  
            for(int u(0) ; u < g.m_nbrLevels ; u++)
            {
	      
	       //CORRECT
                int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
		//int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v] + cn*(g.m_site.size() / g.m_nbrLevels / repetition) );
                
		doubleType term1(0);
                doubleType term2(0);
		
                for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
                {
                        int k(g.m_site[i].m_neighborIdx[j]);

                        // Direct update of occupancy during iterations, no hopping rate scaling
			//TRY
                        //term1 += g.m_site[k].m_occup * g.m_site[i].m_hopFrom[j];
                        //term2 += (1-g.m_site[k].m_occup) * g.m_site[i].m_hopTo[j];
			//term1 += g.m_site[k].m_occup * 1.0;
                        //term2 += (1-g.m_site[k].m_occup) * 1.0;
			
			// Method with hopping rate scaling
			term1 += g.m_site[k].m_occup * (g.m_site[i].m_hopFrom[j] / hoppingScale[countSiteCalculated]);
                        term2 += g.m_site[k].m_occup * ( (g.m_site[i].m_hopTo[j]/ hoppingScale[countSiteCalculated]) - (g.m_site[i].m_hopFrom[j] / hoppingScale[countSiteCalculated]) );


                        // No update during iterations
                        //term1 += g.m_site[k].m_occupT * g.m_site[i].m_hopFrom[j];
                        //term2 += (1-g.m_site[k].m_occupT) * g.m_site[i].m_hopTo[j];
                }
		// Direct update during iterations
                g.m_site[i].m_occupT = g.m_site[i].m_occup;
		
		// Method without hopping scaling
                g.m_site[i].m_occup =  relaxFactor * (term1 / (term1 + term2)) + (1.0 - relaxFactor) * g.m_site[i].m_occupT;

		// Method with hopping rate scaling
		g.m_site[i].m_occup = relaxFactor * term1 / (1.0 - term2) + (1.0-relaxFactor) * g.m_site[i].m_occupT;
		
                // Calculate the error on occupation
                errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
                if(errPTmp > errP) {errP = errPTmp;}
                
                // Direct update uring iterations
                //g.m_site[i].m_occupT = g.m_site[i].m_occup;

                g.m_site[i].m_cP = -( log(1.0 / g.m_site[i].m_occup - 1) * kB * g.m_temperature - (carrierFactor*g.m_site[i].m_energy) );
                pTemp += g.m_site[i].m_occup;
		
		//CORRECT
		
                for(int w(1) ; w <repetition ; w ++)
                {
                    int q(i+w*(g.m_site.size() / g.m_nbrLevels / repetition));
                    g.m_site[q].m_occup = g.m_site[i].m_occup;
                    g.m_site[q].m_cP = g.m_site[i].m_cP;
		    
		    // Direct update uring iterations
                    g.m_site[q].m_occupT = g.m_site[q].m_occup;
                }
                
                
                countSiteCalculated ++;
            }
	  //}//CORRECT
        }
         
        errOccup = errP;
        
	//CORRECT
	pTemp /= g.m_site.size() / repetition;
	//pTemp /= g.m_site.size();
	
	// No update during iterations
        /*for(int h(0) ; h < g.m_site.size() ; h++)
        {
            g.m_site[h].m_occupT = g.m_site[h].m_occup;
        }*/


        // Adjust concentration
        if(fabs(g.m_concentration-pTemp) / g.m_concentration > accuracyCon)
        {
            errOccup = adjustConcentration(pTemp,g);
        }

        nbIterTot ++;

        if(nbIterTot % nbIterUpdate == 0)
        {
            // Adjust concentration
            //adjustConcentration(pTemp,g);

            // Compute mobility
            doubleType mobility(0);
            errPTmp = 0;
            errP = -1e10;
            for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / repetition ; v++)
            {
	      	  	      	  //for(int cn(0) ; cn < repetition  ;cn++) // CORRECT
		  //{ //CORRECT
	      
                for(int u(0) ; u < g.m_nbrLevels ; u++)
                {
		    //CORRECT
                    int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
		    //int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v] + cn*(g.m_site.size() / g.m_nbrLevels / repetition) );
                    
		    for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
                    {
                            int k(g.m_site[i].m_neighborIdx[j]);
			    //TRY
                            mobility += g.m_site[i].m_occup * ( 1 - g.m_site[k].m_occup ) * g.m_site[i].m_hopTo[j] * g.m_site[i].m_displacementField[j];
			    //mobility += g.m_site[i].m_occup * ( 1 - g.m_site[k].m_occup ) * 1.0 * g.m_site[i].m_displacementField[j];
                    }

                    // Calculate the error on occupation
                    //errPTmp = fabs(g.m_site[i].m_occup - g.m_site[i].m_occupT)/g.m_site[i].m_occupT;
                    //if(errPTmp > errP) {errP = errPTmp;}

                }
                
		  //}//CORRECT
            }
            //errOccup = errP;
            
            //CORRECT
            mobility /= pTemp * (g.m_site.size() / (doubleType) repetition) * ( g.m_field.norm() * g.m_field.norm() );
            //mobility /= pTemp * g.m_site.size() * ( g.m_field.norm() * g.m_field.norm() );
	    
	    
            //mobility /= pTemp * (g.m_site.size() / (doubleType) repetition) * ( g.m_field.norm() );
            err = fabs(mobility - mobilityF) / fabs(mobilityF) / (doubleType) nbIterUpdate;
            if (mobilityF == 0) {err = 1.0;}
            mobilityF = mobility;
            errFlux << std::setprecision(outputPrecision) << std::setw(outputLength) << nbIterTot << std::setw(outputLength) << err  << std::setw(outputLength) << errOccup << std::setw(outputLength) << mobilityF << std::setw(outputLength) << pTemp << std::setw(outputLength) << checkTransportConvergence(g) << endl ;
	    

        }

    }while ( (err > accuracyMob || errOccup > accuracyOcc) && nbIterTot <= nbIterMax );

    errFlux << std::setprecision(outputPrecision) << std::setw(outputLength) << nbIterTot << std::setw(outputLength) << err  << std::setw(outputLength) << errOccup << std::setw(outputLength) << mobilityF << std::setw(outputLength) << pTemp << std::setw(outputLength) << checkTransportConvergence(g) << endl ;
    errFlux.close();
    cout << "Done (Mobility is " << mobilityF << ")" << endl;
    return mobilityF;
}


Vect3 findPercolationThreshold(Graph &g, std::string folder)
{
    // All the values in unit log10(H*1000) H in eV
    double threshMin(-30.0);
    double threshMax(5.0);
    double threshErr(0.01);
    Vect3 threshVect;
    Vect3 clusterIdVect;

    for (int i(0) ; i < 3 ; i++)
    {

        double aD(threshMin);
        double bD(threshMax);
        double thresh(0);
        int clusterId;
        bool isInf(false);

        while( (bD-aD) / 2.0 > threshErr || !isInf)
        {

            thresh = (bD+aD) / 2.0;
            std::cout << "threshold: " << thresh << '\n';
            double threshMod = (exp(thresh*log(10))/1000.0);

            g.cluster(threshMod, 1);
            g.percolation();

            double clusterMax(0);

            for(int j(0) ; j < g.m_clusterNb ; j ++)
            {
                double cSize(0);
                switch(i)
                {
                    case 0 : cSize = g.m_clusterDim[j].m_x; break;
                    case 1 : cSize = g.m_clusterDim[j].m_y; break;
                    case 2 : cSize = g.m_clusterDim[j].m_z; break;

                }

                if(cSize > clusterMax)
                {
                    clusterMax = cSize;
                    clusterId = j;
                }
            }

            cout << clusterMax << '\n';

            if ( clusterMax == 1000 ) { aD = thresh; isInf = true;}
            else { bD = thresh; isInf = false;}
        }

        switch(i)
        {
            case 0 : threshVect.m_x = thresh; clusterIdVect.m_x = clusterId; break;
            case 1 : threshVect.m_y = thresh; clusterIdVect.m_y = clusterId; break;
            case 2 : threshVect.m_z = thresh; clusterIdVect.m_z = clusterId; break;

        }

        // Save result of the percolated graph
        string graphPdb(folder+"graph"+std::to_string(i) +".pdb");
        string graphPsf(folder+"graph"+std::to_string(i) +".psf");
        string systemClustered("system_cluster"+std::to_string(i)+".pdb");
        g.writeGraph(graphPdb,graphPsf,false);
        g.writeSystemClustered(systemClustered,folder);
        g.writeCurrent(folder + "occup"+std::to_string(i) +".pdb",folder + "current"+std::to_string(i) +".pdb",folder + "current"+std::to_string(i) +".psf",false);

    }

    // Write the result in a file
    ofstream os(folder + "percolationThreshold.txt");
    if(os.is_open())
    {
        os << "X" << "  " << "Y" << "   " << "Z" << '\n';
        os << setprecision(15) << threshVect.m_x << " " << threshVect.m_y << "    " << threshVect.m_z << '\n';
        os << clusterIdVect.m_x << "    " << clusterIdVect.m_y << "     " << clusterIdVect.m_z << '\n';
        os.close();
    }
    else
    {
        std::cerr << "Can't open percolation threshold file !" << '\n';
    }

    std::cout << "Percolation thresholds: " << threshVect.m_x << "  " << threshVect.m_y << "    " << threshVect.m_z << '\n';
    return threshVect;

}

double findPercolationThresholdN(Graph &g, std::string folder)
{
    // All the values in unit log10(H*1000) H in eV
    double threshMin(-30.0);
    double threshMax(5.0);
    double threshErr(0.01);

    double aD(threshMin);
    double bD(threshMax);
    double thresh(0);
    int clusterId;
    bool isInf(false);

    while( (bD-aD) / 2.0 > threshErr || !isInf)
    {

        thresh = (bD+aD) / 2.0;
        std::cout << "threshold: " << thresh << '\n';
        double threshMod = (exp(thresh*log(10))/1000.0);

        g.cluster(threshMod, 1);
        double clusterMax(0);

        for(int j(0) ; j < g.m_clusterNb ; j ++)
        {
            if(g.m_clusterSiteNb[j] > clusterMax)
            {
                clusterMax = g.m_clusterSiteNb[j];
                clusterId = j;
            }
        }

        cout << clusterMax << '\n';

        if ( clusterMax >= g.m_site.size() / 2.0 ) { aD = thresh; isInf = true;}
        else { bD = thresh; isInf = false;}
    }


    // Save result of the percolated graph
    string graphPdb(folder+"graph" +".pdb");
    string graphPsf(folder+"graph" +".psf");
    string systemClustered("system_cluster.pdb");
    g.writeGraph(graphPdb,graphPsf,false);
    g.writeSystemClustered(systemClustered,folder);
    g.writeCurrent(folder + "occup" +".pdb",folder + "current"+".pdb",folder + "current"+".psf",false);

    // Write the result in a file
    ofstream os(folder + "percolationThreshold.txt");
    if(os.is_open())
    {
        os << "Threshold" << "  " << "ID" << '\n';
        os << setprecision(15) << thresh << " " << clusterId << '\n';
        os.close();
    }
    else
    {
        std::cerr << "Can't open percolation threshold file !" << '\n';
    }

    std::cout << "Percolation threshold: " << thresh << "   cluster Id: " << clusterId << '\n';
    return thresh;
}



void prepareSystem(std::string inputGroFile, std::string inputConfigFile, std::string outputPdbFile, int generateGaussian)
{
    /// Extend and convert into gromacs a set of pdb file in a folder ///
    vector<string> coordFile;
    coordFile.push_back(inputGroFile);

    for (int i(0) ; i < coordFile.size() ; i++)
    {

        // Fragment the system with the fragment.py python script
        //string command("python3 /proj/polymer/users/Programs/Mobility/LOE-CTP-FRAG/Fragmentation/fragment.py " + coordFile[i]);
        // CUBIC mod
	string command("fragment.py " + coordFile[i]);
	//string command("fragment_cubic.py " + coordFile[i]);
        //system("module load python3/3.6.1");
        system(command.c_str());
	string filePdb = coordFile[i].substr(0, coordFile[i].size()-4) + ".pdb";
	remove(filePdb.c_str());
        string filePdbFragmented = coordFile[i].substr(0, coordFile[i].size()-4) + "_fragmented.pdb";
        coordFile[i] = filePdbFragmented;


        string filePdbExtended = coordFile[i].substr(0, coordFile[i].size()-4) + "_big.pdb";
        //string filePdbExtendedParsed = coordFile[i].substr(0, coordFile[i].size()-4) + "_big_parsed.pdb";
        string filePdbExtendedParsed = outputPdbFile;
	extendPdb(coordFile[i]);
        parseBigPdb(filePdbExtended, filePdbExtendedParsed, filePdbFragmented);
        //string topExtended =  fileNameTop.substr(0, fileNameTop.size()-4) + "_big.top";
        //extendTop(fileNameTop);

        string fileNameGro( coordFile[i].substr(0, coordFile[i].size()-4) + ".gro");
        Graph g;
        //g.m_fileNameTop = topExtended;
        g.m_fileNameTop = "";
        g.m_fileNameCoord = filePdbExtendedParsed;
        g.initSystem();
        g.segmentSystem();
        //g.m_system.writeSystemGro(fileNameGro);
	string fileNameOrient("orientation.tcl");
	g.writeFragmentVector(fileNameOrient);

        // Remove all the unused intermediate file
        remove("extendPdb.tcl");
        remove(filePdbExtended.c_str());
        string filePdbExtendedT = filePdbExtended.substr(0, filePdbExtended.size()-4) + "T.pdb";
        remove(filePdbExtendedT.c_str());
	remove("Molefacture_tmpmol.xbgf");
	remove(filePdbFragmented.c_str());
	remove("saturated.pdb");
	remove("sat.tcl");

        // Displace files in folder
        //rename(filePdbExtendedParsed.c_str(), (folder + filePdbExtendedParsed).c_str());
        //rename(topExtended.c_str(), (folder + topExtended).c_str());
        //rename(fileNameGro.c_str(), (folder + fileNameGro).c_str());
	//rename(fileNameOrient.c_str(), (folder + fileNameOrient).c_str());
	
	// Read config file so that we can create links between conjugated segment 
	readConfigParam(inputConfigFile, g);
	
	// Create links in the system
	g.findNeighbors();
	
	// Create gaussian com files
	if(generateGaussian)
	{
	  g.generateGaussianCom();
	}

    }
}


void prepareSystemGaussian(std::string inputPdbFile, std::string inputConfigFile)
{
    Graph g;
    g.m_fileNameTop = "";
    g.m_fileNameCoord = inputPdbFile;
    g.initSystem();
    g.segmentSystem();
    
    // Read config file so that we can create links between conjugated segment 
    readConfigParam(inputConfigFile, g);
    
    // Create links in the system
    g.findNeighbors();
    
    // Create gaussin com files
    g.generateGaussianCom();

}


int plotPercolation(Graph &g, std::string percolationFile, bool dimCluster)
{

    /// Threshold curve parameters ///
    // All the values in unit ln(H), H in eV
    double threshMin(g.m_Hmin);
    double threshMax(g.m_Hmax);
    double threshStep(g.m_dH);
    std::vector<double> pThreshold;
    double thresh(threshMin-threshStep);
    while(thresh < threshMax)
    {
        thresh += threshStep;
        pThreshold.push_back(exp(thresh)*exp(thresh));
    }


    /// Compute the threshold behavior ///
    std::ofstream oFlux(percolationFile);

    if(oFlux)
    {
        if(dimCluster)
        {
            oFlux << "Threshold" << " " << "clusterMaxNorm" << "   " << "clusterMaxXYZ" << " " << "clusterMaxX" << "   " << "clusterMaxY" << "   " << "clusterMaxZ" << "    " << "clusterMaxNb" << std::endl;
        }
        else
        {
            oFlux << "Threshold" << "    " << "clusterMaxNb" << std::endl;
        }


        for (int i(0) ; i < pThreshold.size() ; i ++)
        {	
            g.cluster(pThreshold[i], 1);

            if(dimCluster)
            {
                g.percolation();
            }

            double clusterMaxNorm(0);
            double clusterMaxX(0);
            double clusterMaxY(0);
            double clusterMaxZ(0);
            double clusterMaxXYZ(0);
            double clusterMaxNb(0);

            for(int i(0) ; i < g.m_clusterNb ; i ++)
            {
                if(dimCluster)
                {
                    double cSizeXYZ = max(max(g.m_clusterDim[i].m_x,g.m_clusterDim[i].m_y) , g.m_clusterDim[i].m_z);
                    if(cSizeXYZ > clusterMaxXYZ)
                    {
                        clusterMaxXYZ = cSizeXYZ;
                    }
                    double cSizeX = g.m_clusterDim[i].m_x;
                    if(cSizeX > clusterMaxX)
                    {
                        clusterMaxX = cSizeX;
                    }
                    double cSizeY = g.m_clusterDim[i].m_y;
                    if(cSizeY > clusterMaxY)
                    {
                        clusterMaxY = cSizeY;
                    }
                    double cSizeZ = g.m_clusterDim[i].m_z;
                    if(cSizeZ > clusterMaxZ)
                    {
                        clusterMaxZ = cSizeZ;
                    }

                    double cSizeNorm = g.m_clusterDim[i].norm();
                    if(cSizeNorm > clusterMaxNorm)
                    {
                        clusterMaxNorm = cSizeNorm;
                    }
                }

                if(g.m_clusterSiteNb[i] > clusterMaxNb)
                {
                    clusterMaxNb = g.m_clusterSiteNb[i];
                }
            }
            
            clusterMaxNb /= g.m_site.size();

            if(dimCluster)
            {
                oFlux << log(sqrt(pThreshold[i])) << " " << clusterMaxNorm << "   " << clusterMaxXYZ << " " << clusterMaxX << "   " << clusterMaxY << "   " << clusterMaxZ << "  " << clusterMaxNb << std::endl;
            }
            else
            {
                oFlux << log(sqrt(pThreshold[i])) << " " << clusterMaxNb << std::endl;
            }
        }

    oFlux.close();
    }

    else
    {
        cerr << "Can't open file " << percolationFile << std::endl;
	return -1;
    }
    
    return 0;
}


void computeMobilityDependance(Graph &g, std::vector<Vect3> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder, RateType rate, ConvergenceParameters &convParam)
{
    /// Compute the mobility as a function of field///
    std::ofstream oFlux(folder + "Mobility_dependance.txt");
    if(oFlux)
    {
	oFlux << std::setw(outputLength) << "field(V/cm)" << std::setw(outputLength) << "temperature(K)" << std::setw(outputLength) << "concentration(p/site)" << std::setw(outputLength) << "Mobility(a.u.)" << std::endl;

	for (int i(0) ; i < pField.size() ; i ++)
	{
	    for(int j(0) ; j < pTemperature.size() ; j ++)
	    {
		for (int k(0) ; k < pConcentration.size() ; k ++)
		{
		    // Set up the field
		    Vect3 field(0,0,0);
		    field = pField[i] * 1e-8;
		    std::cout << "Field intensity is : " << field.norm() * 1e8 << " V/cm" << std::endl;
		    g.m_field = field;

		    // Set up the temperature and concentration
		    g.m_temperature = pTemperature[j];
		    g.m_concentration = pConcentration[k];

		    // Compute rate according to these parameters
		    g.computeRate(rate);
		    //g.computeRate(maPheno);

		    // Compute mobility and save results
		    //bool usePrevDistrib(i!=0);
		    bool usePrevDistrib(false);
		    doubleType mob = computeMobility(g,usePrevDistrib, folder, convParam);
		    oFlux << std::setprecision(outputPrecision/2) << std::setw(outputLength/3) << field.m_x * 1e8 << std::setw(outputLength/3) << field.m_y * 1e8 << std::setw(outputLength/3) << field.m_z * 1e8 << std::setw(outputLength) << g.m_temperature << std::setprecision(outputPrecision/2) << std::setw(outputLength) << g.m_concentration << std::setw(outputLength) << mob << std::endl ;

		    // Save carrier distribution
		    ofstream occupFlux(folder + "distrib_F" + std::to_string(field.m_x * 1e8)  + "_" + std::to_string(field.m_y * 1e8) + "_" + std::to_string(field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt" );
		    for (int h(0) ; h < g.m_site.size() ; h++)
		    {
			occupFlux << std::setw(outputLength) << g.m_site[h].m_pos.m_x << std::setw(outputLength) << g.m_site[h].m_pos.m_y << std::setw(outputLength) << g.m_site[h].m_pos.m_z << std::setw(outputLength) << g.m_site[h].m_level << std::setw(outputLength) << g.m_site[h].m_occup << std::endl;
		    }
		    occupFlux.close();

		}
	    }
	}
    oFlux.close();
    }

    else
    {
        cerr << "Can't open file" << '\n';
    }

}


int buildGraph(std::string inputPdbFile, std::string inputConfigFile, std::string outputGraphFile, int useZindo)
{
      /// Create a graph and init its parameters ///
    Graph g;
    g.m_fileNameCoord = inputPdbFile;
    readConfigParam(inputConfigFile, g);
    
    /// Init the system ///
    g.initSystem();
    g.segmentSystem();
    g.findNeighbors();


    /// Output some stats on the linking in the graph ///
    int nbLinkMax(-1e6);
    int nbLinkMin(1e6);
    double nbLinkmean(0);
    for (int i(0) ; i < g.m_site.size() ; i ++)
    {
	int nbLink(g.m_site[i].m_neighborIdx.size());
	if(nbLink > nbLinkMax) {nbLinkMax = nbLink;}
	if(nbLink < nbLinkMin) {nbLinkMin = nbLink;}
	nbLinkmean += nbLink;
    }
    nbLinkmean /= g.m_site.size();
    std::cout << "Max number of links: " << nbLinkMax << "  Min number of links: " << nbLinkMin << "    Mean number of links: " << nbLinkmean << std::endl;
    if(nbLinkMin < 2) 
    {
      std::cerr << "Minimum number of links is < 2, not enought connectivity, consider increasing the distance threshold! Aborting... " << std::endl;
      return -1;
      
    }
    
    /// Check that we have full connectivity in the system ///
    cout << "Check graph connectivity... ";
    g.cluster(1e-70, 0);
    if(g.m_clusterNb > 2)
    {
	cerr << "No full connectivity in the system... consider increasing distance threshold for connection! Aborting... " << endl;
	return -1;
    }
    
    
    
    
    /// Calculate energies and transfer integrals, output some stats ///
    if(useZindo)
    {
      g.getEnergyAndTransfer();
    }
    else
    {
      g.getMillerAbrahams();
    }
    cout << "Min and Max transfer integral : " << '\n';
    int cTranfer(0);
    for (int x(0) ; x < g.m_nbrLevels2Read ; x ++)
    {
	for(int y(0) ; y < g.m_nbrLevels2Read ; y ++)
	{
	    cout << "Level 1:  " << x << "  Level 2:  " << y << "  Min transfer integral:  " << g.m_minTransfer[cTranfer] << "  Max transfer integral:  " << g.m_maxTransfer[cTranfer] << '\n';
	    cTranfer ++;
	}
    }


    /// Save the graph with linking, energies and transfer integrals to be loaded later on ///
    g.parseTransfer();
    int saveGraph = g.saveGraph(outputGraphFile);
    if(saveGraph != 0) {return -2;}
    
    return 0;

}

Graph initGraph(std::string inputPdbFile, std::string inputConfigFile, std::string inputGraphFile)
{
    
    std::cout << "Initializing graph... " << std::endl;
  
    /// Create a graph and init its parameters ///
    Graph g;
    g.m_fileNameCoord = inputPdbFile;
    readConfigParam(inputConfigFile, g);
    
    /// Init the system ///
    g.initSystem();
    g.segmentSystem();
    g.findNeighbors();


    /// Output some stats on the linking in the graph ///
    int nbLinkMax(-1e6);
    int nbLinkMin(1e6);
    double nbLinkmean(0);
    
    
    // DEBUG
    //std::ofstream myfile;
    //myfile.open ("neigh.txt");
    
    for (int i(0) ; i < g.m_site.size() ; i ++)
    {
	int nbLink(g.m_site[i].m_neighborIdx.size());
	if(nbLink > nbLinkMax) {nbLinkMax = nbLink;}
	if(nbLink < nbLinkMin) {nbLinkMin = nbLink;}
	nbLinkmean += nbLink;
	//std::cout << "site id " << i << "  nb links " << nbLink << std::endl;
	//myfile << i << " " << g.m_site[i].m_pos.m_x << " " << g.m_site[i].m_pos.m_y << " " << g.m_site[i].m_pos.m_z << " " << g.m_site[i].m_neighborIdx.size() << std::endl;  
    }
    //myfile.close();
    
    nbLinkmean /= g.m_site.size();
    std::cout << "Max number of links: " << nbLinkMax << "  Min number of links: " << nbLinkMin << "    Mean number of links: " << nbLinkmean << std::endl;
    if(nbLinkMin < 2) 
    {
      std::cerr << "Minimum number of links is < 2, not enought connectivity, consider increasing the distance threshold! Aborting... " << std::endl;
      
    }
    
    /// Check that we have full connectivity in the system ///
    cout << "Check graph connectivity... ";
    g.cluster(1e-70, 0);
    if(g.m_clusterNb > 2)
    {
	cerr << "No full connectivity in the system... consider increasing distance threshold for connection! Aborting... " << endl;
    }
    
    
    /// Read energies and transfer integrals from graph file, output some stats ///
    int loadGraph = g.loadGraph(inputGraphFile);
    
    std::cout << "Initializing done!" << std::endl;
    return g;

}

void readConfigParam(std::string configFile, Graph &g)
{
    ifstream is(configFile);
    if(is.is_open())
    {
        string line;
        double param;
        getline(is,line);

        /// Graphic settings for VMD site orientation vectors ///
        is >> line >> param;
        g.m_vectLength = param;
        is >> line >> param;
        g.m_vectRadius = param;
        is >> line >> param;
        g.m_vectResol = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// Edges threshold parameters ///
        is >> line >> param;
        g.m_dThresholdInter = param;
        is >> line >> param;
        g.m_dThresholdIntra = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// Miller Abrahams parameters ///
        is >> line >> param;
        g.m_w0_intra = param;
        is >> line >> param;
        g.m_w0_inter = param;
        is >> line >> param;
        g.m_w0_level = param;
        is >> line >> param;
        g.m_loc_inter = param;
        is >> line >> param;
        g.m_loc_intra = param;
        is >> line >> param;
        g.m_H_intra = param;
        is >> line >> param;
        g.m_H_level = param;
        getline(is,line);
        getline(is,line);
        getline(is,line);

        /// MO levels information ///
        is >> line >> param;
        g.m_nbrLevels = param;
        is >> line >> line;
        g.m_levelType = line;
        is >> line >> param;
        g.m_nbrLevels2Read = param;
        is >> line >> param;
        g.m_dosBroadening = param;
        is >> line;
        for(int p(0) ; p < g.m_nbrLevels2Read ; p++)
        {
            is >> param;
            g.m_levelEnergies.push_back(param);
        }
        is >> line >> param;
	g.m_nbProcSharedGaussian = param;
	getline(is,line);
        getline(is,line);
        getline(is,line);
	
	
	/// Percolation curve parameters ///
	is >> line >> param;
        g.m_Hmin = param;
        is >> line >> param;
        g.m_Hmax = param;
        is >> line >> param;
        g.m_dH = param;
        
        is.close();

    }
    else
    {
        cerr << "Can't open config file !" << '\n';
        assert(is.is_open());
    }  
}

void readMobilityParam(std::string configFile, std::vector<Vect3> &pField, std::vector<double> &pTemperature, std::vector<double> &pConcentration, ConvergenceParameters &pConvergence)
{
    ifstream is(configFile);
    if(is.is_open())
    {
        string line;
        double param;
        int paramN;
        for (int i(0) ; i < 32 ; i ++)
        {
            getline(is,line);
        }
        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
	    Vect3 field;
            is >> param; field.m_x = param;
	    is >> param; field.m_y = param;
	    is >> param; field.m_z = param;
            pField.push_back(field);

        }

        getline(is,line);
        getline(is,line);
        getline(is,line);

        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
            is >> param;
            pTemperature.push_back(param);
        }

        getline(is,line);
        getline(is,line);
        getline(is,line);
        is >> paramN;
        for(int i(0) ; i < paramN ; i ++)
        {
            is >> param;
            pConcentration.push_back(param);
        }
    
        getline(is,line);
        getline(is,line);
        getline(is,line);
        is >> line>> param; pConvergence.m_accuracyMob = param;
	is >> line>> param; pConvergence.m_accuracyOcc = param;
	is >> line>> param; pConvergence.m_accuracyCon = param;
	is >> line>> param; pConvergence.m_nbIterUpdate = (int) param;
	is >> line>> param; pConvergence.m_nbIterMax = (int) param;

        is.close();

    }
    else
    {
        std::cerr << "Can't open config file!" << '\n';
    }

}


/// Export local current distribution ///
void exportLocalCurrent(Graph &g,  std::vector<Vect3> pField, std::vector<double> pTemperature, std::vector<double> pConcentration, std::string folder, RateType rate)
{
    /// Compute the local current distribution as a function of field///
    for (int i(0) ; i < pField.size() ; i ++)
    {
	for(int j(0) ; j < pTemperature.size() ; j ++)
	{
	    for (int k(0) ; k < pConcentration.size() ; k ++)
	    {
		// Set up the field
		Vect3 field(0,0,0);
		field = pField[i] * 1e-8;
		std::cout << "Field intensity is : " << field.norm() * 1e8 << " V/cm" << std::endl;
		g.m_field = field;

		// Set up the temperature and concentration
		g.m_temperature = pTemperature[j];
		g.m_concentration = pConcentration[k];

		// Compute rate according to these parameters
		g.computeRate(rate);
		//g.computeRate(maPheno);

		// Read local currents
		std::string carrierDistributionFile(folder + "distrib_F" + std::to_string(field.m_x * 1e8)  + "_" + std::to_string(field.m_y * 1e8) + "_" + std::to_string(field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt" );
		g.readCarrierDistribution(carrierDistributionFile);

		/// Pre-compute r.F to speed up calculation ///
		for (int n=0 ; n < g.m_site.size() ; n ++)
		{
			g.m_site[n].m_displacementField.resize(0);
			for (int m(0) ; m < g.m_site[n].m_neighbor.size() ; m ++)
			{
			    g.m_site[n].m_displacementField.push_back(g.m_site[n].m_neighborDisplacement[m]*g.m_field);
			}
		}

		// Compute local current and save results
		//std::string localCurrentFile(folder + "lc_F" + std::to_string(field.m_x * 1e8)  + "_" + std::to_string(field.m_y * 1e8) + "_" + std::to_string(field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".txt" );
		//std::ofstream lcFlux(localCurrentFile);
		std::string localCurrentFile(folder + "lc_F" + std::to_string(field.m_x * 1e8)  + "_" + std::to_string(field.m_y * 1e8) + "_" + std::to_string(field.m_z * 1e8) + "_T" + std::to_string(g.m_temperature) + "_P" + std::to_string(g.m_concentration) + ".bin" );
		std::ofstream lcFlux(localCurrentFile, std::ios::out | std::ios::binary);

		if(lcFlux)
		{
		    double pTemp(0);
		    for (int n(0) ; n < g.m_site.size() ; n++)
		    {
			//for(int m(0) ; m < g.m_nbrLevels ; m++)
			//{
			    int s(n);
			    pTemp += g.m_site[s].m_occup;

			    for (int l(0) ; l < g.m_site[s].m_neighbor.size() ; l ++)
			    {

				    int q(g.m_site[s].m_neighborIdx[l]);
				    //double localCurrent = g.m_site[s].m_occup * ( 1 - g.m_site[q].m_occup ) * g.m_site[s].m_hopTo[l] * g.m_site[s].m_displacementField[l];
				    //double localCurrent = g.m_site[s].m_occup * ( 1 - g.m_site[q].m_occup ) * g.m_site[s].m_hopTo[l] * g.m_site[s].m_neighborDisplacement[l].norm();
				    double localCurrent = g.m_site[s].m_occup * ( 1 - g.m_site[q].m_occup ) * g.m_site[s].m_hopTo[l];
				    //localCurrent /= pTemp * ( g.m_field.norm() );

				    //lcFlux << std::min(s,q) << "	" << std::max(s,q) << "	" << g.m_site[s].m_pos.m_x << "    " << g.m_site[s].m_pos.m_y << "    " << g.m_site[s].m_pos.m_z << "    " << g.m_site[q].m_pos.m_x << "    " << g.m_site[q].m_pos.m_y << "    " << g.m_site[q].m_pos.m_z << "    " << localCurrent << std::endl;
				    //double mI(std::min(s,q));
				    //double mA(std::max(s,q));
				    double mI(s);
				    double mA(q);
				    lcFlux.write((char*)&mI, sizeof(double));
				    lcFlux.write((char*)&mA, sizeof(double));
				    lcFlux.write((char*)&g.m_site[s].m_pos.m_x, sizeof(double));
				    lcFlux.write((char*)&g.m_site[s].m_pos.m_y, sizeof(double));
				    lcFlux.write((char*)&g.m_site[s].m_pos.m_z, sizeof(double));
				    lcFlux.write((char*)&g.m_site[q].m_pos.m_x, sizeof(double));
				    lcFlux.write((char*)&g.m_site[q].m_pos.m_y, sizeof(double));
				    lcFlux.write((char*)&g.m_site[q].m_pos.m_z, sizeof(double));
				    lcFlux.write((char*)&localCurrent, sizeof(double));	
				    double periodic = (double) (int) g.m_site[s].m_periodic[l];
				    lcFlux.write((char*)&periodic, sizeof(double));	
			    }

			//}
		    }

		    lcFlux.close();
		}
		else
		{
		    std::cerr << "Can't open file to write local current!" << '\n';
		}

	    }
	}
    }
}


doubleType checkTransportConvergence(Graph &g)
{
  doubleType maxDeviation(-1e15);
  
  for (int v(0) ; v < g.m_site.size() / g.m_nbrLevels / repetition ; v++)
  {
      for(int u(0) ; u < g.m_nbrLevels ; u++)
      {
	  int i(u*g.m_site.size()/g.m_nbrLevels + g.m_orderMeanField[v]);
	  doubleType inNode(0);
	  doubleType outNode(0);  

	  for (int j(0) ; j < g.m_site[i].m_neighbor.size() ; j ++)
	  {
		  int k(g.m_site[i].m_neighborIdx[j]);
		  inNode += g.m_site[k].m_occup * (1.0 - g.m_site[i].m_occup) * g.m_site[i].m_hopFrom[j];
		  outNode += g.m_site[i].m_occup * (1.0 - g.m_site[k].m_occup) * g.m_site[i].m_hopTo[j];

	  }
	  
	  //doubleType diff( (inNode-outNode) / inNode);
	  doubleType diff( (inNode-outNode) );
	  if(fabs(diff) > maxDeviation) { maxDeviation = fabs(diff);}
      }
  }
  
  return maxDeviation;
  
}



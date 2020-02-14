#include "System.h"

void System::readSystem(std::string fileNameTop, std::string fileNameCoord)
{
    std::cout << "Reading system coordinates... ";

    /// Read system ///
    std::string line;
    double nbr;
    Atom a;
    int fragIdOld(1);
    int chainIdOld(1);
    Polymer poly;
    RigidFragment mono;

    std::ifstream inFlux(fileNameCoord.c_str());
    if (inFlux.is_open())
    {
        int counter(0);
        inFlux >> line;
        inFlux >> nbr;
        m_boxSize.m_x = nbr;
        inFlux >> nbr;
        m_boxSize.m_y = nbr;
        inFlux >> nbr;
        m_boxSize.m_z = nbr;
        getline(inFlux,line);
        m_headerLine = line;

        // Read all atoms and add rigid fragments and polymers to the system
        bool firstAtom(true);
        while (!inFlux.eof())
        {

            getline(inFlux, line);
            if(line.size() >= 10)
            {
                std::string element = line.substr(11,5);
                std::string resName = line.substr(18,3);
                int siteId = std::stoi(line.substr(22,4));
                double x = std::stod(line.substr(26,12));
                double y = std::stod(line.substr(38,8));
                double z = std::stod(line.substr(46,8));
                int fragId = std::stoi(line.substr(54,6));
                int chainId = std::stoi(line.substr(60,6));

                a.m_name = 0;
                a.m_pos.m_x = x;
                a.m_pos.m_y = y;
                a.m_pos.m_z = z;
                element.erase(std::remove(element.begin(), element.end(),' '),element.end());
                if(element == "C") {a.m_type = C;}
                if(element == "O") {a.m_type = O;}
                if(element == "H") {a.m_type = H;}
                if(element == "N") {a.m_type = N;}
                if(element == "S") {a.m_type = S;}

                if(firstAtom)
                {
                    mono.m_fragNbr = fragId-1;
                    mono.m_molNbr = chainId-1;
                    mono.m_name = resName;
                    mono.m_siteNbr = siteId - 1;
                    mono.m_type = electron;
                }

                if(chainId != chainIdOld)
                {
                    poly.m_fragment.push_back(mono);
                    //std::cout << "chg chain " << mono.m_fragNbr << " " << mono.m_molNbr << " " << mono.m_siteNbr << std::endl;
                    m_polymer.push_back(poly);
                    poly.m_fragment.clear();
                    mono.m_atom.clear();
                    mono.m_fragNbr = fragId-1;
                    mono.m_molNbr = chainId-1;
                    mono.m_name = resName;
                    mono.m_siteNbr = siteId - 1;
                    mono.m_type = electron;
                }
                else if(fragId != fragIdOld)
                {
                    poly.m_fragment.push_back(mono);
                    //std::cout << "chg frag " << mono.m_fragNbr << " " << mono.m_molNbr << " " << mono.m_siteNbr << std::endl;
                    mono.m_atom.clear();
                    mono.m_fragNbr = fragId-1;
                    mono.m_molNbr = chainId-1;
                    mono.m_name = resName;
                    mono.m_siteNbr = siteId - 1;
                    mono.m_type = electron;
                }

                mono.m_atom.push_back(a);

                chainIdOld = chainId;
                fragIdOld = fragId;
                counter ++;
                firstAtom = false;

            }
        }

        // When we are done reading the file we need to add the last rigid fragment and polymer to the system
        //std::cout << "last " << mono.m_fragNbr << " " << mono.m_molNbr << " " << mono.m_siteNbr << std::endl;
        poly.m_fragment.push_back(mono);
        m_polymer.push_back(poly);
        m_nbAtom = counter;
        inFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "Done" << '\n';
}


void System::writeSystem(std::string fileNameCoord)
{
    std::cout << "Writing the system... ";
    int counter(1);
    Atom a;

    std::ofstream outFlux(fileNameCoord.c_str());
    if (outFlux.is_open())
    {
        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_boxSize.m_x;
        outFlux << std::setw(9) << m_boxSize.m_y;
        outFlux << std::setw(9) << m_boxSize.m_z;
        outFlux << m_headerLine;
        outFlux << '\n';

        for (int i(0) ; i < m_polymer.size() ; i ++ )
        {
            for (int j(0) ; j < m_polymer[i].m_fragment.size() ; j ++ )
            {
                for(int k(0) ; k < m_polymer[i].m_fragment[j].m_atom.size() ; k++)
                {
                    a = m_polymer[i].m_fragment[j].m_atom[k];

                    outFlux << "ATOM";
                    outFlux << std::setw(7);
                    outFlux << counter;
                    outFlux << std::setw(3);
                    std::string line;
                    if(a.m_type == C) {line = "C";}
                    if(a.m_type == O) {line = "O";}
                    if(a.m_type == H) {line = "H";}
                    if(a.m_type == S) {line = "S";}
                    if(a.m_type == N) {line = "N";}
                    outFlux << line;
                    outFlux << std::setw(8);
                    outFlux << m_polymer[i].m_fragment[j].m_name << " 0";
                    outFlux << std::setw(4);
                    outFlux << m_polymer[i].m_fragment[j].m_siteNbr;
                    outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                    outFlux << a.m_pos.m_x;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_y;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_z;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                    //outFlux << (double) id1;
                    outFlux << m_polymer[i].m_fragment[j].m_fragNbr ;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                    outFlux << (double) m_polymer[i].m_fragment[j].m_type;
                    outFlux << '\n';

                    counter ++;
                }
            }
        }

        outFlux << "TER" << '\n' << "ENDMDL" << '\n' ;
        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void System::writeSystemGro(std::string fileNameGro)
{
    std::cout << "Writing the system (Gromacs format)... ";
    int counter(1);
    Atom a;

    if(fileNameGro == "")
    {
        fileNameGro = "System.gro";
    }

    std::ofstream outFlux(fileNameGro.c_str());
    if (outFlux.is_open())
    {

        outFlux << "SYSTEM";
        outFlux << '\n';
        outFlux << std::setw(5) << m_nbAtom;
        outFlux << '\n';

        for (int i(0) ; i < m_polymer.size() ; i ++ )
        {
            for (int j(0) ; j < m_polymer[i].m_fragment.size() ; j ++ )
            {
                for(int k(0) ; k < m_polymer[i].m_fragment[j].m_atom.size() ; k++)
                {
                    a = m_polymer[i].m_fragment[j].m_atom[k];

                    outFlux << std::setw(5);
                    outFlux << m_polymer[i].m_fragment[j].m_siteNbr+1;
                    outFlux << std::setw(3);
                    outFlux << m_polymer[i].m_fragment[j].m_name;
                    outFlux << std::setw(7);

                    std::string line;
                    if(a.m_type == C) {line = "C";}
                    if(a.m_type == O) {line = "O";}
                    if(a.m_type == H) {line = "H";}
                    if(a.m_type == S) {line = "S";}
                    if(a.m_type == N) {line = "N";}

                    outFlux << line;
                    outFlux << std::setw(5);
                    outFlux << counter%100000;

                    outFlux << std::fixed << std::setprecision(3) << std::setw(8) ;
                    outFlux << a.m_pos.m_x * 0.1;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_y * 0.1;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_z * 0.1;
                    /*outFlux << std::fixed << std::setprecision(4) << std::setw(8);
                    outFlux << 0.0;
                    outFlux << std::setw(8);
                    outFlux << 0.0;
                    outFlux << std::setw(8);
                    outFlux << 0.0;*/

                    outFlux << '\n';

                    counter ++;
                }

            }

        }

        outFlux << std::fixed << std::setprecision(4) << m_boxSize.m_x * 0.1 << " " << m_boxSize.m_y * 0.1 << " " << m_boxSize.m_z * 0.1 << '\n';

        outFlux.close();
    }
    else
    {
        std::cerr << "Unable to open coordinates file" << std::endl;
    }
    std::cout << "done" << '\n';
}


void System::initGeom(double vectLength)
{
    std::cout << "Initialize rigid segments properties... " << std::endl;
    
    std::string orientationFile("orientation.tcl");
    std::ifstream f(orientationFile.c_str());
    bool isOrient; 
    isOrient = f.good();
    
    if(isOrient)
    {
      readFragmentVector(orientationFile, vectLength);

      // Be sure that we are good with PBC
      for ( int i(0) ; i < m_polymer.size() / 8 ; i++ )
      {
	  for ( int j(0) ; j < m_polymer[i].m_fragment.size() ; j ++ )
	  {
	      for(int k(1) ; k < 8 ; k++)
	      {
		int idxRep((m_polymer.size()/8)*k+i);
		// Translation vector of replication
		Vect3 translate;
		if(k==4 || k==5 || k==6 || k==7) {translate.m_x = m_boxSize.m_x/2.0;}
		if(k==2 || k==3 || k==6 || k==7) {translate.m_y = m_boxSize.m_y/2.0;}
		if(k==1 || k==3 || k==5 || k==7) {translate.m_z = m_boxSize.m_z/2.0;}

		m_polymer[idxRep].m_fragment[j].m_center = m_polymer[i].m_fragment[j].m_center + translate;
		m_polymer[idxRep].m_fragment[j].m_orientN = m_polymer[i].m_fragment[j].m_orientN;
		m_polymer[idxRep].m_fragment[j].m_orientP = m_polymer[i].m_fragment[j].m_orientP;
		m_polymer[idxRep].m_fragment[j].m_orientT = m_polymer[i].m_fragment[j].m_orientT;
	      }
	  }
      }
      
    }
    else
    {
      // Modified to take into account 8 times replication by pbc
#pragma omp parallel for
      for ( int i(0) ; i < m_polymer.size() / 8 ; i++ )
      {
	  for ( int j(0) ; j < m_polymer[i].m_fragment.size() ; j ++ )
	  {
	      bool invert(false);
	      if( j % 2 == 0 ) {invert = true;}

	      m_polymer[i].m_fragment[j].computeGeom(invert);
	      
	      for(int k(1) ; k < 8 ; k++)
	      {
		int idxRep((m_polymer.size()/8)*k+i);
		// Translation vector of replication
		Vect3 translate;
		if(k==4 || k==5 || k==6 || k==7) {translate.m_x = m_boxSize.m_x/2.0;}
		if(k==2 || k==3 || k==6 || k==7) {translate.m_y = m_boxSize.m_y/2.0;}
		if(k==1 || k==3 || k==5 || k==7) {translate.m_z = m_boxSize.m_z/2.0;}

		m_polymer[idxRep].m_fragment[j].m_center = m_polymer[i].m_fragment[j].m_center + translate;
		m_polymer[idxRep].m_fragment[j].m_orientN = m_polymer[i].m_fragment[j].m_orientN;
		m_polymer[idxRep].m_fragment[j].m_orientP = m_polymer[i].m_fragment[j].m_orientP;
		m_polymer[idxRep].m_fragment[j].m_orientT = m_polymer[i].m_fragment[j].m_orientT;
	      }
	  }
      }
    }

    std::cout << "Done" << '\n';
}


int newBoxPbc(int box, const Vect3 &pbc)
{
    int newBoxId(box);

    if(pbc.m_x != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,x);
    }
    if(pbc.m_y != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,y);
    }
    if(pbc.m_z != 0)
    {
        newBoxId = lookupTabPbc(newBoxId,z);
    }

    return newBoxId;
}

BoxIdx newBoxNoPbc(int box1, int box2)
{
    bool revert(box2 < box1);
    int b1(box1);
    int b2(box2);
    if(revert)
    {
      int bTemp(b1);
      b1 = b2;
      b2 = bTemp;
    }
    
    int newb1(b1);
    int newb2(b2);
    
    if(b1 == 1 and b2 == 3) {newb1 = 0; newb2 = 2;}
    if(b1 == 1 and b2 == 5) {newb1 = 0; newb2 = 4;}
    if(b1 == 1 and b2 == 7) {newb1 = 0; newb2 = 6;}
    if(b1 == 2 and b2 == 3) {newb1 = 0; newb2 = 1;}
    if(b1 == 2 and b2 == 6) {newb1 = 0; newb2 = 4;}
    if(b1 == 2 and b2 == 7) {newb1 = 0; newb2 = 5;}
    if(b1 == 3 and b2 == 5) {newb1 = 2; newb2 = 4;}
    if(b1 == 3 and b2 == 6) {newb1 = 1; newb2 = 4;}
    if(b1 == 3 and b2 == 7) {newb1 = 0; newb2 = 4;}
    if(b1 == 4 and b2 == 5) {newb1 = 0; newb2 = 1;}
    if(b1 == 4 and b2 == 6) {newb1 = 0; newb2 = 2;}
    if(b1 == 4 and b2 == 7) {newb1 = 0; newb2 = 3;}
    if(b1 == 5 and b2 == 6) {newb1 = 1; newb2 = 2;}
    if(b1 == 5 and b2 == 7) {newb1 = 0; newb2 = 2;}
    if(b1 == 6 and b2 == 7) {newb1 = 0; newb2 = 1;}
    
    if(revert)
    {
      int bTemp(newb1);
      newb1 = newb2;
      newb2 = bTemp;
    }    
    
    BoxIdx boxIdx;
    boxIdx.m_revert = revert;
    boxIdx.m_b1 = newb1;
    boxIdx.m_b2 = newb2;
    
    return boxIdx;
}


int lookupTabPbc(int box, Direction dir)
{
    if(box==0)
    {
        if(dir == x) {return 4;}
        if(dir == y) {return 2;}
        if(dir == z) {return 1;}
    }
    if(box==1)
    {
        if(dir == x) {return 5;}
        if(dir == y) {return 3;}
        if(dir == z) {return 0;}
    }
    if(box==2)
    {
        if(dir == x) {return 6;}
        if(dir == y) {return 0;}
        if(dir == z) {return 3;}
    }
    if(box==3)
    {
        if(dir == x) {return 7;}
        if(dir == y) {return 1;}
        if(dir == z) {return 2;}
    }
    if(box==4)
    {
        if(dir == x) {return 0;}
        if(dir == y) {return 6;}
        if(dir == z) {return 5;}
    }
    if(box==5)
    {
        if(dir == x) {return 1;}
        if(dir == y) {return 7;}
        if(dir == z) {return 4;}
    }
    if(box==6)
    {
        if(dir == x) {return 2;}
        if(dir == y) {return 4;}
        if(dir == z) {return 7;}
    }
    if(box==7)
    {
        if(dir == x) {return 3;}
        if(dir == y) {return 5;}
        if(dir == z) {return 6;}
    }
    
    return -1;
}


void System::readFragmentVector(std::string fileNameFragmentVector, double vectLength)
{
    std::ifstream of(fileNameFragmentVector.c_str());
    std::string dummyString;
    double dummyDouble;
    
    if(of)
    {
        for(int i(0) ; i < m_polymer.size() ; i++)
        {
            for (int j(0) ; j < m_polymer[i].m_fragment.size() ; j++ )
            {
	        Vect3 v1, v2, v3, v4;

                // Read orient N
		of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v1.m_x = dummyDouble;
		of >> dummyDouble; 
		v1.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v1.m_z = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v2.m_x = dummyDouble;
		of >> dummyDouble; 
		v2.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v2.m_z = dummyDouble;	
		for(int k(0) ; k < 23 ; k++)
		{
		  of >> dummyString;
		}
		
                // Read orient P
		of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v1.m_x = dummyDouble;
		of >> dummyDouble; 
		v1.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v1.m_z = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v3.m_x = dummyDouble;
		of >> dummyDouble; 
		v3.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v3.m_z = dummyDouble;	
		for(int k(0) ; k < 23 ; k++)
		{
		  of >> dummyString;
		}
		
                // Read orient T
		of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v1.m_x = dummyDouble;
		of >> dummyDouble; 
		v1.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v1.m_z = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(1);
		dummyDouble = std::stod(dummyString);
		v4.m_x = dummyDouble;
		of >> dummyDouble; 
		v4.m_y = dummyDouble;
		of >> dummyString;
		dummyString = dummyString.substr(0,dummyString.size()-1);
		dummyDouble = std::stod(dummyString);
		v4.m_z = dummyDouble;	
		for(int k(0) ; k < 23 ; k++)
		{
		  of >> dummyString;
		}
		
		v2 = v2 - v1; v2 = v2 / vectLength; 
		v3 = v3 - v1; v3 = v3 / vectLength; 
		v4 = v4 - v1; v4 = v4 / vectLength; 
                
                m_polymer[i].m_fragment[j].m_center = v1;
                m_polymer[i].m_fragment[j].m_orientN = v2;
                m_polymer[i].m_fragment[j].m_orientP = v3;
                m_polymer[i].m_fragment[j].m_orientT = v4;

            }
        }

        for(int i(0) ; i < m_molecule.size() ; i++)
        {
	    Vect3 v1, v2, v3, v4;

	    // Read orient N
	    of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v1.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_z = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v2.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v2.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v2.m_z = dummyDouble;	
	    for(int k(0) ; k < 23 ; k++)
	    {
	      of >> dummyString;
	    }
	    
	    // Read orient P
	    of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v1.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_z = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v3.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v3.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v3.m_z = dummyDouble;	
	    for(int k(0) ; k < 23 ; k++)
	    {
	      of >> dummyString;
	    }
	    
	    // Read orient T
	    /*of >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString >> dummyString;		
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v1.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v1.m_z = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(1);
	    dummyDouble = std::stod(dummyString);
	    v4.m_x = dummyDouble;
	    of >> dummyDouble; 
	    v4.m_y = dummyDouble;
	    of >> dummyString;
	    dummyString = dummyString.substr(0,dummyString.size()-1);
	    dummyDouble = std::stod(dummyString);
	    v4.m_z = dummyDouble;	
	    for(int k(0) ; k < 23 ; k++)
	    {
	      of >> dummyString;
	    }
	    */
	    
	    v2 = v2 - v1; v2 = v2 / vectLength; 
	    v3 = v3 - v1; v3 = v3 / vectLength; 
	    //v4 = v4 - v1; v4 = v4 / vectLength; 
	    
	    m_molecule[i].m_center = v1;
	    m_molecule[i].m_orientN = v2;
	    m_molecule[i].m_orientP = v3;
	    //m_molecule[i].m_orientT = v4;
	    
        }

        of.close();

    }
    else
    {
        std::cerr << "Can't open script file for reading rigid fragment orientations" << '\n';
    }
}

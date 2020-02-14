#include "Parse.h"

using namespace std;

void parseBigPdb(string fileI, string fileO, string fileRef)
{
    cout << "Parsing the big pdb file... ";
    string header;
    string line;

    // Read the reference file
    vector<int> atomIdVect;
    vector<int> siteIdVect;
    vector<int> fragIdVect;
    vector<int> chainIdVect;
    int atomId(0);
    ifstream refF(fileRef);
    getline(refF,line);
    while(!refF.eof())
    {
        getline(refF, line);
        if(line.size() >= 10)
        {
            int siteId = std::stoi(line.substr(22,4));
            int fragId = std::stoi(line.substr(54,6));
            int chainId = std::stoi(line.substr(60,6));
            atomIdVect.push_back(atomId+1);
            siteIdVect.push_back(siteId);
            fragIdVect.push_back(fragId);
            chainIdVect.push_back(chainId);
            atomId ++;

        }
    }
    refF.close();
    // Extend the reference information with the new box
    int nbSiteId(0);
    for(int i(0) ; i < siteIdVect.size() ; i++) { if( siteIdVect[i] > nbSiteId ) {nbSiteId = siteIdVect[i];} }
    //int nbSiteId(siteIdVect[atomId-1]);
    int nbChainId(chainIdVect[atomId-1]);
    for(int i(1) ; i < 8 ; i++)
    {
        for(int j(0) ; j < atomId ; j++)
        {
	    int siteIdReal(-1);
	    if(siteIdVect[j] != -1)
	    {
	      siteIdReal = siteIdVect[j]+nbSiteId*i;
	    }
            int fragIdReal(fragIdVect[j]);
            int chainIdReal(chainIdVect[j]+nbChainId*i);
            int atomIdReal(atomIdVect[j]+atomId*i);
            atomIdVect.push_back(atomIdReal);
            siteIdVect.push_back(siteIdReal);
            fragIdVect.push_back(fragIdReal);
            chainIdVect.push_back(chainIdReal);
        }

    }


    // Read input file from vmd extended box
    ifstream is(fileI);
    ofstream outFlux(fileO);
    // Read header
    getline(is, line);
    header = line;
    // Read file
    vector<string> data;
    while(!is.eof())
    {
        getline(is, line);
        if(line.size() >= 10)
        {
            data.push_back(line);
        }
    }
    is.close();

    /// Write in new file with correct information
    outFlux << header << '\n';
    int counter(1);

    /// Write polymer
    for (int i(0) ; i < data.size() ; i ++)
    {
        line = data[i];
        std::string element = line.substr(11,5);
        std::string resName = line.substr(17,3);
        int siteId = std::stoi(line.substr(22,4));
        double x = std::stod(line.substr(26,12));
        double y = std::stod(line.substr(38,8));
        double z = std::stod(line.substr(46,8));
        double fragId = std::stod(line.substr(54,6));
        double chainId = std::stod(line.substr(60,6));
        element.erase(std::remove(element.begin(), element.end(),' '),element.end());

        outFlux << "ATOM";
        outFlux << std::setw(7);
        outFlux << counter % 100000;
        outFlux << std::setw(5);
        outFlux << element;
        outFlux << std::setw(5);
        outFlux << resName;
        outFlux << std::setw(1);
        outFlux << "0";
        outFlux << std::setw(4);
        outFlux << siteIdVect[i] % 10000;
        outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
        outFlux << x;
        outFlux << std::setw(8);
        outFlux << y;
        outFlux << std::setw(8);
        outFlux << z;
        outFlux << std::fixed << std::setprecision(2) << std::setw(6);
        outFlux << fragIdVect[i] ;
        outFlux << std::fixed << std::setprecision(2) << std::setw(6);
        outFlux << chainIdVect[i] ;
        outFlux << std::setw(10);
        outFlux << element;
        outFlux << '\n';

        counter ++;
    }

    outFlux << "END" << '\n';

    outFlux.close();
    cout << "Done." << '\n';
}


void extendPdb(string fileI)
{
    cout << "Extending the pbd file with periodic boundary conditions... ";

    /// Create and execute the Tcl script that will perform the extension (Note VMD has to be in the PATH) ///
    ofstream script("extendPdb.tcl");
    string fileIShort = fileI.substr(0, fileI.size()-4);
    string filePdbT(fileIShort + "_bigT.pdb");
    string filePdb(fileIShort + "_big.pdb");
    if (script.is_open())
    {
        script << "# load required package" << '\n' << "package require topotools" << '\n' << "# load a molecule"
        << '\n' << "set mol [mol new " << fileI << " type pdb waitfor all]" << '\n' << "# do the magic"
        << '\n' << "set newmol [::TopoTools::replicatemol $mol 2 2 2 ]" << '\n' << "animate write pdb " << filePdbT
        << " $newmol" << '\n' << "quit" ;

        script.close();

        system("vmd -dispdev text -e extendPdb.tcl");
    }
    else
    {
        cout << "Can't open tcl script file to extend pdb volume." << '\n';
    }


    /// Open the new (extended) pdb file and write the header ///
    string line;
    ifstream pdbT(filePdbT);
    ofstream pdb(filePdb);

    if (pdb.is_open() && pdbT.is_open())
    {
        getline(pdbT, line);
        pdb << line << '\n';

        while(!pdbT.eof())
        {
            getline(pdbT, line);
            pdb << line << '\n';

        }

        pdb.close();
        pdbT.close();

    }
    else
    {
        cout << "Can't open extended pdb file to write header." << '\n';
    }

    cout << "Done." << '\n';
}


void extendTop(string fileI)
{
    cout << "Writing extended top file... ";

    string topExtended = fileI.substr(0, fileI.size()-4) + "_big.top";

    ifstream top(fileI);
    ofstream topEx(topExtended);
    string line("");

    if (top.is_open() && topEx.is_open())
    {
        while(line != "[ molecules ]" )
        {
            getline(top, line);
            topEx << line << '\n';
        }
        getline(top, line);
        topEx << line << '\n';
        int nbPdt(0);
        int nbTos(0);
        top >> line >> nbPdt >> line >> nbTos;
        topEx << "pedot    " << 8*nbPdt << '\n' << "tos    " << 8*nbTos << '\n';


        top.close();
        topEx.close();

    }
    else
    {
        cout << "Can't open top file to write." << '\n';
    }

    cout << "Done." << '\n';

}

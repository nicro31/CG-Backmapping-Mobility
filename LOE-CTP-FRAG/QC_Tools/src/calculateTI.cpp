#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <locale>
#include <functional>
#include <algorithm>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <sys/stat.h>
#include <math.h>

#include "MATRIX/matrix.hpp"
#include "IO/io.hpp"
#include "IO/argumentparser.hpp"
#include "QC_FUNCTIONS/qc_functions.hpp"
#include "PARAMETERS/parameters.hpp"
#include "LOG/log.hpp"
#include "IO/FILE_READERS/punreader.hpp"
#include "IO/FILE_READERS/logreader.hpp"

#include "calcJconfig.hpp"
#include "calculateTI.hpp"

#define hartreeToeV 27.21138602

using namespace catnip;
using namespace std;

vector<vector<double> > calculateTI(int argc,const char *argv[], int orbLevel){
  
  /*cout << "Running calc_J VERSION " << calcJ_VERSION_MAJOR << ".";
  cout << calcJ_VERSION_MINOR << endl;
  
  cout <<argc << endl;
  for(int i = 0; i < argc; ++i)
  {
        cout << argv[i] << '\n';
  }*/
  
  vector<vector<double> > resultCollec;

  string line;
  LOG("Preparing parser",1);
  auto ArgPars = prepareParser();
  LOG("Parsing arguments",1);
  ArgPars->parse(argv,argc);
  LOG("Preparing parameter object",1);
  auto par = prepareParameters(ArgPars);
  

  /*cout << "log file for first monomer is:      "+par->getLog1()+'\n';
  cout << "log file for second monomer is:     "+par->getLog2()+'\n';
	cout << "log file for dimer is:              "+par->getLogP()+'\n';
	cout << "pun file for the first monomer is:  "+par->getPun1()+'\n';
	cout << "pun file for the second monomer is: "+par->getPun2()+'\n';
	cout << "pun file for the dimer is:          "+par->getPunP()+'\n';
  */
 
	
	//Open the .pun file find the total number of molecular orbitals
 
  LOG("Reading pun files",1);
  LOG("Reading pun file: "+par->getPunP(),2);
  PunReader pr_P(par->getPunP());
  pr_P.read();
  
  LOG("Reading pun file: "+par->getPun1(),2);
  PunReader pr_1(par->getPun1());
  pr_1.read();
  
  LOG("Reading pun file: "+par->getPun2(),2);
  PunReader pr_2(par->getPun2());
  pr_2.read();

  LOG("Reading log files",1);
  LOG("Reading log file: "+par->getLogP(),2);
  LogReader lr_P(par->getLogP());
  lr_P.read();
  LOG("Reading log file: "+par->getLog1(),2);
  LogReader lr_1(par->getLog1());
  lr_1.read();
  LOG("Reading log file: "+par->getLog2(),2);
  LogReader lr_2(par->getLog2());
  lr_2.read();
  // Load in general coefficients 

  // If there are only alpha orbitals then we can assume that restricted HF 
  // was used

  // No need to worry about beta orbitals


  {
    Matrix * mat_S = lr_P.getOverlapMatrix();

    Matrix * mat_P_Coef = pr_P.getCoefsMatrix(par->getSpinP());
    auto vec_P_OE = lr_P.getOE(par->getSpinP());
    Matrix * mat_P_OE = new Matrix(vec_P_OE);

    int HOMO1 = lr_1.getHOMOLevel(par->getSpin1());
    LOG("Getting "+par->getSpin1()+" of monomer 1",2);
    Matrix * mat_1_Coef = pr_1.getCoefsMatrix(par->getSpin1());
    auto vec_1_OE = lr_1.getOE(par->getSpin1());
    Matrix * mat_1_OE = new Matrix(vec_1_OE);

    int HOMO2 = lr_2.getHOMOLevel(par->getSpin2());
    LOG("Getting "+par->getSpin2()+" of monomer 2",2);
    Matrix * mat_2_Coef = pr_2.getCoefsMatrix(par->getSpin2());
    auto vec_2_OE = lr_2.getOE(par->getSpin2());
    Matrix * mat_2_OE = new Matrix(vec_2_OE);
    
    // Unscramble dimer coef and energies first need to see how the dimer
    // and monomer coefficients line up. To determine how the ceofficients
    // line up we will first look at how the atoms appear in each of the 
    // .gjf files. We will also check to see how many coefficients are 
    // assocaited with each of the atoms by checking the .log files. Given
    // the position of the atoms in the monomer unit and the positions of
    // the atoms in the dimer we can determine how the coefficients need 
    // to be rearranged.     
    auto coord_P = lr_P.getCoords();
    auto coord_1 = lr_1.getCoords();
    auto coord_2 = lr_2.getCoords();

    // Convert coords to matrices
    Matrix coord_P_mat(coord_P);
    Matrix coord_1_mat(coord_1);
    Matrix coord_2_mat(coord_2);

    auto basis_P = lr_P.getBasisFuncCount();  
    auto basis_1 = lr_1.getBasisFuncCount();  
    auto basis_2 = lr_2.getBasisFuncCount();  

    int MO1 = mat_1_OE->get_rows();
    int MO2 = mat_2_OE->get_rows();

    pair<int,int> Orbs1 = { MO1, HOMO1 };
    pair<int,int> Orbs2 = { MO2, HOMO2 };

    LOG("Creating transfercomplex",1);
    TransferComplex TC(
        mat_1_Coef,
        mat_2_Coef,
        mat_P_Coef,
        Orbs1,
        Orbs2,
        mat_S,
        mat_P_OE,
        par->getCounterPoise());

    // Set the transfer complex to counterpoise if it is the case. 

    // If the basis function search returns 0 for any of the components then
    // we cannot automatically determine what the transfer integral is
    if(basis_1.size()!=0 && basis_2.size()!=0 && basis_P.size()!=0){
      LOG("Unscrambling matrices",1);
      TC.unscramble(
          coord_1_mat,
          coord_2_mat,
          coord_P_mat,
          basis_P,
          basis_2);
    }
    
    map<string,string> orbitaltypes;
    map<string,int> orbitalnums;

    /////// N. ROLLAND, modification to calculate transfer integrals for several orbitals at the same time /////////
    // Orbital type and a map of the corresponding number
    // E.g.
    // orbital type    orbital number
    // "mon1" "LUMO"    "mon1" +3
    // "mon2" "HOMO"    "mon2"  0
    //
    //    monomer1 LUMO+3
    //    monomer2 HOMO
    
    
    orbitaltypes["mon1"]=par->getOrbType1();
    orbitaltypes["mon2"]=par->getOrbType2();
    
    int orbInc(1);
    if(orbitaltypes["mon1"] == "HOMO"){orbInc = -1;} 
    
    for(int i(0) ; abs(i) <= orbLevel ; i += orbInc)
    {
      for(int j(0) ; abs(j) <= orbLevel ; j += orbInc)
      {

	par->setOrbNum1(i);
	par->setOrbNum2(j);

    /*cout << endl;
    cout << "Dimer     Spin " << par->getSpinP() << endl;

    cout << "Monomer 1 Spin " << par->getSpin1() << " ";
    if(par->getOrbNum1()==0){
      cout << "Orbital " << par->getOrbType1() << endl;
    }else if(par->getOrbNum1()>0){
      cout << "Orbital " << par->getOrbType1();
      cout << "+" << par->getOrbNum1() << endl;
    }else{
      cout << "Orbital " << par->getOrbType1();
      cout << par->getOrbNum1() << endl;
    }

    cout << "Monomer 2 Spin " << par->getSpin2() << " ";
    if(par->getOrbNum2()==0){
      cout << "Orbital " << par->getOrbType2() << endl;
    }else if(par->getOrbNum2()>0){
      cout << "Orbital " << par->getOrbType2();
      cout << "+" << par->getOrbNum2() << endl;
    }else{
      cout << "Orbital " << par->getOrbType2();
      cout << par->getOrbNum2() << endl;
    }*/
    
	  /// Get energies ///
	  //cout << "Getting energies of monomer 1 orbital " << orbitaltypes["mon1"] << i << " and monomer 2 orbital " << orbitaltypes["mon2"] << j << "... ";
	  double E1(mat_1_OE->get_elem( HOMO1 + (orbitaltypes["mon1"] == "LUMO") + i ) );
	  double E2(mat_2_OE->get_elem( HOMO2 + (orbitaltypes["mon2"] == "LUMO") + j ) );
	  //cout << "Done!" << endl;

	  /// Get transfer integral ///
	  //cout << "Calculating transfer integral from orbital " << orbitaltypes["mon1"] << i << " to " << orbitaltypes["mon2"] << j << "... " << std::endl; 
	  orbitalnums["mon1"]=i;
	  orbitalnums["mon2"]=j;
	  vector<double> result = TC.calcJ(orbitaltypes,orbitalnums);
	  result[0] = E1*hartreeToeV;
	  result[1] = E2*hartreeToeV;
	  //cout << "E_mon1 = " << E1*hartreeToeV << "eV" << endl << "E_mon2 = " << E2*hartreeToeV << "eV" << endl;
	  resultCollec.push_back(result);
	  //cout << "Done!" << std::endl;
      }      
    }
    
    //orbitalnums["mon1"]=par->getOrbNum1();
    //orbitalnums["mon2"]=par->getOrbNum2();
    //LOG("Calculating transfer integral",1);
    //TC.calcJ(orbitaltypes,orbitalnums);
    
    
    delete mat_P_OE;
    delete mat_1_OE;
    delete mat_2_OE;

  }

	return resultCollec;
}

#ifndef SYSTEM_H_INCLUDED
#define SYSTEM_H_INCLUDED
#include "Units.h"

struct BoxIdx
{
  bool m_revert;
  int m_b1;
  int m_b2;
};

struct System
{
    std::string m_headerLine;
    Vect3 m_boxSize;
    int m_nbAtom;
    std::vector<Polymer> m_polymer;
    std::vector<RigidFragment> m_molecule;
    Topo m_topo;

    void readSystem(std::string fileNameTop, std::string fileNameCoord);
    void initGeom(double vectLength); // Compute rigid fragment properties (center and orientation)
    void writeSystem(std::string fileNameCoord);
    void writeSystemGro(std::string fileNameGro = "");
    void readFragmentVector(std::string fileNameFragmentVector, double vectLength);


};


int newBoxPbc(int box, const Vect3 &pbc);

BoxIdx newBoxNoPbc(int box1, int box2);

int lookupTabPbc(int box, Direction dir);


#endif // SYSTEM_H_INCLUDED

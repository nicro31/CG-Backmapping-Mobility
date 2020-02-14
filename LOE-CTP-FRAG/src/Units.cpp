#include "Units.h"

AtomType mass2Type(float mass)
{
    if( mass == (float) 12.0100 ) { return C; }
    if( mass == (float) 16.0000 ) { return O; }
    if( mass == (float) 32.0600 ) { return S; }
    if( mass == (float) 1.0080 ) { return H; }
    if( mass == (float) 14.007 ) { return N; }

    return Z;
}

void RigidFragment::computeGeom(bool invert)
{
        Vect3 center;
        Vect3 orientN;
        Vect3 orientP;
        Vect3 orientT;
	
	std::vector<Vect3> inPoints;

	int nbAtom(0);
	for(int i(0) ; i < m_atom.size() ; i++)
	{
	  //if(m_atom[i].m_type != H)
	  if(true)
	  {
	    nbAtom ++;
	    Vect3 inV(m_atom[i].m_pos.m_x, m_atom[i].m_pos.m_y, m_atom[i].m_pos.m_z);
	    inPoints.push_back(inV);
	  }
	}
	
	orientation(inPoints, center, orientN, orientT, orientP);
	
        orientP = orientP / orientP.norm();
        orientN = orientN / orientN.norm();
        orientT = orientT / orientT.norm();

        m_center = center;
        m_orientN = orientN;
        m_orientP = orientP;
        m_orientT = orientT;

}


bool crossingRaySphere(const Vect3 &rf1,const Vect3 &rf2,const Vect3 &rfInfluence,const double &sphereR, const Vect3 &boxSize)
{
    Vect3 vrf1(0,0,0);
    Vect3 vrf2 = Vect3::displacementPbc(rf1, rf2, boxSize);
    Vect3 vrfInfluence = Vect3::displacementPbc(rf1, rfInfluence, boxSize);

    // We have to check if the rfInfluence is in between the segment limit or not
    bool isInside(true);
    Vect3 v = Vect3::displacementPbc(rf2, rfInfluence, boxSize);
    double angle1( acos( (vrf2/vrf2.norm()) * (vrfInfluence/vrfInfluence.norm()) ) );
    double angle2( acos( (vrf2/vrf2.norm()) * (v/v.norm()) ) );
    if( (angle1 > pi / 2.0) || (angle2 < pi / 2.0) ) {isInside = false;}

    double d(0);

    if(isInside)
    {
        Vect3 v1((vrfInfluence-vrf1)^(vrfInfluence-vrf2));
        Vect3 v2(vrf2-vrf1);
        d = (v1.norm() / v2.norm());
    }
    else
    {
        d = std::min(vrfInfluence.norm(),v.norm());
    }

    return(d<=sphereR);
}


bool operator==(Vect3 const& v1, Vect3 const& v2)
{
    if(v1.m_x == v2.m_x && v1.m_y == v2.m_y && v1.m_z == v2.m_z)
    {
        return true;
    }
    else
    {
        return false;
    }
}


void RigidFragment::orientation(const std::vector<Vect3> &inputPoints, Vect3 &center, Vect3 &orientN, Vect3 &orientT, Vect3 &orientP)
{
  
    int nbAtom(inputPoints.size());
    ApproxMVBB::Matrix3Dyn points(3, nbAtom);
	
    for(int i(0) ; i < nbAtom ; i++)
    {
	points.col(i) = ApproxMVBB::Vector3(inputPoints[i].m_x,inputPoints[i].m_y,inputPoints[i].m_z);
    }
      
      
    ApproxMVBB::OOBB oobb = ApproxMVBB::approximateMVBB(points,
							  0.001,
							  500,
							  5, /*increasing the grid size decreases speed */
							  0,
							  5);

    // To make all points inside the OOBB :
    ApproxMVBB::Matrix33 A_KI = oobb.m_q_KI.matrix().transpose();  // faster to store the transformation matrix first
    auto size                 = points.cols();
    for(unsigned int i = 0; i < size; ++i)
    {
	oobb.unite(A_KI * points.col(i));
    }

    // To make the box have a minimum extent of greater 0.1:
    // see also oobb.expandToMinExtentRelative(...)
    oobb.expandToMinExtentAbsolute(0.01);

    // Get oobb corners points
    bool coordinateSystemIsI(true);
    ApproxMVBB::Vector3List p;
    p = oobb.getCornerPoints();
      
    // Calculate oobb center
    double xCenter(0);
    double yCenter(0);
    double zCenter(0);
    for(int i(0) ; i < 8 ; i++)
    {
	xCenter += p[i][0];
	yCenter += p[i][1];
	zCenter += p[i][2];
    }
    xCenter /= 8; yCenter /= 8; zCenter /= 8;
    
    // Calculate oobb vectors
    std::vector<double> normal(3);
    std::vector<double> small(3);
    std::vector<double> large(3);
    std::vector<double> vx(3);
    std::vector<double> vy(3);
    std::vector<double> vz(3);
    
    vx[0] = p[1][0] - p[0][0];
    vx[1] = p[1][1] - p[0][1];
    vx[2] = p[1][2] - p[0][2];
    double normx(sqrt(vx[0]*vx[0] + vx[1]*vx[1] + vx[2]*vx[2]));
    
    vy[0] = p[2][0] - p[0][0];
    vy[1] = p[2][1] - p[0][1];
    vy[2] = p[2][2] - p[0][2];
    double normy(sqrt(vy[0]*vy[0] + vy[1]*vy[1] + vy[2]*vy[2]));
    
    vz[0] = p[4][0] - p[0][0];
    vz[1] = p[4][1] - p[0][1];
    vz[2] = p[4][2] - p[0][2];
    double normz(sqrt(vz[0]*vz[0] + vz[1]*vz[1] + vz[2]*vz[2]));
    
    double normMax( std::max( std::max(normx, normy), normz ) );
    double normMin( std::min( std::min(normx, normy), normz ) );
    
    if( normx == normMin )
    { 
	normal = vx;
	if( normy == normMax ) { large = vy; small = vz;}
	else { large = vz ; small = vy;}
    }
    if( normy == normMin )
    { 
	normal = vy;
	if( normx == normMax ) { large = vx; small = vz;}
	else { large = vz ; small = vx;}
    }
    if( normz == normMin )
    { 
	normal = vz;
	if( normy == normMax ) { large = vy; small = vx;}
	else { large = vx ; small = vy;}
    }

    center.m_x = xCenter; center.m_y = yCenter; center.m_z = zCenter;
    orientN.m_x = normal[0]; orientN.m_y = normal[1]; orientN.m_z = normal[2]; 
    orientT.m_x = large[0]; orientT.m_y = large[1]; orientT.m_z = large[2];
    orientP.m_x = small[0]; orientP.m_y = small[1]; orientP.m_z = small[2];
	  
}

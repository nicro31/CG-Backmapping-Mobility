#include "Graph.h"


void Graph::initSystem()
{
    std::string fileNameTop = m_fileNameTop;
    std::string fileNameCoord = m_fileNameCoord;
    std::string fileNameSave("partition.pdb");

    m_system.readSystem(fileNameTop, fileNameCoord);
    m_system.initGeom(m_vectLength);
    //m_system.writeSystem(fileNameSave);

}

void Graph::segmentSystem()
{
    std::cout << "Creating conjugated segments... ";

    /// Segmentation from site idx in the pdb file ///
    conjSegment site;
    Vect3 pos;
    int siteNbrOld(-1);
    for (int i(0) ; i < m_system.m_polymer.size() ; i ++)
    {
        for(int j(0) ; j < m_system.m_polymer[i].m_fragment.size(); j++)
        {
            int siteNbr(m_system.m_polymer[i].m_fragment[j].m_siteNbr);
	    //watch(siteNbr);
	    if(siteNbr != -2)
	    {
	      if( siteNbr != siteNbrOld && i+j != 0 )
	      {
		  pos = pos / site.m_fragment.size();
		  site.m_pos = pos;
		  m_site.push_back(site);
		  conjSegment newSite; site = newSite;
		  Vect3 newPos; pos = newPos;
	      }
	      site.m_fragment.push_back(&m_system.m_polymer[i].m_fragment[j]);
	      pos += m_system.m_polymer[i].m_fragment[j].m_center;	      
	      siteNbrOld = siteNbr;
	    }
        }
    }
    pos = pos / site.m_fragment.size();
    site.m_pos = pos;
    m_site.push_back(site);

    /// Init site geometry from rigid fragments geometry ///
    for (int i(0) ; i < m_site.size() ; i ++)
    {
        Vect3 orientN;
        Vect3 orientP;

        for (int j(0) ; j < m_site[i].m_fragment.size() ; j ++)
        {
            orientN += m_site[i].m_fragment[j]->m_orientN;
            orientP += m_site[i].m_fragment[j]->m_orientP;
        }

        orientN = orientN / (double) m_site[i].m_fragment.size(); orientN = orientN / orientN.norm();
        orientP = orientP / (double) m_site[i].m_fragment.size(); orientP = orientP / orientP.norm();

        m_site[i].m_orientN = orientN;
        m_site[i].m_orientP = orientP;

    }

    std::cout << "Done" << '\n' ;
}

void Graph::findNeighbors()
{
    std::cout << "Looking for conjugated segments neighbors... ";
    double dThresholdInter = m_dThresholdInter;
    double dThresholdIntra = m_dThresholdIntra;
    double dThresholdDopantInter = m_dThresholdDopantInter;
    double dThresholdDopantIntra = m_dThresholdDopantIntra;

    /// Find conjugated segment neighbor based on a simple distance threshold between the rigid fragment of each conjugated segment ///
    bool periodic;
    
    // Init a random generator
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);


    for (int i(0) ; i < m_site.size() ; i ++)
    {
        //std::cout << i << " / " << m_site.size() << '\n';

        for (int j(0) ; j < m_site.size() ; j ++)
        {
            if(i != j)
            {

            std::vector<Vect3> pbcType; // Vector with the pbc type that linked both sites

            bool first(true); // Is there a bound between those conjugated segment already ?

            for( int k(0) ; k  < m_site[i].m_fragment.size() ; k++ )
            {
                for ( int l(0) ; l < m_site[j].m_fragment.size() ; l ++)
                {

                    // Calculate distance between the rigid fragments
                    Vect3 pbc;
                    Vect3 rfDisplacement( Vect3::displacementPbc((m_site[i].m_fragment[k]->m_center),(m_site[j].m_fragment[l]->m_center), m_system.m_boxSize, pbc  ) );
                    double d( rfDisplacement.norm() );

                    // Check if sites are linked by periodicity
                    if(pbc.m_x != 0 || pbc.m_y != 0 || pbc.m_z != 0)
                    {
                        periodic = true;
                    }
                    else
                    {
                        periodic = false;
                    }
                                      

                    // Calculate displacement vector between sites, with appropriate pbc
                    Vect3 siteDisplacement( Vect3::siteDisplacementPbc(m_site[i].m_pos,m_site[j].m_pos, m_system.m_boxSize, pbc) );

                    // Determine the hopping type
                    HopType h;
                    if( m_site[i].m_fragment[k]->m_molNbr ==  m_site[j].m_fragment[l]->m_molNbr)
                    {
                        h = intra;
                    }
                    else
                    {
                        h = inter;
                    }
		    /// BIMODAL TEMPERATURE DEPENDENT DISTRIBUTION ///
		    /*double T0(100.0);
		    double T(300.0);
		    double p_intra(exp(-T/T0));
		    if(dis(gen) <= p_intra)
		    {
		      h = intra;
		    }
		    else
		    {
		      h = inter;
		    }*/
		    

                    // Add the site to the neighbor list if needed
                    int diffFragNb = abs( m_site[i].m_fragment[k]->m_fragNbr - m_site[j].m_fragment[l]->m_fragNbr );
                    if( (d <= dThresholdInter && h == inter) || (d <= dThresholdIntra && h == intra) || (h == intra && dThresholdIntra == -1 && diffFragNb == 1) )
                    {

                        // Determine if this is the first time that we link both sites with this pbc type
                        first = true;
                        for (int p(0) ; p < pbcType.size() ; p ++)
                        {
                            if(pbc == pbcType[p]) {first = false;}
                        }


                        if(first)
                        {


                            m_site[i].m_neighbor.push_back(&m_site[j]);
                            m_site[i].m_neighborIdx.push_back(j);
                            m_site[i].m_periodic.push_back(periodic);

                            std::vector<Bound> bVect;
                            Bound b;
                            b.m_id1 = k;
                            b.m_id2 = l;
                            bVect.push_back(b);
                            m_site[i].m_boundFragment.push_back(bVect);

                            std::vector<HopType> hVect;
                            hVect.push_back(h);
                            m_site[i].m_hopType.push_back(hVect);

                            std::vector<int> iVect;
                            iVect.push_back(-1);
                            m_site[i].m_ionBridge.push_back(iVect);

                            m_site[i].m_neighborPbc.push_back(pbc);
                            m_site[i].m_neighborDisplacement.push_back(siteDisplacement);

                            pbcType.push_back(pbc);


                        }
                        else
                        {
                            Bound b;
                            b.m_id1 = k;
                            b.m_id2 = l;
                            m_site[i].m_boundFragment[m_site[i].m_boundFragment.size()-1].push_back(b);

                            m_site[i].m_hopType[m_site[i].m_hopType.size()-1].push_back(h);

                            m_site[i].m_ionBridge[m_site[i].m_ionBridge.size()-1].push_back(-1);

                            m_site[i].m_periodic[m_site[i].m_periodic.size()-1] = m_site[i].m_periodic[m_site[i].m_periodic.size()-1] && periodic;

                        }

                    }
                }

            }


            }

        }
    }


    /// Find conjugated segment neighbors based on sphere search around counterions ///
    /*for (int i(0) ; i < m_system.m_molecule.size() ; i ++ )
    {
        std::vector<RigidFragment*> IonNeighbor;
        std::vector<double> distance;

        // Look for all rigid fragments in the sphere of radius dTresholdIons around the counter ion
        for (int j(0) ; j < m_site.size() ; j ++)
        {
            double dMin(1e10);
            bool isNeighboor(false);

            for( int k(0) ; k  < m_site[j].m_fragment.size() ; k++ )
            {
                double d(Vect3::distancePbc(m_system.m_molecule[i].m_center, m_site[j].m_fragment[k]->m_center, m_system.m_boxSize, periodic));

                if( d <= std::max(dThresholdDopantInter,dThresholdDopantIntra) )
                {
                    IonNeighbor.push_back(m_site[j].m_fragment[k]);
                    isNeighboor = true;
                    if(d < dMin) {dMin = d;}

                }
            }

            if(isNeighboor)
            {
                distance.push_back(dMin);
            }

        }

        // Update neighbor list
        for (int j(0) ; j < IonNeighbor.size() ; j ++)
        {
            for( int k(0) ; k  < IonNeighbor.size() ; k++ )
            {
                if( IonNeighbor[j]->m_siteNbr !=  IonNeighbor[k]->m_siteNbr )
                {
                    int site1(IonNeighbor[j]->m_siteNbr);
                    int site2(IonNeighbor[k]->m_siteNbr);

                    // Check if sites are linked by periodicity
                    Vect3::distancePbc(m_site[site1].m_pos, m_site[site2].m_pos , m_system.m_boxSize, periodic);

                    // Conjugated segments are already linked ?
                    bool first(true);
                    int idx;
                    for (int l(0) ; l < m_site[site1].m_neighborIdx.size() ; l++ )
                    {
                        if (m_site[site1].m_neighborIdx[l] == site2 )
                        {
                            first = false;
                            idx = l;

                        }
                    }

                    HopType h;
                    if( IonNeighbor[j]->m_molNbr ==  IonNeighbor[k]->m_molNbr)
                    {
                        h = dopantIntra;
                    }
                    else
                    {
                        h = dopantInter;
                    }

                    if( (distance[j] <= dThresholdDopantInter && distance[k] <= dThresholdDopantInter && h == dopantInter) || (distance[j] <= dThresholdDopantIntra && distance[k] <= dThresholdDopantIntra && h == dopantIntra) )
                    {
                        // Add links
                        if(first)
                        {
                            m_site[site1].m_neighbor.push_back(&m_site[site2]);
                            m_site[site1].m_neighborIdx.push_back(site2);
                            m_site[site1].m_periodic.push_back(periodic);

                            std::vector<Bound> bVect;
                            Bound b;
                            b.m_id1 = IonNeighbor[j]->m_fragNbr;
                            b.m_id2 = IonNeighbor[k]->m_fragNbr;
                            bVect.push_back(b);
                            m_site[site1].m_boundFragment.push_back(bVect);

                            std::vector<HopType> hVect;
                            hVect.push_back(h);
                            m_site[site1].m_hopType.push_back(hVect);

                            std::vector<int> iVect;
                            iVect.push_back(i);
                            m_site[site1].m_ionBridge.push_back(iVect);

                        }
                        else
                        {
                            Bound b;
                            b.m_id1 = IonNeighbor[j]->m_fragNbr;
                            b.m_id2 = IonNeighbor[k]->m_fragNbr;
                            m_site[site1].m_boundFragment[idx].push_back(b);

                            m_site[site1].m_hopType[idx].push_back(h);

                            m_site[site1].m_ionBridge[idx].push_back(i);

                            m_site[site1].m_periodic[idx] = m_site[site1].m_periodic[idx] && periodic;

                        }

                        first = false;
                    }

                }
            }
        }


    }*/

    /// Remove wrong neighbors, i.e. neighbors that are influenced by other polymer chain ///
    /*int removedEdges(0);
    if (m_sphereInfluenceRad > 0)
    {
        #pragma omp parallel for reduction(+:removedEdges)
        for (int i = 0 ;  i < m_site.size() ; i++ )
        {
            //std::cout << i << " / " << m_site.size() << '\n';

            for (int j(0) ; j < m_site[i].m_boundFragment.size() ; j++)
            {
                bool realNeighbors(true);

                for (int k(0) ; k < m_site[i].m_boundFragment[j].size() ; k++)
                {
                    int rf1( m_site[i].m_boundFragment[j][k].m_id1 );
                    int rf2( m_site[i].m_boundFragment[j][k].m_id2 );

                    for(int n(0) ; n < m_site.size() ; n++)
                    {
                        if(n != i && n != m_site[i].m_neighborIdx[j])
                        {
                            for (int p(0) ; p < m_site[i].m_fragment.size() ; p++ )
                            {
                                realNeighbors = realNeighbors && !crossingRaySphere(m_site[i].m_fragment[rf1]->m_center ,m_site[i].m_neighbor[j]->m_fragment[rf2]->m_center,m_site[n].m_fragment[p]->m_center,m_sphereInfluenceRad, m_system.m_boxSize);
                            }
                        }
                    }

                }

                if(!realNeighbors)
                {
                    m_site[i].m_neighbor.erase(m_site[i].m_neighbor.begin() + j);
                    m_site[i].m_neighborIdx.erase(m_site[i].m_neighborIdx.begin() + j);
                    m_site[i].m_periodic.erase(m_site[i].m_periodic.begin() + j);
                    m_site[i].m_boundFragment.erase(m_site[i].m_boundFragment.begin() + j);
                    m_site[i].m_hopType.erase(m_site[i].m_hopType.begin() + j);
                    m_site[i].m_ionBridge.erase(m_site[i].m_ionBridge.begin() + j);

                    // Dangerous with long chain...
                    m_site[i].m_neighborDisplacement.erase(m_site[i].m_neighborDisplacement.begin() + j);
                    m_site[i].m_neighborPbc.erase(m_site[i].m_neighborPbc.begin() + j);

                    j --;

                    removedEdges ++;
                }
            }
        }
    }
    std::cout << "  " << removedEdges << " edges removed..." ;*/





    m_edgeNbAll = 0;
    m_edgeNbIntra = 0;
    m_edgeNbInter = 0;
    m_edgeNbDopantIntra = 0;
    m_edgeNbDopantInter = 0;

    m_edgeNbAll_pbc = 0;
    m_edgeNbIntra_pbc = 0;
    m_edgeNbInter_pbc = 0;
    m_edgeNbDopantIntra_pbc = 0;
    m_edgeNbDopantInter_pbc = 0;

    for (int i(0) ; i < m_site.size() ; i++)
    {
        for (int j(0) ; j < m_site[i].m_neighbor.size() ; j ++)
        {

            bool edgeAll(false);
            bool edgeIntra(false);
            bool edgeInter(false);
            bool edgeDopantIntra(false);
            bool edgeDopantInter(false);

            for (int k(0) ; k < m_site[i].m_hopType[j].size() ; k ++)
            {
                edgeAll = true;
                edgeIntra = edgeIntra || ( m_site[i].m_hopType[j][k] == intra );
                edgeInter = edgeInter || ( m_site[i].m_hopType[j][k] == inter );
                edgeDopantIntra = edgeDopantIntra || ( m_site[i].m_hopType[j][k] == dopantIntra );
                edgeDopantInter = edgeDopantInter || ( m_site[i].m_hopType[j][k] == dopantInter );
            }

            if (edgeAll)
            {
                m_edgeNbAll ++;
                if(!m_site[i].m_periodic[j]) { m_edgeNbAll_pbc ++; }
            }
            if (edgeIntra)
            {
                m_edgeNbIntra ++;
                if(!m_site[i].m_periodic[j]) { m_edgeNbIntra_pbc ++; }
            }
            if (edgeInter)
            {
                m_edgeNbInter ++;
                if(!m_site[i].m_periodic[j]) { m_edgeNbInter_pbc ++; }
            }
            if (edgeDopantIntra)
            {
                m_edgeNbDopantIntra++;
                if(!m_site[i].m_periodic[j]) { m_edgeNbDopantIntra_pbc ++; }
            }
            if (edgeDopantInter)
            {
                m_edgeNbDopantInter++;
                if(!m_site[i].m_periodic[j]) { m_edgeNbDopantInter_pbc ++; }
            }
        }

    }

    std::cout << "Done" << '\n' ;

}

void Graph::writeGraph(std::string fileNameGraphPdb, std::string fileNameGraphPsf, bool periodic)
{
    std::cout << "Writing the graph... ";
    Atom a;

    std::ofstream outFlux(fileNameGraphPdb.c_str());
    if (outFlux.is_open())
    {
        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_system.m_boxSize.m_x;
        outFlux << std::setw(9) << m_system.m_boxSize.m_y;
        outFlux << std::setw(9) << m_system.m_boxSize.m_z;
        outFlux << m_system.m_headerLine;
        outFlux << '\n';


        for (int i(0) ; i < m_site.size() ; i ++ )
        {
            a.m_type = C;
            a.m_pos = m_site[i].m_pos;

            outFlux << "ATOM";
            outFlux << std::setw(7);
            outFlux << i+1;
            outFlux << std::setw(3);
            outFlux << a.m_type;
            outFlux << std::setw(8);
            outFlux << "RES";
            outFlux << std::setw(4);
            outFlux << i+1;
            outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

            outFlux << a.m_pos.m_x;
            outFlux << std::setw(8);
            outFlux << a.m_pos.m_y;
            outFlux << std::setw(8);
            outFlux << a.m_pos.m_z;
            outFlux << std::fixed << std::setprecision(2) << std::setw(6);

            //outFlux << log10(m_site[i].m_occup);
            outFlux << m_site[i].m_energy;
            outFlux << std::fixed << std::setprecision(2) << std::setw(6);
            outFlux << (double) m_site[i].m_clusterId;
            outFlux << '\n';


        }

        outFlux.close();

    }

    else
    {
        std::cerr << "Unable to open graph pdb file" << std::endl;
    }

    for (int n(1) ; n < 6 ; n ++)
    {
        std::string fileType;
        int edgeNb;

        if(periodic)
        {
            switch (n)
            {
                case 1 : fileType = "All"; edgeNb = m_edgeNbAll;
                break;
                case 2 : fileType = "Intra"; edgeNb = m_edgeNbIntra;
                break;
                case 3 : fileType = "Inter"; edgeNb = m_edgeNbInter;
                break;
                case 4 : fileType = "DopantIntra"; edgeNb = m_edgeNbDopantIntra;
                break;
                case 5 : fileType = "DopantInter"; edgeNb = m_edgeNbDopantInter;
                break;
            }
        }
        else
        {
            switch (n)
            {
                case 1 : fileType = "All"; edgeNb = m_edgeNbAll_pbc;
                break;
                case 2 : fileType = "Intra"; edgeNb = m_edgeNbIntra_pbc;
                break;
                case 3 : fileType = "Inter"; edgeNb = m_edgeNbInter_pbc;
                break;
                case 4 : fileType = "DopantIntra"; edgeNb = m_edgeNbDopantIntra_pbc;
                break;
                case 5 : fileType = "DopantInter"; edgeNb = m_edgeNbDopantInter_pbc;
                break;
            }
        }

        std::string nameTmp(fileNameGraphPsf);
        nameTmp.pop_back();
        nameTmp.pop_back();
        nameTmp.pop_back();
        nameTmp.pop_back();
        nameTmp += fileType + ".psf";

        std::ofstream outFlux2(nameTmp.c_str());
        if (outFlux2.is_open())
        {
            outFlux2 << "PSF" << '\n' << '\n' ;
            outFlux2 << "       1 !NTITLE" << '\n';
            outFlux2 << " REMARKS connectivity of the graph" << '\n' << '\n';
            outFlux2 << std::setw(8);
            outFlux2 << m_site.size();
            outFlux2 << std::setw(7);
            outFlux2 << "!NATOM" << '\n' ;


            for (int i(0) ; i < m_site.size() ; i ++ )
            {
                outFlux2 << std::setw(8);
                outFlux2 << i+1;
                outFlux2 << std::setw(4);
                outFlux2 << "RES";
                outFlux2 << std::setw(3);
                outFlux2 << 1;
                outFlux2 << std::setw(7);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
                outFlux2 << 1.0;
                outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
                outFlux2 << 1.0;
                outFlux2 << std::setw(12);
                outFlux2 << 0;
                outFlux2 << '\n';

            }

            outFlux2 << '\n';
            outFlux2 << std::setw(8);
            outFlux2 << edgeNb;
            outFlux2 << std::setw(14);
            outFlux2 << "!NBOND: bonds" << '\n' ;

            int counter(0);
            for ( int i(0) ; i < m_site.size() ; i ++ )
            {
                for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++ )
                {
                    bool toWrite(false);

                    for (int k(0) ; k < m_site[i].m_hopType[j].size() ; k ++)
                    {
                        switch (n)
                        {
                            case 1 : toWrite = true;
                            break;
                            case 2 : toWrite = toWrite || ( m_site[i].m_hopType[j][k] == intra );
                            break;
                            case 3 : toWrite = toWrite || ( m_site[i].m_hopType[j][k] == inter );
                            break;
                            case 4 : toWrite = toWrite || ( m_site[i].m_hopType[j][k] == dopantIntra );
                            break;
                            case 5 : toWrite = toWrite || ( m_site[i].m_hopType[j][k] == dopantInter );
                            break;
                        }

                    }

                    if( (toWrite && !periodic && !m_site[i].m_periodic[j]) || (toWrite && periodic) )
                    {
                        outFlux2 << std::setw(8);
                        outFlux2 << i+1;
                        outFlux2 << std::setw(8);
                        outFlux2 << m_site[i].m_neighborIdx[j]+1;
                        counter ++;

                        if (counter % 4 == 0) { outFlux2 << '\n'; }
                    }

                }
            }

            outFlux2.close();
        }
        else
        {
            std::cerr << "Unable to open graph psf file" << std::endl;
        }

    }

    std::cout << "done" << '\n';

}

void Graph::cluster(double threshold, int clusteringMode)
{
    std::cout << "Clustering the graph... ";
    m_clusterSiteNb.erase( m_clusterSiteNb.begin(),  m_clusterSiteNb.end());
    m_clusterDim.erase(m_clusterDim.begin(), m_clusterDim.end());
    m_clusterNb = 0;

    /// Cluster the graph based on the threshold value, and store a clusterId for each site ///
    std::vector<int> siteGlob;
    int clusterId(0);

    while (siteGlob.size() < m_site.size())
    {
        bool isPresent(true);
        std::vector<int> site;

        int counter(0);
        while (isPresent && counter < m_site.size())
        {
            isPresent = false;

            for (int j(0) ; j < siteGlob.size() ; j ++)
            {
                isPresent = isPresent || (siteGlob[j] == counter);
            }

            counter ++;
        }

        site.push_back(counter - 1);

        int idx(0);

        while(idx < site.size())
        {
            int idxSite(site[idx]);

            for (int i(0) ; i < m_site[idxSite].m_neighborIdx.size() ; i ++)
            {
                if( ( clusteringMode == 1 && (m_site[idxSite].m_hopTo[i] >= threshold || m_site[idxSite].m_hopFrom[i] >= threshold) ) || clusteringMode == 0 )
                {
                    isPresent = false;

                    for (int j(0) ; j < site.size() ; j ++)
                    {
                        isPresent = isPresent || (m_site[idxSite].m_neighborIdx[i] == site[j]);
                    }

                    if(!isPresent)
                    {
                        site.push_back(m_site[idxSite].m_neighborIdx[i]);
                    }
                }
            }

            idx++;

        }

        for (int i(0) ; i < site.size() ; i++)
        {
            m_site[site[i]].m_clusterId = clusterId;
            siteGlob.push_back(site[i]);
        }


        m_clusterSiteNb.push_back(site.size());
        clusterId ++;

    }

    m_clusterNb = clusterId;
    std::cout << "Done (" << clusterId <<" clusters found)" << '\n';
}


double Graph::percolation()
{

    std::cout << "Compute percolation ratio... ";
    m_clusterDim.resize(0);

    /// Define percolation ratio as the bigger cluster dimension (x, y or z) divided by the box size ///
    for (int i(0) ; i < m_clusterNb ; i ++)
    {
        m_clusterDim.push_back(Vect3(0,0,0));
    }

    for (int i(0) ; i < m_clusterNb ; i ++)
    {
        //std::cout << "i = " << i << '\n';

        Direction d;

        for (int j(0) ; j < 3 ; j ++)
        {
            switch(j)
            {
                case 0 : d = x;
                break;
                case 1 : d = y;
                break;
                case 2 : d = z;
                break;
            }

            std::vector<std::vector<Edge>> g = createJohnsonGraph(i, d);


            // run Johnson's algorithm to get all pairs shortest paths
            AllSP asp = johnson(g);

            // find the "bigger shortest path", ie. the path with highest displacement,
            // amongst all shortest path pairs.
            double bigger = -1e30;
            for (int k = 1; k < g.size(); k++)
            {
                for (int l = 1; l < g.size(); l++)
                {
                    if (fabs(asp[k][l]) > bigger)
                    {
                       bigger = fabs(asp[k][l]);
                    }
                }

            }


            switch(j)
            {
                case 0 : m_clusterDim[i].m_x = bigger;
                break;
                case 1 : m_clusterDim[i].m_y = bigger;
                break;
                case 2 : m_clusterDim[i].m_z = bigger;
                break;
            }


          //std::cout << "Bigger shortest path = " << bigger << std::endl;
        }

    }

    std::cout << "Done" << '\n';

}


std::vector<std::vector<Edge>> Graph::createJohnsonGraph(int clusterId, Direction d)
{
  std::vector<std::vector<Edge>> g;
  int n(m_clusterSiteNb[clusterId]);
  g.resize(n+1);

  std::vector<int> site;

  for(int i(0) ; i < m_site.size() ; i ++)
  {
       if (m_site[i].m_clusterId == clusterId)
       {
           site.push_back(i);
           //if(clusterId==75){std::cout << "ghh  " << i << '\n'; }

       }
  }

  //std::cout << "here    " <<  n << "   " << site.size() << "    " << clusterId << "     " << m_clusterNb << '\n';


  for(int i(0) ; i < site.size() ; i ++)
  {
      int idx1(site[i]);
      int v1;
      int v2;
      double c;

      for (int j(0) ; j < site.size() ; j++ )
      {
          int idx2(site[j]);

          bool isNeighboor(false);

          std::vector<Vect3> displacementVect;

          for (int k(0) ; k < m_site[idx1].m_neighborIdx.size() ; k++)
          {
              isNeighboor = isNeighboor || (m_site[idx1].m_neighborIdx[k] == idx2);
              if(m_site[idx1].m_neighborIdx[k] == idx2) {displacementVect.push_back(m_site[idx1].m_neighborDisplacement[k]);}

          }

          if ( isNeighboor )
          {

              int v1=i+1;
              int v2=j+1;
              c = 0;

              for (int p(0) ; p < displacementVect.size() ; p ++)
              {

              Vect3 displacement = displacementVect[p];

              switch(d)
              {
                    case 0 : c = displacement.m_x;
                    break;
                    case 1 : c = displacement.m_y;
                    break;
                    case 2 : c = displacement.m_z;
                    break;
              }

              assert(g.size()>v1);
              g[v1].push_back({v2, c});
              //std::cout << v1 << "  " << v2 << "    " << c << '\n';

              }
          }



      }
  }


  return g;
}

void Graph::computeRate(RateType rate)
{
    std::cout << "Computing hopping rates... ";

    for (int i(0) ; i < m_site.size() ; i++)
    {
        m_site[i].m_hopTo.resize(0);
        m_site[i].m_hopFrom.resize(0);

        for(int j(0) ; j < m_site[i].m_neighbor.size() ; j++)
        {
            HopType h(m_site[i].m_hopType[j][0]);

            bool periodic;
            Vect3 displacement = m_site[i].m_neighborDisplacement[j];
            doubleType d(displacement.norm());
            if(rate == maPheno)
            {
                doubleType distanceMin(1e8);
                for(int k(0) ; k < m_site[i].m_fragment.size() ; k++)
                {
                    bool periodic(true);
                    int id1;
                    int id2;

                    for(int l(0) ; l < m_site[i].m_neighbor[j]->m_fragment.size() ; l++)
                    {
                        doubleType dd(Vect3::distancePbc(m_site[i].m_fragment[k]->m_center,(m_site[i].m_neighbor[j]->m_fragment[l])->m_center,m_system.m_boxSize,periodic));
                        if( dd < distanceMin)
                        {
                            distanceMin = dd;
                            id1 = k;
                            id2 = l;
                        }
                    }
                }
                d = distanceMin;
            }

            doubleType hopTo(0);
            doubleType hopFrom(0);
            doubleType transfer(m_site[i].m_transferIntegral[j][0]);
	    
	    int nbSiteNoLevel(m_site.size()/m_nbrLevels);
	    int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );
	    int idxNeighbor( m_site[i].m_neighborIdx[j] );
	    int idOrb2( (idxNeighbor - (idxNeighbor % nbSiteNoLevel)) / nbSiteNoLevel );
	    
	    //TRY
     /*if(displacement.m_y != 0)
     {
	      if(idOrb1 != idOrb2 && idOrb1 == 0)
	      {
		transfer = -0.01;
	      }
	      if(idOrb1 != idOrb2 && idOrb1 == 1)
	      {
		transfer = -0.02;
	      }
     }*/
	    
            Vect3 neighborPbc;
	    
	    doubleType carrierFactor(0);
	    if(m_levelType == "HOMO") {carrierFactor = -1;}
	    if(m_levelType == "LUMO") {carrierFactor = 1;}
	    doubleType dE(m_site[i].m_neighbor[j]->m_energy - m_site[i].m_energy);
	    doubleType dEnergy_ij(carrierFactor * dE + carrierFactor * (m_field * displacement));
	    //doubleType dEnergy_ji(-carrierFactor * dE + carrierFactor * (m_field * displacement));
	    doubleType dEnergy_ji(-dEnergy_ij);
	    
	    doubleType w0;
	    doubleType decayTransfer;
	    if( h == inter )
	    {
	      w0 = m_w0_inter;
	      if(rate == maPheno) {decayTransfer = exp( - 2.0 * d / m_loc_inter );}
	      if(rate == maTransfer) {decayTransfer = transfer*transfer;}
	    }
	    if( h == intra )
	    {
	      w0 = m_w0_intra;
	      if(rate == maPheno) {decayTransfer = exp( - 2.0 * d / m_loc_intra );}
	      if(rate == maTransfer) {decayTransfer = transfer*transfer;}
	    }
	    if( h == level )
	    {
	      w0 = m_w0_level;
	      decayTransfer = transfer*transfer;
	    }
	    
	    if ( dEnergy_ij >= 0 )
	    {
		hopTo += ( m_w0_inter * decayTransfer * exp( - dEnergy_ij / (kB * m_temperature) ) );
		hopFrom += ( m_w0_inter * decayTransfer );
	    }
	    else
	    {
		hopTo += ( m_w0_inter * decayTransfer );
		hopFrom += ( m_w0_inter * decayTransfer * exp( - dEnergy_ji / (kB * m_temperature) ) );
	    }

            m_site[i].m_hopTo.push_back(hopTo);
            m_site[i].m_hopFrom.push_back(hopFrom);
	    
	    //std::cout << hopTo << "	" << hopFrom << "	" << m_field.norm() << "	" << displacement.norm() << std::endl;

        }
    }
    
    std::cout << "Done" << std::endl;
}


void Graph::writeCurrent(std::string fileNameOccupancyPdb, std::string fileNameCurrentPdb, std::string fileNameCurrentPsf, bool periodic)
{
    std::cout << "Writing the current distribution... ";
    Atom a;

    std::ofstream outFlux0(fileNameOccupancyPdb.c_str());
    if (outFlux0.is_open())
    {

        outFlux0 << "CRYST1";
        outFlux0 << std::setw(9) << m_system.m_boxSize.m_x;
        outFlux0 << std::setw(9) << m_system.m_boxSize.m_y;
        outFlux0 << std::setw(9) << m_system.m_boxSize.m_z;
        outFlux0 << m_system.m_headerLine;
        outFlux0 << '\n';

        for (int i(0) ; i < m_site.size() ; i ++ )
        {
            a.m_type = C;
            a.m_pos = m_site[i].m_pos;

            outFlux0 << "ATOM";
            outFlux0 << std::setw(7);
            outFlux0 << i+1;
            outFlux0 << std::setw(3);
            outFlux0 << a.m_type;
            outFlux0 << std::setw(8);
            outFlux0 << "RES";
            outFlux0 << std::setw(4);
            outFlux0 << i+1;
            outFlux0 << std::fixed << std::setprecision(3) << std::setw(12) ;

            outFlux0 << a.m_pos.m_x;
            outFlux0 << std::setw(8);
            outFlux0 << a.m_pos.m_y;
            outFlux0 << std::setw(8);
            outFlux0 << a.m_pos.m_z;
            outFlux0 << std::fixed << std::setprecision(2) << std::setw(6);

            outFlux0 << m_site[i].m_occup;
            outFlux0 << std::fixed << std::setprecision(2) << std::setw(6);
            outFlux0 << (double) m_site[i].m_clusterId;
            outFlux0 << '\n';


        }

        outFlux0.close();

    }

    else
    {
        std::cerr << "Unable to open occupancy pdb file" << std::endl;
    }

    int nbEdges(0);
    std::ofstream outFlux(fileNameCurrentPdb.c_str());
    int counterG(0);
    if (outFlux.is_open())
    {

        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_system.m_boxSize.m_x;
        outFlux << std::setw(9) << m_system.m_boxSize.m_y;
        outFlux << std::setw(9) << m_system.m_boxSize.m_z;
        outFlux << m_system.m_headerLine;
        outFlux << '\n';


        int counter(0);
        for (int i(0) ; i < m_site.size() ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++ )
            {
                a.m_type = C;
                a.m_pos = m_site[i].m_pos;
                int k(m_site[i].m_neighborIdx[j]);

                //double current(m_site[i].m_occup * ( 1 - m_site[k].m_occup ) * m_site[i].m_hopTo[j] - m_site[k].m_occup * ( 1 - m_site[i].m_occup ) * m_site[i].m_hopFrom[j] );
                double current(std::max(fabs(m_site[i].m_hopTo[j]),fabs(m_site[i].m_hopFrom[j])));
                current = log(current * 1000.0) / log(10.0) ;

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter+1;
                outFlux << std::setw(3);
                outFlux << a.m_type;
                outFlux << std::setw(8);
                outFlux << "RES";
                outFlux << std::setw(4);
                outFlux << 0;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << a.m_pos.m_x;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_y;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_z;
                outFlux << std::setw(6) << std::fixed << std::setprecision(2) << std::setw(6);

                outFlux << (float) current ;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux <<  m_site[i].m_clusterId;
                outFlux << '\n';
                counter ++;



                a.m_type = C;
                a.m_pos = m_site[i].m_neighbor[j]->m_pos;

                outFlux << "ATOM";
                outFlux << std::setw(7);
                outFlux << counter+1;
                outFlux << std::setw(3);
                outFlux << a.m_type;
                outFlux << std::setw(8);
                outFlux << "RES";
                outFlux << std::setw(4);
                outFlux << 0;
                outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                outFlux << a.m_pos.m_x;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_y;
                outFlux << std::setw(8);
                outFlux << a.m_pos.m_z;
                outFlux << std::setw(6) << std::fixed << std::setprecision(2) << std::setw(6);

                outFlux << (float) current ;
                outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                outFlux << m_site[k].m_clusterId;
                outFlux << '\n';
                counter++;

                if ( periodic || ( !periodic &&  !m_site[i].m_periodic[j] ) )
                {
                    nbEdges ++;
                }

            }

        }
        counterG = counter;

        outFlux.close();

    }

    else
    {
        std::cerr << "Unable to open current pdb file" << std::endl;
    }

    std::ofstream outFlux2(fileNameCurrentPsf.c_str());
    if (outFlux2.is_open())
    {
        outFlux2 << "PSF" << '\n' << '\n' ;
        outFlux2 << "       1 !NTITLE" << '\n';
        outFlux2 << " REMARKS connectivity of the graph" << '\n' << '\n';
        outFlux2 << std::setw(8);
        outFlux2 << counterG;
        outFlux2 << std::setw(7);
        outFlux2 << "!NATOM" << '\n' ;

        int counter(0);
        for (int i(0) ; i < m_site.size() ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
            {
                outFlux2 << std::setw(8);
                outFlux2 << counter+1;
                outFlux2 << std::setw(4);
                outFlux2 << "RES";
                outFlux2 << std::setw(3);
                outFlux2 << 1;
                outFlux2 << std::setw(7);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
                outFlux2 << 1.0;
                outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
                outFlux2 << 1.0;
                outFlux2 << std::setw(12);
                outFlux2 << 0;
                outFlux2 << '\n';
                counter++;

                outFlux2 << std::setw(8);
                outFlux2 << counter+1;
                outFlux2 << std::setw(4);
                outFlux2 << "RES";
                outFlux2 << std::setw(3);
                outFlux2 << 1;
                outFlux2 << std::setw(7);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(5);
                outFlux2 << "RES";
                outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
                outFlux2 << 1.0;
                outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
                outFlux2 << 1.0;
                outFlux2 << std::setw(12);
                outFlux2 << 0;
                outFlux2 << '\n';
                counter++;
            }
        }

        outFlux2 << '\n';
        outFlux2 << std::setw(8);
        outFlux2 << nbEdges;
        outFlux2 << std::setw(14);
        outFlux2 << "!NBOND: bonds" << '\n' ;

        counter = 0;
        int counterReal(0);
        for ( int i(0) ; i < m_site.size() ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++ )
            {
                if( (!periodic && !m_site[i].m_periodic[j]) ||  periodic )
                {
                    outFlux2 << std::setw(8);
                    outFlux2 << counter+1;
                    outFlux2 << std::setw(8);
                    outFlux2 << counter+2;

                    counterReal += 2;
                    if ( (counterReal / 2) % 4 == 0) { outFlux2 << '\n'; }
                }

                    counter ++;
                    counter ++;

            }
        }

        outFlux2.close();
    }
    else
    {
        std::cerr << "Unable to open current psf file" << std::endl;
    }

    std::cout << "done" << '\n';

}


void Graph::gdm(double sigma, double mu, unsigned int seed, bool coarsed)
{
    // Init gaussian
    std::default_random_engine generator;
    generator.seed(seed);
    std::normal_distribution<double> distribution(mu,sigma);

    if(!coarsed)
    {
        for (int v(0) ; v < m_site.size() / m_nbrLevels / repetition ; v++)
        {
            for(int u(0) ; u < m_nbrLevels ; u++)
            {
                double energy(0);
		
                /// Draw energies from truncated normal distribution ///
                do
                {
                    energy = (distribution(generator) + m_levelEnergies[u]);
                }
                //Truncation
                //while(energy < sigma + m_levelEnergies[u]);
                //No truncation
                while(false);

                int i(u*m_site.size()/m_nbrLevels + m_orderMeanField[v]);

                m_site[i].m_energy = energy;

                for(int w(1) ; w <8 ; w ++)
                {
                    int q(i+w*(m_site.size() / m_nbrLevels / repetition));
                    m_site[q].m_energy = energy;
                }
            }
        }
    }
    else
    {
        for(int i(0) ; i < m_site.size() ; i++)
        {
                double energy(0);
                /// Draw energies from truncated normal distribution ///
                while(energy < sigma)
                {
                    energy = (distribution(generator) + m_levelEnergies[0]);
                }
                m_site[i].m_energy = energy;
        }
    }

}


void Graph::writeBoundGeometry(std::string fileNameGeometry)
{
    int nbBin(100);
    double binLength(90.0 / (double) nbBin);
    std::vector<std::vector<double>> hist;
    hist.resize(nbBin);
    for (int i(0) ; i < nbBin ; i++)
    {
        for (int j(0) ; j < nbBin ; j++)
        {
                hist[i].push_back(0);
        }
    }

    for(int i(0) ; i < m_site.size() ; i++)
    {
        for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
        {
            HopType h;
            h = intra;

            for (int k(0) ; k < m_site[i].m_hopType[j].size() ; k++)
            {
                if(m_site[i].m_hopType[j][k] == dopantInter) {h = dopantInter;}
            }

            //if(h == dopantInter )
            if(true)
            {

                Vect3 bound(m_site[i].m_neighbor[j]->m_pos - m_site[i].m_pos);
                bound = bound / bound.norm();

                double angleN(acos(m_site[i].m_orientN*bound)*180.0/pi);
                if(angleN > 90.0) {angleN = 180.0 - angleN;}
                double angleP(acos(m_site[i].m_orientP*bound)*180.0/pi);
                if(angleP > 90.0) {angleP= 180.0 - angleP;}

                /*
                double angleN(acos(m_site[i].m_orientN*m_site[i].m_neighbor[j]->m_orientN)*180.0/pi);
                if(angleN > 90.0) {angleN = 180.0 - angleN;}
                double angleP(acos(m_site[i].m_orientP*m_site[i].m_neighbor[j]->m_orientP)*180.0/pi);
                if(angleP > 90.0) {angleP= 180.0 - angleP;}
                */

                int idxI = floor(angleN / binLength);
                //int idxJ = floor(angleP / binLength);
                int idxJ(0);

                /*if(idxI == hist.size()-1  && idxJ == hist.size()-1)
                {
                    std::cout << "0 case : " << m_site[i].m_fragment[0]->m_molNbr << "      " << m_site[i].m_neighbor[j]->m_fragment[0]->m_molNbr << '\n';
                    std::cout << acos(m_site[i].m_orientN*m_site[i].m_orientP)*180.0/pi << '\n';
                    std::cout << m_site[i].m_orientN.m_x << "   " <<m_site[i].m_orientN.m_y  << "   " << m_site[i].m_orientN.m_z << '\n' ;
                    std::cout << m_site[i].m_orientP.m_x << "   " <<m_site[i].m_orientP.m_y  << "   " << m_site[i].m_orientP.m_z << '\n' ;

                    std::cout << acos(m_site[m_site[i].m_neighborIdx[j]].m_orientN*m_site[m_site[i].m_neighborIdx[j]].m_orientP)*180.0/pi << '\n';
                    std::cout << m_site[m_site[i].m_neighborIdx[j]].m_orientN.m_x << "   " <<m_site[m_site[i].m_neighborIdx[j]].m_orientN.m_y  << "   " << m_site[m_site[i].m_neighborIdx[j]].m_orientN.m_z << '\n' ;
                    std::cout << m_site[m_site[i].m_neighborIdx[j]].m_orientP.m_x << "   " <<m_site[m_site[i].m_neighborIdx[j]].m_orientP.m_y  << "   " << m_site[m_site[i].m_neighborIdx[j]].m_orientP.m_z << '\n' ;
                }*/

                hist[idxI][idxJ] = hist[idxI][idxJ] + 1;

            }
        }
    }

    std::ofstream of(fileNameGeometry.c_str());

    if(of)
    {
        of << "  ";

        for(int i(0) ; i < hist[0].size() ; i++)
        {
            of << i*binLength << "  ";
        }

        of << '\n';

        for(int i(0) ; i < hist.size() ; i++)
        {
            of << i*binLength << "  ";

            for (int j(0) ; j < hist[i].size() ; j++ )
            {
                of << hist[i][j] << "  ";

            }

            of << '\n';
        }
    }
    else
    {
        std::cerr << "Can't open geometry file" << '\n';
    }
}



void Graph::writeBoundGeometryRF(std::string fileNameGeometry)
{
    int nbBin(20);
    double binLength(90.0 / (double) nbBin);
    std::vector<std::vector<double>> hist;
    hist.resize(nbBin);
    for (int i(0) ; i < nbBin ; i++)
    {
        for (int j(0) ; j < nbBin ; j++)
        {
                hist[i].push_back(0);
        }
    }

    for(int i(0) ; i < m_site.size() ; i++)
    {
        for (int j(0) ; j < m_site[i].m_boundFragment.size() ; j++ )
        {
            for (int k(0) ; k < m_site[i].m_boundFragment[j].size() ; k++)
            {

                if(true)
                //if(m_site[i].m_hopType[j][k] == inter )
                {
                    int id1(m_site[i].m_boundFragment[j][k].m_id1);
                    int id2(m_site[i].m_boundFragment[j][k].m_id2);

                    Vect3 bound(m_site[i].m_neighbor[j]->m_fragment[id2]->m_center - m_site[i].m_fragment[id1]->m_center);
                    bound = bound / bound.norm();

                    double angleN(acos(m_site[i].m_fragment[id1]->m_orientN*bound)*180.0/pi);
                    if(angleN > 90.0) {angleN = 180.0 - angleN;}
                    double angleP(acos(m_site[i].m_fragment[id1]->m_orientP*bound)*180.0/pi);
                    if(angleP > 90.0) {angleP= 180.0 - angleP;}


                    int idxI = floor(angleN / binLength);
                    int idxJ = floor(angleP / binLength);
                    //idxJ = 0;

                    hist[idxI][idxJ] = hist[idxI][idxJ] + 1;

                    /*if(idxI == 0 )
                    {
                        std::cout << "0 case : " << m_site[i].m_fragment[id1]->m_molNbr << "      " << id1 << "     " << m_site[i].m_neighbor[j]->m_fragment[id2]->m_molNbr << "      " << id2 << '\n';

                    }*/

                }

            }
        }
    }

    std::ofstream of(fileNameGeometry.c_str());

    if(of)
    {
        of << "  ";

        for(int i(0) ; i < hist[0].size() ; i++)
        {
            of << i*binLength << "  ";
        }

        of << '\n';

        for(int i(0) ; i < hist.size() ; i++)
        {
            of << i*binLength << "  ";

            for (int j(0) ; j < hist[i].size() ; j++ )
            {
                of << hist[i][j] << "  ";

            }

            of << '\n';
        }

        // Reverse chart
        /*for(int i(hist[0].size()-1) ; i >= 0 ; i--)
        {
            of << i*binLength << "  ";
        }

        of << '\n';

        for(int i(hist.size() - 1) ; i >= 0 ; i--)
        {
            of << i*binLength << "  ";

            for (int j(hist[i].size() - 1) ; j >= 0 ; j-- )
            {
                of << hist[i][j] << "  ";

            }

            of << '\n';
        }*/
    }
    else
    {
        std::cerr << "Can't open geometry file" << '\n';
    }
}


void Graph::writeBoundGeometryRFDist(std::string fileNameGeometry)
{
    int nbBinI(100);
    int nbBinJ(100);
    double binLengthI(90.0 / (double) nbBinI);
    double binLengthJ(2*m_dThresholdDopantInter / (double) nbBinJ);
    std::vector<std::vector<double>> hist;
    hist.resize(nbBinI);
    for (int i(0) ; i < nbBinI ; i++)
    {
        for (int j(0) ; j < nbBinJ ; j++)
        {
                hist[i].push_back(0);
        }
    }

    for(int i(0) ; i < m_site.size() ; i++)
    {
        for (int j(0) ; j < m_site[i].m_boundFragment.size() ; j++ )
        {
            for (int k(0) ; k < m_site[i].m_boundFragment[j].size() ; k++)
            {

                if(true)
                //if(m_site[i].m_hopType[j][k] == inter )
                {
                    int id1(m_site[i].m_boundFragment[j][k].m_id1);
                    int id2(m_site[i].m_boundFragment[j][k].m_id2);

                    bool periodic;

                    Vect3 bound = Vect3::displacementPbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize);
                    bound = bound / bound.norm();

                    double angleN(acos(m_site[i].m_fragment[id1]->m_orientN*bound)*180.0/pi);
                    if(angleN > 90.0) {angleN = 180.0 - angleN;}


                    double d = Vect3::distancePbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize , periodic);

                    int idxI = floor(angleN / binLengthI);
                    int idxJ = floor(d / binLengthJ);

                    hist[idxI][idxJ] = hist[idxI][idxJ] + 1;

                }

            }
        }
    }

    std::ofstream of(fileNameGeometry.c_str());

    if(of)
    {
        of << "  ";

        for(int i(0) ; i < hist[0].size() ; i++)
        {
            of << i*binLengthJ << "  ";
        }

        of << '\n';

        for(int i(0) ; i < hist.size() ; i++)
        {
            of << i*binLengthI << "  ";

            for (int j(0) ; j < hist[i].size() ; j++ )
            {
                of << hist[i][j] << "  ";

            }

            of << '\n';
        }

    }
    else
    {
        std::cerr << "Can't open geometry file" << '\n';
    }
}


void Graph::writeFragmentVector(std::string fileNameFragmentVector)
{
    std::ofstream of(fileNameFragmentVector.c_str());

    if(of)
    {
        for(int i(0) ; i < m_system.m_polymer.size() ; i++)
        {
            for (int j(0) ; j < m_system.m_polymer[i].m_fragment.size() ; j++ )
            {
                Vect3 v1(m_system.m_polymer[i].m_fragment[j].m_center);
                Vect3 v2(m_system.m_polymer[i].m_fragment[j].m_orientN);
                Vect3 v3(m_system.m_polymer[i].m_fragment[j].m_orientP);
                Vect3 v4(m_system.m_polymer[i].m_fragment[j].m_orientT);

                // Write orient N
                of << "graphics top color 0" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 0" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} {" <<
                v1.m_x + v2.m_x*m_vectLength*1.20 << " " << v1.m_y + v2.m_y*m_vectLength*1.20 << " " << v1.m_z + v2.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';

                // Write orient P
                of << "graphics top color 1" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 1" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} {" <<
                v1.m_x + v3.m_x*m_vectLength*1.20 << " " << v1.m_y + v3.m_y*m_vectLength*1.20 << " " << v1.m_z + v3.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';

                // Write orient T
                of << "graphics top color 7" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v4.m_x*m_vectLength << " " << v1.m_y + v4.m_y*m_vectLength << " " << v1.m_z + v4.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 7" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v4.m_x*m_vectLength << " " << v1.m_y + v4.m_y*m_vectLength << " " << v1.m_z + v4.m_z*m_vectLength << "} {" <<
                v1.m_x + v4.m_x*m_vectLength*1.20 << " " << v1.m_y + v4.m_y*m_vectLength*1.20 << " " << v1.m_z + v4.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';


            }
        }

        for(int i(0) ; i < m_system.m_molecule.size() ; i++)
        {

            Vect3 v1(m_system.m_molecule[i].m_center);
            Vect3 v2(m_system.m_molecule[i].m_orientN);
            Vect3 v3(m_system.m_molecule[i].m_orientP);

            // Write orient N
            of << "graphics top color 0" << '\n';
            of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
            v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} radius " <<
            m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

            of << "graphics top color 0" << '\n';
            of << "graphics top cone {" <<
            v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} {" <<
            v1.m_x + v2.m_x*m_vectLength*1.20 << " " << v1.m_y + v2.m_y*m_vectLength*1.20 << " " << v1.m_z + v2.m_z*m_vectLength*1.20 << "} radius " <<
            2*m_vectRadius << " resolution " << m_vectResol << '\n';

            // Write orient P
            of << "graphics top color 1" << '\n';
            of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
            v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} radius " <<
            m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

            of << "graphics top color 1" << '\n';
            of << "graphics top cone {" <<
            v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} {" <<
            v1.m_x + v3.m_x*m_vectLength*1.20 << " " << v1.m_y + v3.m_y*m_vectLength*1.20 << " " << v1.m_z + v3.m_z*m_vectLength*1.20 << "} radius " <<
            2*m_vectRadius << " resolution " << m_vectResol << '\n';
        }

        of.close();

    }
    else
    {
        std::cerr << "Can't open script file for writing rigid fragment orientations" << '\n';
    }
}

void Graph::readFragmentVector(std::string fileNameFragmentVector)
{
    std::ifstream of(fileNameFragmentVector.c_str());
    std::string dummyString;
    double dummyDouble;
    
    if(of)
    {
        for(int i(0) ; i < m_system.m_polymer.size() ; i++)
        {
            for (int j(0) ; j < m_system.m_polymer[i].m_fragment.size() ; j++ )
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
		
		v2 = v2 - v1; v2 = v2 / m_vectLength; 
		v3 = v3 - v1; v3 = v3 / m_vectLength; 
		v4 = v4 - v1; v4 = v4 / m_vectLength; 
                
                m_system.m_polymer[i].m_fragment[j].m_center = v1;
                m_system.m_polymer[i].m_fragment[j].m_orientN = v2;
                m_system.m_polymer[i].m_fragment[j].m_orientP = v3;
                m_system.m_polymer[i].m_fragment[j].m_orientT = v4;

            }
        }

        for(int i(0) ; i < m_system.m_molecule.size() ; i++)
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
	    
	    v2 = v2 - v1; v2 = v2 / m_vectLength; 
	    v3 = v3 - v1; v3 = v3 / m_vectLength; 
	    v4 = v4 - v1; v4 = v4 / m_vectLength; 
	    
	    m_system.m_molecule[i].m_center = v1;
	    m_system.m_molecule[i].m_orientN = v2;
	    m_system.m_molecule[i].m_orientP = v3;
	    //m_system.m_molecule[i].m_orientT = v4;
	    
        }

        of.close();

    }
    else
    {
        std::cerr << "Can't open script file for reading rigid fragment orientations" << '\n';
    }
}

void Graph::writeSiteVector(std::string fileNameSiteVector)
{
    std::ofstream of(fileNameSiteVector.c_str());

    if(of)
    {
        for(int i(0) ; i < m_site.size() ; i++)
        {

                Vect3 v1(m_site[i].m_pos);
                Vect3 v2(m_site[i].m_orientN);
                Vect3 v3(m_site[i].m_orientP);

                // Write orient N
                of << "graphics top color 0" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 0" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} {" <<
                v1.m_x + v2.m_x*m_vectLength*1.20 << " " << v1.m_y + v2.m_y*m_vectLength*1.20 << " " << v1.m_z + v2.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';

                // Write orient P
                of << "graphics top color 1" << '\n';
                of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} radius " <<
                m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

                of << "graphics top color 1" << '\n';
                of << "graphics top cone {" <<
                v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} {" <<
                v1.m_x + v3.m_x*m_vectLength*1.20 << " " << v1.m_y + v3.m_y*m_vectLength*1.20 << " " << v1.m_z + v3.m_z*m_vectLength*1.20 << "} radius " <<
                2*m_vectRadius << " resolution " << m_vectResol << '\n';

        }

        /*for(int i(0) ; i < m_system.m_tos.size() ; i++)
        {

            Vect3 v1(m_system.m_tos[i].m_center);
            Vect3 v2(m_system.m_tos[i].m_orientN);
            Vect3 v3(m_system.m_tos[i].m_orientP);

            // Write orient N
            of << "graphics top color 0" << '\n';
            of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
            v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} radius " <<
            m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

            of << "graphics top color 0" << '\n';
            of << "graphics top cone {" <<
            v1.m_x + v2.m_x*m_vectLength << " " << v1.m_y + v2.m_y*m_vectLength << " " << v1.m_z + v2.m_z*m_vectLength << "} {" <<
            v1.m_x + v2.m_x*m_vectLength*1.20 << " " << v1.m_y + v2.m_y*m_vectLength*1.20 << " " << v1.m_z + v2.m_z*m_vectLength*1.20 << "} radius " <<
            2*m_vectRadius << " resolution " << m_vectResol << '\n';

            // Write orient P
            of << "graphics top color 1" << '\n';
            of << "graphics top cylinder {" << v1.m_x << " " << v1.m_y << " " << v1.m_z << "} {" <<
            v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} radius " <<
            m_vectRadius << " resolution " << m_vectResol << " filled yes" << '\n';

            of << "graphics top color 1" << '\n';
            of << "graphics top cone {" <<
            v1.m_x + v3.m_x*m_vectLength << " " << v1.m_y + v3.m_y*m_vectLength << " " << v1.m_z + v3.m_z*m_vectLength << "} {" <<
            v1.m_x + v3.m_x*m_vectLength*1.20 << " " << v1.m_y + v3.m_y*m_vectLength*1.20 << " " << v1.m_z + v3.m_z*m_vectLength*1.20 << "} radius " <<
            2*m_vectRadius << " resolution " << m_vectResol << '\n';
        }*/

        of.close();

    }
    else
    {
        std::cerr << "Can't open script file for writing rigid fragment orientations" << '\n';
    }
}

void Graph::importTransferIntegral(std::vector<std::string> fileNameTransfer)
{
    std::cout << "Importing rates... " ;

    /// Prepare arrays ///
    //int n(m_site.size()*8);
    int nSite(m_site.size());
    int n(nSite);


    if(m_transferMat != 0)
    {
        int n(m_site.size());
        for (int i = 0; i < n; i++)
        {
            delete[] m_transferMat[i];
        }
        delete[] m_transferMat;
        m_transferMat = 0;

        m_maxTransfer.resize(0);
        m_minTransfer.resize(0);
    }

    m_transferMat = new std::vector<double>*[nSite];
    for (int i = 0; i < nSite; i++)
    {
        m_transferMat[i] = new std::vector<double>[nSite];

        /*for (int j(0) ; j < nSite ; j++)
        {
            m_transferMat[i][j] = 0;
        }*/
    }



    /// Read transfer integrals ///
    for (int k(0) ; k < fileNameTransfer.size() ; k++)
    {
        double maxTransfer(-1e15);
        double minTransfer(1e15);
        std::ifstream iflux(fileNameTransfer[k].c_str());
        if (iflux)
        {
            for (int i(0) ; i < n ; i ++)
            {
                for (int j(0) ; j < n ; j ++)
                {
                    double transferMat;
                    iflux >> transferMat;

                    if(i < nSite && j < nSite)
                    {
                        (m_transferMat[i][j]).push_back(fabs(transferMat));
                        if (fabs(transferMat) > maxTransfer) {maxTransfer = fabs(transferMat);}
                        if (fabs(transferMat) < minTransfer && fabs(transferMat) != 0) {minTransfer = fabs(transferMat);}

                    }
                }
            }
            iflux.close();

            m_minTransfer.push_back(minTransfer);
            m_maxTransfer.push_back(maxTransfer);
        }
        else
        {
            std::cerr << "Can't open rate file." << '\n';
        }

        // Fill upper diagonal of the matrix
        /*for (int i(0) ; i < n ; i ++)
        {
            for (int j(0) ; j < n ; j ++)
            {
                if( j < i)
                {
                    (m_transferMat[i][j])[(m_transferMat[i][j]).size()-1] = m_transferMat[j][i][(m_transferMat[i][j]).size()-1];

                }
            }
        }*/
    }

    std::cout << "Done" << '\n';

}


void Graph::writeTransferIntegral(std::string fileNameTransferPdb, std::string fileNameTransferPsf, bool periodic, int l1, int l2)
{
    std::cout << "Writing the transfer integral distribution... ";
    Atom a;


    int nbEdges(0);
    std::ofstream outFlux(fileNameTransferPdb.c_str());
    int counterG(0);
    if (outFlux.is_open())
    {
        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_system.m_boxSize.m_x;
        outFlux << std::setw(9) << m_system.m_boxSize.m_y;
        outFlux << std::setw(9) << m_system.m_boxSize.m_z;
        outFlux << m_system.m_headerLine;
        outFlux << '\n';

        int counter(0);
        for (int i(0) ; i < m_site.size()  ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++ )
            {

                // Find index for transfer matrix, and indicate if site idx have to be reversed to get the value of the tranfer integral
                int level1(m_site[i].m_level);
                int level2(m_site[i].m_neighbor[j]->m_level);
                int idxTransfer(m_nbrLevels*level1 + level2);

                // Retrieve transfer integral value
                double transfer(0);
                int idx1(i % (m_site.size()/m_nbrLevels) );
                int idx2(m_site[i].m_neighborIdx[j] % (m_site.size()/m_nbrLevels) );
                transfer = m_transferMat[idx1][idx2][idxTransfer];

                if(level1 == l1 && level2 == l2)
                {
                    a.m_type = C;
                    a.m_pos = m_site[i].m_pos;
                    int k(m_site[i].m_neighborIdx[j]);
                    double t(-10);
                    if(transfer != 0 ) {t = log10(transfer*1000); }
                    else {std::cout << "Zero transfer !" << '\n';}

                    outFlux << "ATOM";
                    outFlux << std::setw(7);
                    outFlux << counter+1;
                    outFlux << std::setw(3);
                    outFlux << idxTransfer;
                    outFlux << std::setw(8);
                    outFlux << "RES";
                    outFlux << std::setw(4);
                    outFlux << "site";
                    outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                    outFlux << a.m_pos.m_x;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_y;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_z;
                    //outFlux << std::setw(6) << std::scientific << std::setprecision(0) << std::setw(6);
                    outFlux << std::fixed << std::setw(6);

                    outFlux << (float) t ;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                    outFlux <<  0.0;
                    outFlux << '\n';
                    counter ++;



                    a.m_type = C;
                    a.m_pos = m_site[i].m_neighbor[j]->m_pos;

                    outFlux << "ATOM";
                    outFlux << std::setw(7);
                    outFlux << counter+1;
                    outFlux << std::setw(3);
                    outFlux << idxTransfer;
                    outFlux << std::setw(8);
                    outFlux << "RES";
                    outFlux << std::setw(4);
                    outFlux << "site";
                    outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                    outFlux << a.m_pos.m_x;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_y;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_z;
                    //outFlux << std::setw(6) << std::scientific << std::setprecision(0) << std::setw(6);
                    outFlux << std::fixed << std::setw(6);

                    outFlux << (float) t ;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                    outFlux << 0.0;
                    outFlux << '\n';
                    counter++;

                    if ( periodic || ( !periodic &&  !m_site[i].m_periodic[j] ) )
                    {
                        nbEdges ++;
                    }

                }
            }

        }
        counterG = counter;

        outFlux.close();

    }

    else
    {
        std::cerr << "Unable to open transfer pdb file" << std::endl;
    }

    std::ofstream outFlux2(fileNameTransferPsf.c_str());
    if (outFlux2.is_open())
    {
        outFlux2 << "PSF" << '\n' << '\n' ;
        outFlux2 << "       1 !NTITLE" << '\n';
        outFlux2 << " REMARKS connectivity of the graph" << '\n' << '\n';
        outFlux2 << std::setw(8);
        outFlux2 << counterG;
        outFlux2 << std::setw(7);
        outFlux2 << "!NATOM" << '\n' ;

        int counter(0);
        for (int i(0) ; i < m_site.size()  ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
            {

                // Find index for transfer matrix, and indicate if site idx have to be reversed to get the value of the tranfer integral
                int level1(m_site[i].m_level);
                int level2(m_site[i].m_neighbor[j]->m_level);
                int idxTransfer(m_nbrLevels*level1 + level2);

                // Retrieve transfer integral value
                double transfer(0);
                int idx1(i % (m_site.size()/m_nbrLevels) );
                int idx2(m_site[i].m_neighborIdx[j] % (m_site.size()/m_nbrLevels) );
                transfer = m_transferMat[idx1][idx2][idxTransfer];

                if(level1 == l1 && level2 == l2)
                {
                    outFlux2 << std::setw(8);
                    outFlux2 << counter+1;
                    outFlux2 << std::setw(4);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(3);
                    outFlux2 << 1;
                    outFlux2 << std::setw(7);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(5);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(5);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
                    outFlux2 << 1.0;
                    outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
                    outFlux2 << 1.0;
                    outFlux2 << std::setw(12);
                    outFlux2 << 0;
                    outFlux2 << '\n';
                    counter++;

                    outFlux2 << std::setw(8);
                    outFlux2 << counter+1;
                    outFlux2 << std::setw(4);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(3);
                    outFlux2 << 1;
                    outFlux2 << std::setw(7);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(5);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(5);
                    outFlux2 << "RES";
                    outFlux2 << std::setw(12) << std::fixed << std::setprecision(6);
                    outFlux2 << 1.0;
                    outFlux2 << std::setw(14) << std::fixed << std::setprecision(4);
                    outFlux2 << 1.0;
                    outFlux2 << std::setw(12);
                    outFlux2 << 0;
                    outFlux2 << '\n';
                    counter++;

                }
            }
        }

        outFlux2 << '\n';
        outFlux2 << std::setw(8);
        outFlux2 << nbEdges;
        outFlux2 << std::setw(14);
        outFlux2 << "!NBOND: bonds" << '\n' ;

        counter = 0;
        int counterReal(0);
        for ( int i(0) ; i < m_site.size(); i ++ )
        {
            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++ )
            {

                // Find index for transfer matrix, and indicate if site idx have to be reversed to get the value of the tranfer integral
                int level1(m_site[i].m_level);
                int level2(m_site[i].m_neighbor[j]->m_level);
                int idxTransfer(m_nbrLevels*level1 + level2);

                // Retrieve transfer integral value
                double transfer(0);
                int idx1(i % (m_site.size()/m_nbrLevels) );
                int idx2(m_site[i].m_neighborIdx[j] % (m_site.size()/m_nbrLevels) );
                transfer = m_transferMat[idx1][idx2][idxTransfer];

                if(level1 == l1 && level2 == l2)
                {

                    if( (!periodic && !m_site[i].m_periodic[j]) ||  periodic )
                    {
                        outFlux2 << std::setw(8);
                        outFlux2 << counter+1;
                        outFlux2 << std::setw(8);
                        outFlux2 << counter+2;

                        counterReal += 2;
                        if ( (counterReal / 2) % 4 == 0) { outFlux2 << '\n'; }
                    }

                        counter ++;
                        counter ++;

                }

            }
        }

        outFlux2.close();
    }
    else
    {
        std::cerr << "Unable to open transfer psf file" << std::endl;
    }

    std::cout << "done" << '\n';

}


void Graph::importEnergies(std::string fileNameEnergies)
{
    std::cout << "Importing energies... " ;

    /// Prepare arrays ///
    //int n(m_site.size()*8);
    int nSite(m_site.size());
    int n(nSite);
    for (int i = 0; i < n; i++)
    {
        m_site[i].m_energyLevels.resize(0);
    }


    /// Read energies (and write them in standard scale, i.e. inverse from gorilla output, so LUMO > HOMO e.g. )///
    std::ifstream iflux(fileNameEnergies.c_str());
    if (iflux)
    {
        for (int i(0) ; i < n ; i ++)
        {
            for (int j(0) ; j < m_nbrLevels2Read ; j ++)
            {
                double energy;
                iflux >> energy;

                /// Solve issue in energy variation between pbc boxes ///
                if(i >= n/8)
                {
                    int idx(i%(n/8));
                    energy = -m_site[idx].m_energyLevels[j];
                }

                m_site[i].m_energyLevels.push_back(-energy);
                if(j == 0 ) {m_site[i].m_energy = -energy;}
            }
        }
        iflux.close();
    }
    else
    {
        std::cerr << "Can't open energy file." << '\n';
    }

    std::cout << "Done" << '\n';
}


void Graph::writeHistogramTransferDist(std::string fileNameHistogram, int nbBinTranfer, int nbBinDist, int level)
{
    std::cout << "Write Transfer integrals histogram... ";

    int nbBinI(nbBinTranfer);
    int nbBinJ(nbBinDist);

    double interI(log(m_maxTransfer[level]) - log(m_minTransfer[level]));
    double binLengthI( interI / (double) nbBinI);
    double binLengthJ(30.0 / (double) nbBinJ);
    std::vector<std::vector<int>> hist;
    hist.resize(nbBinI);
    for (int i(0) ; i < nbBinI ; i++)
    {
        for (int j(0) ; j < nbBinJ ; j++)
        {
                hist[i].push_back(0);
        }
    }

    //std::cout << "max tranfert = " << log(m_maxTransfer) << '\n';
    //std::cout << "min tranfert = " << log(m_minTransfer) << '\n';

    for(int i(0) ; i < m_site.size() ; i++)
    {
        //std::cout << "i " << i << '\n';

        for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
        {

            HopType h;
            h = intra;

            for (int l(0) ; l < m_site[i].m_hopType[j].size() ; l++)
            {
                if(m_site[i].m_hopType[j][l] == inter) {h = inter;}
            }

            if(h == inter )
            {
                int k(m_site[i].m_neighborIdx[j]);

                if( k >= i)
                {

                    bool periodic;

                    if(m_transferMat[i][k][level] != 0)
                    {

                        double transfer(log(m_transferMat[i][k][level]));
                        double d = Vect3::distancePbc(m_site[i].m_pos, m_site[k].m_pos, m_system.m_boxSize , periodic);

                        int idxI = floor( (transfer - log(m_minTransfer[level])) / binLengthI );
                        int idxJ = floor(d / binLengthJ);

                        if(idxI == hist.size()) {idxI --;}
                        if(idxJ >= hist[0].size()) {idxJ =  hist[0].size() - 1;}

                        hist[idxI][idxJ] = hist[idxI][idxJ] + 1;



                    }

                    else
                    {
                        std::cout << "Found a zero transfer integral !  "  << i << "    "  << k << "    " <<  (float) m_transferMat[i][k][level] << "    " <<  (float) m_transferMat[k][i][level] << '\n';

                    }

                }



            }

        }
    }


    std::ofstream of(fileNameHistogram.c_str());

    if(of)
    {
        of << "  ";

        for(int i(0) ; i < hist[0].size() ; i++)
        {
            of << i*binLengthJ << "  ";
        }

        of << '\n';

        for(int i(0) ; i < hist.size() ; i++)
        {
            of << log(m_minTransfer[level]) + i*binLengthI << "  ";

            for (int j(0) ; j < hist[i].size() ; j++ )
            {
                of << hist[i][j] << "  ";

            }

            of << '\n';
        }

        of.close();

    }
    else
    {
        std::cerr << "Can't open transfer histogram file" << '\n';
    }

    std::cout << "Done" << '\n';

}



void Graph::writeRawTranfer(std::string fileNameRawTransfer)
{
    std::cout << "Writing raw transfer...";

    std::ofstream os(fileNameRawTransfer.c_str());

    if(os)
    {

        for (int i(0) ; i < m_site.size()  ; i ++)
        {

            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j ++)
            {


                HopType h = m_site[i].m_hopType[j][0];

                Vect3 displacement = m_site[i].m_neighborDisplacement[j];
                double d(displacement.norm());

                double transfer(0);
                Vect3 neighborPbc;

                // Find index for transfer matrix, and indicate if site idx have to be reversed to get the value of the tranfer integral
                int level1(m_site[i].m_level);
                int level2(m_site[i].m_neighbor[j]->m_level);
                int idxTransfer(m_nbrLevels*level1 + level2);
                // Retrieve transfer integral value
                int idx1(i % (m_site.size()/m_nbrLevels) );
                int idx2(m_site[i].m_neighborIdx[j] % (m_site.size()/m_nbrLevels) );
                transfer = m_transferMat[idx1][idx2][idxTransfer];
                //std::cout << "Transfer normal is :" << transfer << '\n';

                // Correct for index if sites are linked by pbc (case fake pbc in Gorilla)
                neighborPbc = m_site[i].m_neighborPbc[j];
                int nbSiteNoPbc(m_site.size()/m_nbrLevels/8);


                if ( !neighborPbc.isNull() && h == inter )
                {
                    int molId1(idx1%nbSiteNoPbc);
                    int boxId1( (idx1 - molId1) / nbSiteNoPbc );
                    int molId2(idx2%nbSiteNoPbc);
                    int boxId2( (idx2 - molId2) / nbSiteNoPbc );

                    int newBoxId1(newBoxPbc(boxId1,neighborPbc));
                    int newBoxId2(newBoxPbc(boxId2,neighborPbc));

                    int newMolId1(newBoxId1*nbSiteNoPbc + molId1);
                    int newMolId2(newBoxId2*nbSiteNoPbc + molId2);

                    transfer = m_transferMat[newMolId1][newMolId2][idxTransfer];

                    /*if(transfer == 0)
                    {
                        std::cout << "It is 0 and it was : " << m_transferMat[idx1][idx2][idxTransfer] << '\n';
                        std::cout << "index old : " << idx1 << "    " << idx2 << "index new : " << newMolId1 << "   " << newMolId2 <<'\n';
                    }*/


                }

                if(h == intra)
                {
                    transfer = 1;
                }

                os << log10(1000*transfer*transfer) << "  " << i << "   " << m_site[i].m_neighborIdx[j] ;
                os << '\n';
            }
        }

        os.close();
        std::cout << "Done." << '\n';

    }
    else
    {
        std::cerr << "Can't open raw transfer file." << '\n';
    }

}

void Graph::writeRawRate(std::string fileNameRawRate)
{
    std::cout << "Writing raw rates...";

    std::ofstream os(fileNameRawRate.c_str());

    if(os)
    {

        for (int i(0) ; i < m_site.size() ; i ++)
        {
            for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j ++)
            {

                os << m_site[i].m_hopTo[j] << "     " << i << "     " << m_site[i].m_neighborIdx[j] << '\n';

            }

        }

        os.close();
        std::cout << "Done." << '\n';

    }
    else
    {
        std::cerr << "Can't open raw rate file." << '\n';
    }
}

void Graph::writeRawEnergies(std::string fileNameRawEnergies)
{
    std::cout << "Writing raw energies...";

    std::ofstream os(fileNameRawEnergies.c_str());

    if(os)
    {
        for (int i(0) ; i < m_site.size() ; i ++)
        {
	    for(int j(0) ; j < m_nbrLevels2Read ; j++)
	    {
	      os << m_site[i].m_energyLevels[j] << "  ";
	    }
	    os << std::endl;
        }

        os.close();
        std::cout << "Done." << '\n';

    }
    else
    {
        std::cerr << "Can't open raw energies file." << '\n';
    }
}


int Graph::writeGeometry(std::string geometryFile) 
{
    std::cout << "Writing geometry... ";

    std::ofstream osGeometry(geometryFile.c_str());
    int nbSiteNoLevel(m_site.size()/m_nbrLevels);

    if(osGeometry)
    {
        for(int i(0) ; i < m_site.size() ; i++)
        {
            for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
            {
	                      
	      int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );
	      int idxNeighbor( m_site[i].m_neighborIdx[j] );
	      int idOrb2( (idxNeighbor - (idxNeighbor % nbSiteNoLevel)) / nbSiteNoLevel );
	      osGeometry << "Site" << i << " " << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb1 << " ";
	      osGeometry << "Site" << idxNeighbor << " " << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb2 << " ";

	      bool periodic;

	      Vect3 bound = Vect3::displacementPbc(m_site[i].m_pos, m_site[i].m_neighbor[j]->m_pos, m_system.m_boxSize);
	      bound = bound / bound.norm();

	      double angleN(acos(m_site[i].m_orientN*bound)*180.0/pi);
	      //if(angleN > 90.0) {angleN = 180.0 - angleN;}

	      double angleP(acos(m_site[i].m_orientP*bound)*180.0/pi);
	      //if(angleP > 90.0) {angleP= 180.0 - angleP;}

	      double dangleP(acos(m_site[i].m_orientP*m_site[i].m_neighbor[j]->m_orientP)*180.0/pi);
	      //if(dangleP > 90.0) {dangleP= 180.0 - dangleP;}

	      double dangleN(acos(m_site[i].m_orientN*m_site[i].m_neighbor[j]->m_orientN)*180.0/pi);
	      //if(dangleN > 90.0) {dangleN= 180.0 - dangleN;}

	      double d = Vect3::distancePbc(m_site[i].m_pos, m_site[i].m_neighbor[j]->m_pos, m_system.m_boxSize , periodic);
	      double dMin(1e15);
	      double dMinEdgeGlob(0);
	      Vect3 boundLoc;
	      Vect3 orientNLoc1;
	      Vect3 orientNLoc2;
	      Vect3 orientPLoc1;
	      Vect3 orientPLoc2;
	      Vect3 orientTLoc1;
	      Vect3 orientTLoc2;

	      // Calculate the same properties locally (between closest repeat units)
	      int id1(0) ; int id2(0);
	      for (int k(0) ; k < m_site[i].m_fragment.size() ; k++)
	      {
		  for(int l(0) ; l < m_site[i].m_neighbor[j]->m_fragment.size() ; l++)
		  {

		      int id1(k);
		      int id2(l);

		      // Centroid distance
		      double dTemp = Vect3::distancePbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize , periodic);

		      // Edge to Edge distance
		      double dMinEdge(1e15);
		      for(int u(0) ; u < m_site[i].m_fragment[k]->m_atom.size() ; u++)
		      {
			  for(int v(0) ; v < m_site[i].m_neighbor[j]->m_fragment[l]->m_atom.size() ; v++)
			  {
			      if(m_site[i].m_fragment[id1]->m_atom[u].m_type != H && m_site[i].m_neighbor[j]->m_fragment[id2]->m_atom[v].m_type != H)
			      {
				  double dTempEdge = Vect3::distancePbc(m_site[i].m_fragment[id1]->m_atom[u].m_pos, m_site[i].m_neighbor[j]->m_fragment[id2]->m_atom[v].m_pos, m_system.m_boxSize , periodic);
				  if(dTempEdge < dMinEdge)
				  {
				      dMinEdge = dTempEdge;
				  }
			      }
			  }
		      }

		      if( dTemp < dMin)
		      {
			  dMin = dTemp;
			  dMinEdgeGlob = dMinEdge;
			  boundLoc = Vect3::displacementPbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize);
			  boundLoc = boundLoc / boundLoc.norm();
			  orientNLoc1 = m_site[i].m_fragment[id1]->m_orientN;
			  orientNLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientN;
			  orientPLoc1 = m_site[i].m_fragment[id1]->m_orientP;
			  orientPLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientP;
			  orientTLoc1 = m_site[i].m_fragment[id1]->m_orientT;
			  orientTLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientT;
		      }
		  }
	      }
	      
	      double N1(orientNLoc1*boundLoc);
	      double N2(orientNLoc2*boundLoc);
	      double P(orientPLoc1*orientPLoc2);

	      double angleNLoc(acos(orientNLoc1*boundLoc)*180.0/pi);
	      //if(angleNLoc > 90.0) {angleNLoc = 180.0 - angleNLoc;}

	      double anglePLoc(acos(orientPLoc1*boundLoc)*180.0/pi);
	      //if(anglePLoc > 90.0) {anglePLoc = 180.0 - anglePLoc;}

	      double danglePLoc(acos(orientPLoc1*orientPLoc2)*180.0/pi);
	      //if(danglePLoc > 90.0) {danglePLoc = 180.0 - danglePLoc;}

	      double dangleNLoc(acos(orientNLoc1*orientNLoc2)*180.0/pi);
	      //if(dangleNLoc > 90.0) {dangleNLoc = 180.0 - dangleNLoc;}

	      double dangleTLoc(acos(orientTLoc1*orientTLoc2)*180.0/pi);
	      //if(dangleNLoc > 90.0) {dangleNLoc = 180.0 - dangleNLoc;}

	      /*Vect3 orientTLoc1 = orientNLoc1^orientPLoc1;
	      orientTLoc1 = orientTLoc1 / orientTLoc1.norm();
	      Vect3 orientTLoc2 = orientNLoc2^orientPLoc2;
	      orientTLoc2 = orientTLoc2 / orientTLoc2.norm();
	      double twist(acos(orientTLoc1*orientTLoc2)*180.0/pi);
	      if(twist > 90.0) {twist = 180.0 - twist;}
	      double dihedral(180.0 - danglePLoc);*/

	      // New definition following Slater orbitals
	      Vect3 base;
	      base.m_x = (boundLoc.m_x - orientNLoc1.m_x * (orientNLoc1*boundLoc)) / (sin(acos(orientNLoc1*boundLoc)));
	      base.m_y = (boundLoc.m_y - orientNLoc1.m_y * (orientNLoc1*boundLoc)) / (sin(acos(orientNLoc1*boundLoc)));
	      base.m_z = (boundLoc.m_z - orientNLoc1.m_z * (orientNLoc1*boundLoc)) / (sin(acos(orientNLoc1*boundLoc)));
	      Vect3 orientNLoc2Base = orientNLoc1 * (orientNLoc2*orientNLoc1) + base * (orientNLoc2*base);
	      orientNLoc2Base = orientNLoc2Base / orientNLoc2Base.norm();
	      double theta1(pi/2-acos(orientNLoc1*boundLoc));
	      double theta2(pi/2-acos(orientNLoc2Base*boundLoc));

	      Vect3 baseOrt = base^orientNLoc1;
	      baseOrt = baseOrt / baseOrt.norm();
	      Vect3 orientNLoc2BaseOrt = orientNLoc1 * (orientNLoc2*orientNLoc1) + baseOrt * (orientNLoc2*baseOrt);
	      orientNLoc2BaseOrt = orientNLoc2BaseOrt / orientNLoc2BaseOrt.norm();
	      double phi(acos(orientNLoc1*orientNLoc2BaseOrt));


	      //outFlux << d << "   " << dMin << "   " << angleN << "    " << angleP  << "    " << dangleN << "    " << dangleP << "   " << angleNLoc << "    " << anglePLoc  << "    " << dangleNLoc << "    " << danglePLoc <<  "    " << transfer << "     " << i << "    " << m_site[i].m_neighborIdx[j] << "   " << dangleTLoc <<'\n';

	      // Output new def with Slater orbitals
	      //osGeometry << dMin << "	" << dMinEdgeGlob << "	" << theta1 << "	" << theta2 << "	" << phi;
	      // Ouput with simple geometry
	      osGeometry << dMin << " " << dMinEdgeGlob << " " << N1 << " " << N2 << " " << P << " " << m_site[i].m_transferIntegral[j][0];
	      
	      
	      osGeometry << std::endl;
            }
        }
        
        osGeometry.close();
    }
    else
    {
        std::cerr << "Can't open " << geometryFile << " file!" << std::endl;
	return -1;
    }

    std::cout << "Done." << std::endl;
    return 0;
}


void Graph::addMolLevels2Site()
{
    std::cout << "Adding the different molecular orbitals to the site list...";

    /// Add as many sites as levels considered on each molecule ///
    int nbSite(m_site.size());
    for(int i(0) ; i < m_nbrLevels ; i ++)
    {
        for (int j(0) ; j < nbSite ; j++)
        {
            if(i == 0)
            {
                m_site[j].m_level = 0;
            }
            else
            {
                conjSegment site;
                site = m_site[j];
                site.m_level = i;
                m_site.push_back(site);
            }

        }
    }

    /// Ensure that vector rescaling has not affected the pointer values ///
    for(int i(0) ; i < m_site.size() ; i++)
    {
        for(int j(0) ; j < m_site[i].m_neighbor.size() ; j++)
        {
            m_site[i].m_neighbor[j] = &m_site[m_site[i].m_neighborIdx[j]];
        }
    }


    /// Add new site edges in the graph, intra-molecular, inter-molecular and inter-level ///
    int sCount(0);
    int nbSiteOld(m_site.size() / m_nbrLevels);
    for(int i(0) ; i < m_nbrLevels ; i ++)
    {
        for (int j(0) ; j < nbSiteOld ; j++)
        {
            // Set the appropriate energy to this MO
            m_site[sCount].m_energy = m_site[sCount].m_energyLevels[i];

            // Update the neighbors list based on the new levels to take into account
            //Inter and intra-molecular
            int nbNeighbor(m_site[sCount].m_neighbor.size());
            for (int k(0) ; k < nbNeighbor ; k ++)
            {
                int idx(m_site[sCount].m_neighborIdx[k]);
		
		std::vector<doubleType> transferSave = m_site[sCount].m_transferIntegral[k];
		int transferIdx(i*m_nbrLevels2Read);
		doubleType transfer(transferSave[transferIdx]);
		std::vector<doubleType> transferVect(1, transfer);
		m_site[sCount].m_transferIntegral[k] = transferVect; 

                for (int n(1) ; n < m_nbrLevels ; n++)
                {
                    m_site[sCount].m_neighbor.push_back(&m_site[n*nbSiteOld+idx]);
                    m_site[sCount].m_neighborIdx.push_back(n*nbSiteOld+idx);

                    m_site[sCount].m_periodic.push_back(m_site[sCount].m_periodic[k]);
                    std::vector<HopType> hVect = m_site[sCount].m_hopType[k];
                    m_site[sCount].m_hopType.push_back(hVect);
                    std::vector<Bound> bVect = m_site[sCount].m_boundFragment[k];
                    m_site[sCount].m_boundFragment.push_back(bVect);
                    std::vector<int> iVect = m_site[sCount].m_ionBridge[k];
                    m_site[sCount].m_ionBridge.push_back(iVect);

                    m_site[sCount].m_neighborPbc.push_back(m_site[sCount].m_neighborPbc[k]);
                    m_site[sCount].m_neighborDisplacement.push_back(m_site[sCount].m_neighborDisplacement[k]);
		    
		    // Keep the transfer integral corresponding to this level transition 
		    transferIdx = i*m_nbrLevels2Read + n;
		    transfer = transferSave[transferIdx];
		    transferVect[0] = transfer;
		    m_site[sCount].m_transferIntegral.push_back(transferVect);

                }
            }


            // Inter-level
            for (int n(0) ; n < m_nbrLevels ; n++)
            {
                if(n != m_site[sCount].m_level)
                {
                    m_site[sCount].m_neighbor.push_back(&m_site[n*nbSiteOld+j]);
                    m_site[sCount].m_neighborIdx.push_back(n*nbSiteOld+j);

                    m_site[sCount].m_periodic.push_back(false);
                    HopType h = level;
                    std::vector<HopType> hVect;
                    hVect.push_back(h);
                    m_site[sCount].m_hopType.push_back(hVect);

                    Vect3 zeroVect(0,0,0);
                    m_site[sCount].m_neighborPbc.push_back(zeroVect);
                    m_site[sCount].m_neighborDisplacement.push_back(zeroVect);
		    
		    std::vector<doubleType> transfer(1,m_H_level);
		    m_site[sCount].m_transferIntegral.push_back(transfer);

                    // A priori we don't need to add information about m_boundFragment, m_ionBridge for further calculations
                    // m_HopTo and m_hopFrom will be set by the computeRates method
                }

            }



            sCount ++;

        }

    }

    std::cout << "Done." << '\n';

}


void Graph::orderSites(double latticeParam)
{
    std::cout << "Ordering site... ";

    m_orderMeanField.resize(0);

    if(latticeParam > 0)
    {
        double xMin(1e15);
        double xMax(-1e15);
        double yMin(1e15);
        double yMax(-1e15);
        double zMin(1e15);
        double zMax(-1e15);

        for (int i(0) ; i < m_site.size()/8 ; i++)
        {
            Vect3 v(m_site[i].m_pos);
            if(v.m_x < xMin) {xMin = v.m_x;}
            if(v.m_x > xMax) {xMax = v.m_x;}
            if(v.m_y < yMin) {yMin = v.m_y;}
            if(v.m_y > yMax) {yMax = v.m_y;}
            if(v.m_z < zMin) {zMin = v.m_z;}
            if(v.m_z > zMax) {zMax = v.m_z;}
        }

        double dimX(xMax-xMin);
        double dimY(yMax-yMin);
        double dimZ(zMax-zMin);
        double dimMax(std::max(std::max(dimX,dimY),dimZ));

        int nbCase(ceil(dimMax/latticeParam));
        std::vector<int> ****triMat = new std::vector<int>***[nbCase];

        for(int i(0) ; i < nbCase ; i++)
        {
            triMat[i] = new std::vector<int>**[nbCase];

            for (int j(0) ; j < nbCase ; j++)
            {
                triMat[i][j] = new std::vector<int>* [nbCase];

                for (int k(0) ; k < nbCase ; k++)
                {
                    triMat[i][j][k] = new std::vector<int>();
                }
            }
        }


        for (int i(0) ; i < m_site.size()/8 ; i++)
        {
            Vect3 v(m_site[i].m_pos);
            int idx( floor((v.m_x-xMin)/latticeParam) );
            int idy( floor((v.m_y-yMin)/latticeParam) );
            int idz( floor((v.m_z-zMin)/latticeParam) );

            (triMat[idx][idy][idz])->push_back(i);
        }


        for(int i(0) ; i < nbCase ; i++)
        {
            for (int j(0) ; j < nbCase ; j++)
            {
                for (int k(0) ; k < nbCase ; k++)
                {
                    for(int n(0) ; n < (triMat[i][j][k])->size() ; n++)
                    {
                        m_orderMeanField.push_back((*triMat[i][j][k])[n]);
                    }

                    delete triMat[i][j][k];
                }

                delete[] triMat[i][j];
            }

            delete[] triMat[i];
        }
        delete[] triMat;

    }

    if(latticeParam == 0)
    {
        for (int i(0) ; i < m_site.size()/repetition ; i++)
        {
            m_orderMeanField.push_back(i);
        }
    }

    if(latticeParam == -1)
    {
        for (int i(0) ; i < m_site.size()/repetition ; i++)
        {
            m_orderMeanField.push_back(i);
        }

        std::random_shuffle( m_orderMeanField.begin(), m_orderMeanField.end());
    }

    std::cout << "Done." << '\n';


}

void Graph::parseTransfer()
{
    /// Warning !!!! Use before addMolLevels2Site ///
    std::cout << "Parsing transfer integrals... ";
    int nbSiteNoPbc(m_site.size()/8);
    // Fix pbc
    for( int i(0) ; i < m_site.size() ; i++)
    { 	
	int molId1(i%nbSiteNoPbc);
	int boxId1( (i - molId1) / nbSiteNoPbc );
	
	// Set correct energy levels
	if(boxId1 != 0) {m_site[i].m_energyLevels = m_site[molId1].m_energyLevels;}
         
        for (int j(0) ; j < m_site[i].m_neighborIdx.size() ; j ++)
        {
	  int idxNeighbor(m_site[i].m_neighborIdx[j]);
	  int molId2(idxNeighbor%nbSiteNoPbc);
	  int boxId2( (idxNeighbor - molId2) / nbSiteNoPbc );
	  int parseCase(-1);

	  /// Check if both sites are on the same chain; if so, it is possible to assign user-defined transfer ///
	  if(m_site[i].m_fragment[0]->m_molNbr == m_site[idxNeighbor].m_fragment[0]->m_molNbr && m_H_intra >= 0)
	  {
	    std::vector<doubleType> transfer(m_nbrLevels2Read*m_nbrLevels2Read, m_H_intra);
	    m_site[i].m_transferIntegral[j] = transfer;
	    parseCase=0;
	  }
	  
	  else
	  { 
	    // First case: sites both in the original box, but molId2 < molId1
	    if(boxId1 == 0 && boxId2 == 0 && molId2 < molId1)
	    {
	      // Find idx in neighbor that corresponds to current site
	      int idxInNeighbor(0);
	      for(int k(0) ; k < m_site[molId2].m_neighborIdx.size() ; k++)
	      {
		if(m_site[molId2].m_neighborIdx[k] == molId1) {idxInNeighbor = k; break;}
	      }
	      m_site[i].m_transferIntegral[j] = m_site[molId2].m_transferIntegral[idxInNeighbor];
	      parseCase=1;
	    }

	    // Second case : sites in same box, but not the original one (note here that the issue molId2 < molId1 is already corrected with previous case, because of the site ordering)
	    if( (boxId1 == boxId2) && (boxId1 != 0) )
	    {
	      // Find idx in molId1 in original box that corresponds to site molId2 in original box
	      int idxInMolId1(0);
	      for(int k(0) ; k < m_site[molId1].m_neighborIdx.size() ; k++)
	      {
		if(m_site[molId1].m_neighborIdx[k] == molId2) {idxInMolId1 = k; break;}
	      }
	      m_site[i].m_transferIntegral[j] = m_site[molId1].m_transferIntegral[idxInMolId1];
	      parseCase=2;
	    }
  
	    // Third case : sites in different box, and not linked by pbc (i.e. transfer already calculated for one direction)
	    if( boxId1 != boxId2  && !m_site[i].m_periodic[j] )
	    {
	      BoxIdx bidx = newBoxNoPbc(boxId1, boxId2);
	      int newBoxId1(bidx.m_b1);
	      int newBoxId2(bidx.m_b2);
	      bool revert(bidx.m_revert);

	      int newMolId1(newBoxId1*nbSiteNoPbc + molId1);
	      int newMolId2(newBoxId2*nbSiteNoPbc + molId2);

	      if(revert)
	      {
		int newMolTemp(newMolId1);
		newMolId1 = newMolId2;
		newMolId2 = newMolTemp;
	      }
	      
	      // Find idx in site 2 that links to site 1
	      int idx(0);
	      for(int k(0) ; k < m_site[newMolId1].m_neighborIdx.size() ; k++)
	      {
		if(m_site[newMolId1].m_neighborIdx[k] == newMolId2) {idx = k; break;}
	      }
	      m_site[i].m_transferIntegral[j] = m_site[newMolId1].m_transferIntegral[idx];
	      parseCase=3;
	    }
  
	    // Fourth case : sites linked by pbc
	    if(m_site[i].m_periodic[j] == true)
	    {
	      
	      Vect3 neighborPbc = m_site[i].m_neighborPbc[j];
	      int newBoxId1(newBoxPbc(boxId1,neighborPbc));
	      int newBoxId2(newBoxPbc(boxId2,neighborPbc));
	      int newMolId1(newBoxId1*nbSiteNoPbc + molId1);
	      int newMolId2(newBoxId2*nbSiteNoPbc + molId2);
	      
	      BoxIdx bidx = newBoxNoPbc(newBoxId1, newBoxId2);
	      newBoxId1 = bidx.m_b1;
	      newBoxId2 = bidx.m_b2;
	      bool revert(bidx.m_revert);

	      newMolId1 = newBoxId1*nbSiteNoPbc + molId1;
	      newMolId2 = newBoxId2*nbSiteNoPbc + molId2;

	      if(revert)
	      {
		int newMolTemp(newMolId1);
		newMolId1 = newMolId2;
		newMolId2 = newMolTemp;
	      }

	      // Find idx in new site 1 that links to new site 2
	      int idx(0);
	      for(int k(0) ; k < m_site[newMolId1].m_neighborIdx.size() ; k++)
	      {
		if(m_site[newMolId1].m_neighborIdx[k] == newMolId2) {idx = k; break;}
	      }
	      m_site[i].m_transferIntegral[j] = m_site[newMolId1].m_transferIntegral[idx];
	      parseCase=4;
	      
	    }
	  }
	  
	  if( m_site[i].m_transferIntegral[j][0] == 666 ) {std::cout << "BUG " << i << "	" << idxNeighbor << "	" << parseCase << std::endl;}
	  
        }
    }   
    std::cout << "Done!" << std::endl;

}


void Graph::writeSystemClustered(std::string fileNameCoord, std::string folder)
{
    std::cout << "Writing the clustered system... ";
    int counter(1);
    Atom a;

    std::ofstream outFlux(folder+fileNameCoord);
    if (outFlux.is_open())
    {
        outFlux << "CRYST1";
        outFlux << std::setw(9) << m_system.m_boxSize.m_x;
        outFlux << std::setw(9) << m_system.m_boxSize.m_y;
        outFlux << std::setw(9) << m_system.m_boxSize.m_z;
        outFlux << m_system.m_headerLine;
        outFlux << '\n';

        for (int i(0) ; i < m_site.size()/m_nbrLevels ; i ++ )
        {
            for (int j(0) ; j < m_site[i].m_fragment.size() ; j ++ )
            {
                for(int k(0) ; k < m_site[i].m_fragment[j]->m_atom.size() ; k++)
                {

                    a = m_site[i].m_fragment[j]->m_atom[k];

                    outFlux << "ATOM";
                    outFlux << std::setw(7);
                    outFlux << counter;
                    outFlux << std::setw(3);
                    outFlux << a.m_name;
                    outFlux << std::setw(8);
                    outFlux << "RES A";
                    outFlux << std::setw(4);
                    outFlux << i+1;
                    outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;

                    outFlux << a.m_pos.m_x;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_y;
                    outFlux << std::setw(8);
                    outFlux << a.m_pos.m_z;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);

                    //outFlux << (double) id1;
                    outFlux << m_site[i].m_clusterId ;
                    outFlux << std::fixed << std::setprecision(2) << std::setw(6);
                    outFlux << m_site[i].m_fragment[j]->m_fragNbr;
                    outFlux << '\n';

                    counter ++;

                }

            }
        }
    }
    else
    {
        std::cerr << "Can't open coordinate file !" << '\n';
    }

    std::cout << "Done." << '\n';
}


void Graph::writeboundDistribution(std::string boundDistrib)
{
    //m_system.boundDistribution(boundDistrib);
}


void Graph::writeTransferDistanceAngleNN(std::string fileName, int l1, int l2, RateType rate)
{
    std::cout << "Writing (distance, angleN, transfer)... ";

    std::ofstream outFlux(fileName.c_str());

    if(outFlux)
    {
        for(int i(0) ; i < m_site.size() ; i++)
        {
            for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
            {

                HopType h = m_site[i].m_hopType[j][0];

                double transfer(0);
                int level1(-1);
                int level2(-1);
                if(rate == maTransfer)
                {
                    // Find index for transfer matrix, and indicate if site idx have to be reversed to get the value of the tranfer integral
                    level1 = m_site[i].m_level;
                    level2 = m_site[i].m_neighbor[j]->m_level;
                    int idxTransfer(m_nbrLevels*level1 + level2);
                    // Retrieve transfer integral value
                    int idx1(i % (m_site.size()/m_nbrLevels) );
                    int idx2(m_site[i].m_neighborIdx[j] % (m_site.size()/m_nbrLevels) );
                    transfer = m_transferMat[idx1][idx2][idxTransfer];
                    //std::cout << "Transfer normal is :" << transfer << '\n';

                    // Correct for index if sites are linked by pbc (case fake pbc in Gorilla)
                    Vect3 neighborPbc;
                    neighborPbc = m_site[i].m_neighborPbc[j];
                    int nbSiteNoPbc(m_site.size()/m_nbrLevels/repetition);
                    if ( !neighborPbc.isNull() && h == inter )
                    {
                        int molId1(idx1%nbSiteNoPbc);
                        int boxId1( (idx1 - molId1) / nbSiteNoPbc );
                        int molId2(idx2%nbSiteNoPbc);
                        int boxId2( (idx2 - molId2) / nbSiteNoPbc );

                        int newBoxId1(newBoxPbc(boxId1,neighborPbc));
                        int newBoxId2(newBoxPbc(boxId2,neighborPbc));

                        int newMolId1(newBoxId1*nbSiteNoPbc + molId1);
                        int newMolId2(newBoxId2*nbSiteNoPbc + molId2);

                        transfer = m_transferMat[newMolId1][newMolId2][idxTransfer];

                        if(transfer == 0)
                        {
                            std::cout << "It is 0 and it was : " << m_transferMat[idx1][idx2][idxTransfer] << '\n';
                            std::cout << "index old : " << idx1 << "    " << idx2 << "index new : " << newMolId1 << "   " << newMolId2 <<'\n';
                        }

                    }
                }

                if(rate == maPheno)
                {
                    transfer = m_site[i].m_hopTo[j];
                }


                //if(level1 == l1 && level2 == l2)
                //if(idx1 != idx2)
                if( h == inter )
                {
                    bool periodic;

                    Vect3 boundLoc;
                    Vect3 orientNLoc1;
                    Vect3 orientNLoc2;
                    Vect3 orientPLoc1;
                    Vect3 orientPLoc2;
                    Vect3 orientTLoc1;
                    Vect3 orientTLoc2;

                    // Calculate the same properties locally (between closest repeat units)
                    for (int k(0) ; k < m_site[i].m_fragment.size() ; k++)
                    {

                        for(int l(0) ; l < m_site[i].m_neighbor[j]->m_fragment.size() ; l++)
                        {

                            int id1(k);
                            int id2(l);

                            // Centroid distance
                            double dTemp = Vect3::distancePbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize , periodic);

                            // Edge to Edge distance
                            double dMinEdge(1e15);
                            for(int u(0) ; u < m_site[i].m_fragment[k]->m_atom.size() ; u++)
                            {
                                for(int v(0) ; v < m_site[i].m_neighbor[j]->m_fragment[l]->m_atom.size() ; v++)
                                {
                                    if(m_site[i].m_fragment[id1]->m_atom[u].m_type != H && m_site[i].m_neighbor[j]->m_fragment[id2]->m_atom[v].m_type != H)
                                    {
                                        double dTempEdge = Vect3::distancePbc(m_site[i].m_fragment[id1]->m_atom[u].m_pos, m_site[i].m_neighbor[j]->m_fragment[id2]->m_atom[v].m_pos, m_system.m_boxSize , periodic);
                                        if(dTempEdge < dMinEdge)
                                        {
                                            dMinEdge = dTempEdge;
                                        }
                                    }
                                }
                            }

                            boundLoc = Vect3::displacementPbc(m_site[i].m_fragment[id1]->m_center, m_site[i].m_neighbor[j]->m_fragment[id2]->m_center, m_system.m_boxSize);
                            boundLoc = boundLoc / boundLoc.norm();

                            orientNLoc1 = m_site[i].m_fragment[id1]->m_orientN;
                            orientNLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientN;
                            orientPLoc1 = m_site[i].m_fragment[id1]->m_orientP;
                            orientPLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientP;
                            orientTLoc1 = m_site[i].m_fragment[id1]->m_orientT;
                            orientTLoc2 = m_site[i].m_neighbor[j]->m_fragment[id2]->m_orientT;

                            double angleFactor = fabs( (orientNLoc1*boundLoc) * (orientNLoc2*boundLoc) );


                            outFlux << dMinEdge << "  " << dTemp << "  " << angleFactor << " ";
                        }

                    }

                    outFlux << i << "   " << m_site[i].m_neighborIdx[j] << " " << level1 << "  " << level2 << "  " << transfer << '\n';
                }


            }
        }

        outFlux.close();
    }
    else
    {
        std::cerr << "Can't open file to write (distance, angleN, transfer)" << std::endl;
    }

    std::cout << "Done." << std::endl;
}

void Graph::readCarrierDistribution(std::string fileNameDistrib)
{
    std::cout << "Reading carrier distribution... ";

    std::ifstream inFlux(fileNameDistrib.c_str());

    if(inFlux)
    {
        double tmp(0);
        double p(0);

        for(int i(0) ; i < m_site.size() ; i++)
        {
            inFlux >> tmp >> tmp >> tmp >> tmp >> p;
            m_site[i].m_occup = p;
        }
    }
    else
    {
        std::cerr << "Can't read distribution file!" << '\n';
    }

    std::cout << "Done" << '\n';
}


void Graph::generateGaussianCom()
{
    std::cout << "generating Gaussian com files... ";
    int nbSiteNoPbc(m_site.size()/8);

    /// Generate com file for sites only ///
    for(int i(0) ; i < nbSiteNoPbc ; i++)
    {
        // Create folder where the calculation will be done
        std::string folder = std::string("site") + std::to_string(i);
        #ifdef _WIN32 || _WIN64
        _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
        #else
        mkdir(folder.c_str(),0777);
        #endif // linux
	
	// Generate com files for the conjugated units
	std::string fileName = folder+ "/" + folder + ".com";
	std::ofstream outFlux(fileName.c_str());
	if(outFlux)
	{
	  outFlux << "%NProcShared=" << m_nbProcSharedGaussian << std::endl;
	  outFlux << "%Mem=" << m_nbProcSharedGaussian*2 << "GB" << std::endl;
	  outFlux << "#zindo SCF=Tight nosymm punch=mo pop=full iop(3/33=1)" << std::endl << std::endl;
	  outFlux << folder << std::endl << std::endl;
	  outFlux << "0 1" << std::endl;
	  
	  // Write atom coordinates
	  for (int j(0) ; j < m_site[i].m_fragment.size() ; j ++ )
	  {
	      for(int k(0) ; k < m_site[i].m_fragment[j]->m_atom.size() ; k++)
	      {

		  Atom a = m_site[i].m_fragment[j]->m_atom[k];
		  std::string line;
		  if(a.m_type == C) {line = "C";}
		  if(a.m_type == O) {line = "O";}
		  if(a.m_type == H) {line = "H";}
		  if(a.m_type == S) {line = "S";}
		  if(a.m_type == N) {line = "N";}
		  outFlux << line;
		  outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
		  outFlux << a.m_pos.m_x;
		  outFlux << std::setw(8);
		  outFlux << a.m_pos.m_y;
		  outFlux << std::setw(8);
		  outFlux << a.m_pos.m_z;
		  outFlux << std::endl;
	      }
	      
	  }
	  outFlux << std::endl;
	  outFlux.close();
	}
	else 
	{
	  std::cerr << "Can't generate com file " << fileName << std::endl;
	}
    }
    
    /// Additionnal com files for a few sites will be needed ///
    std::vector<int> additionnalCom;
	
    /// Generate com file for site dimers ///
    for(int i(0) ; i < m_site.size() ; i++)
    {
      	int molId1(i%nbSiteNoPbc);
	int boxId1( (i - molId1) / nbSiteNoPbc );
      
	for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
	{
	    int idxNeighbor(m_site[i].m_neighborIdx[j]);
	    int molId2(idxNeighbor%nbSiteNoPbc);
	    int boxId2( (idxNeighbor - molId2) / nbSiteNoPbc );
	    
	    // First check if we actually need to calculate this dimer (due to box replication)
	    bool toCalculate(false);
	    if(boxId1 == 0 && idxNeighbor > i) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 2) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 6) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 5) {toCalculate = true;}
	    if(boxId1 == 3 && boxId2 == 4) {toCalculate = true;}
	    
	    // No calculation if sites are linked by pbc
	    if(m_site[i].m_periodic[j])
	      toCalculate = false;
	    
	    if(toCalculate)
            {
	      // Generate com files for the dimer formed by the two conjugated units
	      // Create folder where the calculation will be done
	      std::string folder = std::string("site") + std::to_string(i) + std::string("_") + std::string("site") + std::to_string(idxNeighbor);
	      #ifdef _WIN32 || _WIN64
	      _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
	      #else
	      mkdir(folder.c_str(),0777);
	      #endif // linux
	      std::string fileName = folder+ "/" + folder + ".com";
	      std::ofstream outFlux(fileName.c_str());
	      if(outFlux)
	      {
		outFlux << "%NProcShared=" << m_nbProcSharedGaussian << std::endl;
		outFlux << "%Mem=" << m_nbProcSharedGaussian*2 << "GB" << std::endl;
		outFlux << "#zindo SCF=Tight nosymm punch=mo pop=full iop(3/33=1)" << std::endl << std::endl;
		outFlux << folder << std::endl << std::endl;
		outFlux << "0 1" << std::endl;
		
		// Write atom coordinates
		for (int m(0) ; m < m_site[i].m_fragment.size() ; m ++ )
		{
		    for(int k(0) ; k < m_site[i].m_fragment[m]->m_atom.size() ; k++)
		    {

			Atom a = m_site[i].m_fragment[m]->m_atom[k];
			std::string line;
			if(a.m_type == C) {line = "C";}
			if(a.m_type == O) {line = "O";}
			if(a.m_type == H) {line = "H";}
			if(a.m_type == S) {line = "S";}
			if(a.m_type == N) {line = "N";}
			outFlux << line;
			outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
			outFlux << a.m_pos.m_x;
			outFlux << std::setw(8);
			outFlux << a.m_pos.m_y;
			outFlux << std::setw(8);
			outFlux << a.m_pos.m_z;
			outFlux << std::endl;
		    }
		    
		}
		for (int m(0) ; m < m_site[idxNeighbor].m_fragment.size() ; m ++ )
		{
		    for(int k(0) ; k < m_site[idxNeighbor].m_fragment[m]->m_atom.size() ; k++)
		    {

			Atom a = m_site[idxNeighbor].m_fragment[m]->m_atom[k];
			std::string line;
			if(a.m_type == C) {line = "C";}
			if(a.m_type == O) {line = "O";}
			if(a.m_type == H) {line = "H";}
			if(a.m_type == S) {line = "S";}
			if(a.m_type == N) {line = "N";}
			outFlux << line;
			outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
			outFlux << a.m_pos.m_x;
			outFlux << std::setw(8);
			outFlux << a.m_pos.m_y;
			outFlux << std::setw(8);
			outFlux << a.m_pos.m_z;
			outFlux << std::endl;
		    }
		    
		}
		outFlux << std::endl;
		outFlux.close();
	      }
	      else 
	      {
		std::cerr << "Can't generate com file " << fileName << std::endl;
	      }
	      
	      
	      
	      
	      // Generate com file for the neighbor conjugated unit if this one is not in the main box (i.e. it is in a replicated box)
	      /*if( idxNeighbor >= m_site.size() / 8)
	      {		
		std::string folder = std::string("site") + std::to_string(idxNeighbor);
		#ifdef _WIN32 || _WIN64
		_mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
		#else
		mkdir(folder.c_str(),0777);
		#endif // linux
		
		// Generate com files for the conjugated units
		std::string fileName = folder+ "/" + folder + ".com";
		std::ofstream outFlux(fileName.c_str());
		if(outFlux)
		{
		  outFlux << "%NProcShared=" << m_nbProcSharedGaussian << std::endl;
		  outFlux << "#zindo SCF=Tight nosymm punch=mo pop=full iop(3/33=1)" << std::endl << std::endl;
		  outFlux << folder << std::endl << std::endl;
		  outFlux << "0 1" << std::endl;
		  
		  // Write atom coordinates
		  for (int m(0) ; m < m_site[idxNeighbor].m_fragment.size() ; m ++ )
		  {
		      for(int k(0) ; k < m_site[idxNeighbor].m_fragment[m]->m_atom.size() ; k++)
		      {

			  Atom a = m_site[idxNeighbor].m_fragment[m]->m_atom[k];
			  std::string line;
			  if(a.m_type == C) {line = "C";}
			  if(a.m_type == O) {line = "O";}
			  if(a.m_type == H) {line = "H";}
			  if(a.m_type == S) {line = "S";}
			  if(a.m_type == N) {line = "N";}
			  outFlux << line;
			  outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
			  outFlux << a.m_pos.m_x;
			  outFlux << std::setw(8);
			  outFlux << a.m_pos.m_y;
			  outFlux << std::setw(8);
			  outFlux << a.m_pos.m_z;
			  outFlux << std::endl;
		      }
		      
		  }
		  outFlux << std::endl;
		  outFlux.close();
		}
		else 
		{
		  std::cerr << "Can't generate com file " << fileName << std::endl;
		}
	      }*/    
	      
	      // Check if we need com file for one of the single conjugated units
	      if(boxId1 != 0)
	      {
		bool isPresent(false);
		for(int p(0) ; p < additionnalCom.size() ; p++)
		{
		  if(i == additionnalCom[p])
		  {
		    isPresent = true;
		    break;
		  }
		}
		if(!isPresent)
		{
		  additionnalCom.push_back(i);
		}
	      }
	      if(boxId2 != 0)
	      {
		bool isPresent(false);
		for(int p(0) ; p < additionnalCom.size() ; p++)
		{
		  if(idxNeighbor == additionnalCom[p])
		  {
		    isPresent = true;
		    break;
		  }
		}
		if(!isPresent)
		{
		  additionnalCom.push_back(idxNeighbor);
		}
	      }
	      
	    }
	}
	
    }
    
    /// Generate com file for additional sites ///
    for(int p(0) ; p < additionnalCom.size() ; p++)
    {
	int i(additionnalCom[p]);
	
        // Create folder where the calculation will be done
        std::string folder = std::string("site") + std::to_string(i);
        #ifdef _WIN32 || _WIN64
        _mkdir((std::string(_getcwd(NULL,0)) + "\\"+folder).c_str());
        #else
        mkdir(folder.c_str(),0777);
        #endif // linux
	
	// Generate com files for the conjugated units
	std::string fileName = folder+ "/" + folder + ".com";
	std::ofstream outFlux(fileName.c_str());
	if(outFlux)
	{
	  outFlux << "%NProcShared=" << m_nbProcSharedGaussian << std::endl;
	  outFlux << "%Mem=" << m_nbProcSharedGaussian*2 << "GB" << std::endl;
	  outFlux << "#zindo SCF=Tight nosymm punch=mo pop=full iop(3/33=1)" << std::endl << std::endl;
	  outFlux << folder << std::endl << std::endl;
	  outFlux << "0 1" << std::endl;
	  
	  // Write atom coordinates
	  for (int j(0) ; j < m_site[i].m_fragment.size() ; j ++ )
	  {
	      for(int k(0) ; k < m_site[i].m_fragment[j]->m_atom.size() ; k++)
	      {

		  Atom a = m_site[i].m_fragment[j]->m_atom[k];
		  std::string line;
		  if(a.m_type == C) {line = "C";}
		  if(a.m_type == O) {line = "O";}
		  if(a.m_type == H) {line = "H";}
		  if(a.m_type == S) {line = "S";}
		  if(a.m_type == N) {line = "N";}
		  outFlux << line;
		  outFlux << std::fixed << std::setprecision(3) << std::setw(12) ;
		  outFlux << a.m_pos.m_x;
		  outFlux << std::setw(8);
		  outFlux << a.m_pos.m_y;
		  outFlux << std::setw(8);
		  outFlux << a.m_pos.m_z;
		  outFlux << std::endl;
	      }
	      
	  }
	  outFlux << std::endl;
	  outFlux.close();
	}
	else 
	{
	  std::cerr << "Can't generate com file " << fileName << std::endl;
	}
    }
    
    std::cout << " Done!" << std::endl;


}


void Graph::getEnergyAndTransfer()
{
  
   std::cout << "Calculating transfer integrals and gathering orbitals energies... " << std::endl;

    /// Prepare arrays ///
    //int n(m_site.size()*8);
    int nbSiteNoPbc(m_site.size()/8);
    int nSite(m_site.size());
    int n(nSite);
    
    /// Init arrays to save transfer integrals ///
    m_maxTransfer.resize(0);
    m_minTransfer.resize(0);
    for(int i(0) ; i < m_nbrLevels2Read*m_nbrLevels2Read ; i++)
    {
      m_maxTransfer.push_back(-1e15);
      m_minTransfer.push_back(1e15);
    }
    
    /// Init arrays to save levels energies ///
    for (int i = 0; i < n; i++)
    {
	m_site[i].m_energyLevels.clear();
        m_site[i].m_energyLevels.resize(m_nbrLevels2Read);
    }
   
    
    
    /// Calculate transfer integral for all links existing in the system ///
    double progress(0.0);
#pragma omp parallel shared(progress)
{
std::vector<doubleType> minTransferLoc;
std::vector<doubleType> maxTransferLoc;
for(int i(0) ; i < m_nbrLevels2Read*m_nbrLevels2Read ; i++)
{
  maxTransferLoc.push_back(-1e15);
  minTransferLoc.push_back(1e15);
}


#pragma omp for schedule(dynamic,1)
   for(int i(0) ; i < m_site.size() ; i++)
    {
	/// Site position ///
      	int molId1(i%nbSiteNoPbc);
	int boxId1( (i - molId1) / nbSiteNoPbc );
	
	/// Couple of temporary paths ///
	std::string folder1;
	std::string fileName1;
	std::string folder2;
	std::string fileName2;
	std::string folderDimer;
	std::string fileNameDimer;      
      
        // Set folder of site 1 and path to .pun file
	//CORRECT
        //folder1 = std::string("site") + std::to_string(molId1);
	folder1 = std::string("site") + std::to_string(i);
	fileName1 = folder1+ "/" + folder1 + ".pun";
	
	if(i < nbSiteNoPbc)
	{
	  std::ifstream outFlux1(fileName1.c_str());
	  int open1(0);
	  if(outFlux1)
	  {
	    outFlux1.close();
	  }
	  else 
	  {
	    std::cerr << "Can't open " << fileName1 << std::endl;
	    assert(open1);
	  }
	}
      
	for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
	{	
	    int idxNeighbor(m_site[i].m_neighborIdx[j]);
	    // Neighbor position
	    int molId2(idxNeighbor%nbSiteNoPbc);
	    int boxId2( (idxNeighbor - molId2) / nbSiteNoPbc );
	    
	    // First check if this dimer has actually been calculated (due to box replication)
	    bool toCalculate(false);
	    if(boxId1 == 0 && idxNeighbor > i) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 2) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 6) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 5) {toCalculate = true;}
	    if(boxId1 == 3 && boxId2 == 4) {toCalculate = true;}
	    
	    // No calculation if sites are linked by pbc
	    if(m_site[i].m_periodic[j])
	      toCalculate = false;

	    
	    if(toCalculate)
            {
	      // Vector of transfer integral
	      std::vector<doubleType> transferIntegrals;
	      
	      // Set folder of site 2 and path to .pun file
	      //CORRECT
	      //folder2 = std::string("site") + std::to_string(molId2);
	      folder2 = std::string("site") + std::to_string(idxNeighbor);
	      fileName2 = folder2+ "/" + folder2 + ".pun";		
		
	      // Set folder of dimer and path to .pun file
	      folderDimer = std::string("site") + std::to_string(i) + std::string("_") + std::string("site") + std::to_string(idxNeighbor);
	      fileNameDimer = folderDimer+ "/" + folderDimer + ".pun";
	      std::ifstream outFluxDimer(fileNameDimer.c_str());
	      int openDimer(0);
	      if(outFluxDimer)
	      {
		outFluxDimer.close();
	      }
	      else 
	      {
		std::cerr << "Can't open " << fileNameDimer << std::endl;
		assert(openDimer);
	      }	  
	      //std::cout << fileName1 << "	" << fileName2 << "	" << fileNameDimer << std::endl;
	      
	      /// Finally launch the transfer integral calculation ///
	      std::vector<std::vector<double> > transferResult;
	      //std::string command = std::string("calc_J -p_1 ") + fileName1 + std::string(" -p_2 ") + fileName2 + std::string(" -p_P ") + fileNameDimer + std::string(" -orb_ty_1 ") + std::string("HOMO")+ std::string(" -orb_ty_2 ") + std::string("HOMO");
	      const char* command[11];
	      std::string arg0("calc_J");
	      std::string arg1("-p_1");
	      std::string arg2("-p_2");
	      std::string arg3("-p_P");
	      std::string arg4("-orb_ty_1");
	      std::string arg5(m_levelType);
	      std::string arg6("-orb_ty_2");
	      std::string arg7(m_levelType);
	      command[0] = arg0.c_str();
	      command[1] = arg1.c_str();
	      command[2] = fileName1.c_str();
	      command[3] = arg2.c_str();
	      command[4] = fileName2.c_str();
	      command[5] = arg3.c_str();
	      command[6] = fileNameDimer.c_str();
	      command[7] = arg4.c_str();
	      command[8] = arg5.c_str();
	      command[9] = arg6.c_str();
	      command[10] = arg7.c_str();
	      
	      // Check if calculation results already exist
	      std::cout << fileNameDimer ;
	      bool calcExist(false);
	      std::string resFile = folderDimer+ "/" + folderDimer + ".res";
	      std::ifstream resFileIn(resFile);
	      calcExist = resFileIn.is_open();
	      
	      if(!calcExist)
	      {
		transferResult = calculateTI(11,command, m_nbrLevels2Read-1);
		std::ofstream resFileOut(resFile);
		if(!resFileOut.is_open())
		{
		  std::cerr << "Issue to write transfer integrals to result file " << resFile << std::endl;
		}
		else
		{
		  for (int k(0) ; k < transferResult.size() ; k++)
		  {
		    for(int l(0) ; l < 5 ; l ++)
		    {
		      resFileOut << transferResult[k][l] << std::endl;
		    }
		  }
		  resFileOut.close();
		}
	      }
	      else
	      {
		  for (int k(0) ; k < m_nbrLevels2Read*m_nbrLevels2Read ; k++)
		  {
		    std::vector<double> resultFromFile;
		    for(int l(0) ; l < 5 ; l ++)
		    {
		      double dTemp;
		      resFileIn >> dTemp;
		      resultFromFile.push_back(dTemp);
		    }
		    transferResult.push_back(resultFromFile);
		  }
		  resFileIn.close();
	      }
	      
	      //std::cout <<arg0<< " "<<arg1<<" "<<arg2<<" "<<arg3<<" "<<arg4<< " "<<arg5<<" "<<arg6<<" "<<arg7<<std::endl;
	      //std::cout << fileName1 << "	" << fileName2 << "	" << fileNameDimer << std::endl;
	      std::cout << " done... " << std::endl;

	      
	      /// Write transfer integral in the graph transferMat ///
	      for (int k(0) ; k < transferResult.size() ; k++)
	      {
		doubleType transferMat(transferResult[k][3]);

		transferIntegrals.push_back(transferMat);
		if(transferMat == 0)
		  std::cout << "Zero transfer: " << i << "	" << idxNeighbor  << std::endl;
		if (fabs(transferMat) > maxTransferLoc[k]) {maxTransferLoc[k] = fabs(transferMat);}
		if (fabs(transferMat) < minTransferLoc[k]) {minTransferLoc[k] = fabs(transferMat);}
		
		int idOrb2(k%m_nbrLevels2Read);
		int idOrb1( (k - idOrb2) / m_nbrLevels2Read );
		m_site[i].m_energyLevels[idOrb1] = transferResult[k][0];
		m_site[idxNeighbor].m_energyLevels[idOrb2] = transferResult[k][1];
		
		//std::cout << i << "	" << idxNeighbor << "	" << idOrb1 << "	" << idOrb2 << "	" << transferMat << std::endl;

	      }
	      
	      m_site[i].m_transferIntegral.push_back(transferIntegrals);
	      
	    }
	    
	    else
	    {
	      // Fill with dummy transfer integrals
	      std::vector<doubleType> dummyTransfer(m_nbrLevels2Read*m_nbrLevels2Read, dummyNumberDouble);
	      m_site[i].m_transferIntegral.push_back(dummyTransfer);
	    }
	}
	
#pragma omp critical
{
	progress = progress + 100.0/((double)m_site.size());
	std::cout << "Progressing... " << progress << "%" << std::endl;
}
	
    }
    
    
#pragma omp critical 
    {
      for(int i(0) ; i < m_maxTransfer.size() ; i++)
      {
        if (minTransferLoc[i] < m_minTransfer[i]) 
	{
            m_minTransfer[i] = minTransferLoc[i];
        }
        if(maxTransferLoc[i] > m_maxTransfer[i])
	{
	    m_maxTransfer[i] = maxTransferLoc[i];
	}
      }
    }    
}
  std::cout << "Done" << std::endl;  
}

void Graph::getMillerAbrahams()
{
  
   std::cout << "Calculating transfer integrals from Miller-Abrahams... " << std::endl;

    /// Prepare arrays ///
    //int n(m_site.size()*8);
    int nbSiteNoPbc(m_site.size()/8);
    int nSite(m_site.size());
    int n(nSite);
    
    /// Init arrays to save transfer integrals ///
    m_maxTransfer.resize(0);
    m_minTransfer.resize(0);
    for(int i(0) ; i < m_nbrLevels2Read*m_nbrLevels2Read ; i++)
    {
      m_maxTransfer.push_back(-1e15);
      m_minTransfer.push_back(1e15);
    }
    
    /// Init arrays to save levels energies ///
    for (int i = 0; i < n; i++)
    {
	m_site[i].m_energyLevels.clear();
        m_site[i].m_energyLevels.resize(m_nbrLevels2Read);
    }
   
    
    
    /// Calculate transfer integral for all links existing in the system ///
    double progress(0.0);
#pragma omp parallel shared(progress)
{
std::vector<doubleType> minTransferLoc;
std::vector<doubleType> maxTransferLoc;
for(int i(0) ; i < m_nbrLevels2Read*m_nbrLevels2Read ; i++)
{
  maxTransferLoc.push_back(-1e15);
  minTransferLoc.push_back(1e15);
}


#pragma omp for schedule(dynamic,1)
   for(int i(0) ; i < m_site.size() ; i++)
    {
	/// Site position ///
      	int molId1(i%nbSiteNoPbc);
	int boxId1( (i - molId1) / nbSiteNoPbc );
      
	for (int j(0) ; j < m_site[i].m_neighbor.size() ; j++ )
	{	
	    int idxNeighbor(m_site[i].m_neighborIdx[j]);
	    // Neighbor position
	    int molId2(idxNeighbor%nbSiteNoPbc);
	    int boxId2( (idxNeighbor - molId2) / nbSiteNoPbc );
	    
	    // First check if this dimer has actually been calculated (due to box replication)
	    bool toCalculate(false);
	    if(boxId1 == 0 && idxNeighbor > i) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 2) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 1 && boxId2 == 6) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 4) {toCalculate = true;}
	    if(boxId1 == 2 && boxId2 == 5) {toCalculate = true;}
	    if(boxId1 == 3 && boxId2 == 4) {toCalculate = true;}
	    
	    // No calculation if sites are linked by pbc
	    if(m_site[i].m_periodic[j])
	      toCalculate = false;
	    
	    if(toCalculate)
            {
	      // Vector of transfer integral
	      std::vector<doubleType> transferIntegrals;
	      
	      /// Finally launch the transfer integral calculation with Miller-Abrahams///
	      std::vector<std::vector<doubleType> > transferResult;
	      doubleType transferMA(0.0);
	      HopType h(m_site[i].m_hopType[j][0]);
	      bool periodic;
	      Vect3 displacement = m_site[i].m_neighborDisplacement[j];
	      doubleType d(displacement.norm());
	      doubleType distanceMin(1e8);
	      for(int k(0) ; k < m_site[i].m_fragment.size() ; k++)
	      {
		  bool periodic(true);
		  int id1;
		  int id2;

		  for(int l(0) ; l < m_site[i].m_neighbor[j]->m_fragment.size() ; l++)
		  {
		      doubleType dd(Vect3::distancePbc(m_site[i].m_fragment[k]->m_center,(m_site[i].m_neighbor[j]->m_fragment[l])->m_center,m_system.m_boxSize,periodic));
		      if( dd < distanceMin)
		      {
			  distanceMin = dd;
			  id1 = k;
			  id2 = l;
		      }
		  }
	      }
	      d = distanceMin;
	      if( h == inter )
	      {
		transferMA = exp( - d / m_loc_inter );
	      }
	      if( h == intra )
	      {
		transferMA = exp( - d / m_loc_intra );
	      }
	      if( h == level )
	      {
		transferMA = 0;
	      }

	      std::vector<doubleType> transferTemp;
	      transferTemp.push_back(0.0); transferTemp.push_back(0.0); transferTemp.push_back(0.0); transferTemp.push_back(transferMA); transferTemp.push_back(0.0);
	      transferResult.push_back(transferTemp);

	      
	      /// Write transfer integral in the graph transferMat ///
	      for (int k(0) ; k < transferResult.size() ; k++)
	      {
		doubleType transferMat(transferResult[k][3]);

		transferIntegrals.push_back(transferMat);
		if (fabs(transferMat) > maxTransferLoc[k]) {maxTransferLoc[k] = fabs(transferMat);}
		if (fabs(transferMat) < minTransferLoc[k]) {minTransferLoc[k] = fabs(transferMat);}
		
		int idOrb2(k%m_nbrLevels2Read);
		int idOrb1( (k - idOrb2) / m_nbrLevels2Read );
		m_site[i].m_energyLevels[idOrb1] = transferResult[k][0];
		m_site[idxNeighbor].m_energyLevels[idOrb2] = transferResult[k][1];	

	      }
	      
	      m_site[i].m_transferIntegral.push_back(transferIntegrals);
	      
	    }
	    
	    else
	    {
	      // Fill with dummy transfer integrals
	      std::vector<doubleType> dummyTransfer(m_nbrLevels2Read*m_nbrLevels2Read, dummyNumberDouble);
	      m_site[i].m_transferIntegral.push_back(dummyTransfer);
	    }
	}
	
#pragma omp critical
{
	progress = progress + 100.0/((double)m_site.size());
	std::cout << "Progressing... " << progress << "%" << std::endl;
}
	
    }
    
    
#pragma omp critical 
    {
      for(int i(0) ; i < m_maxTransfer.size() ; i++)
      {
        if (minTransferLoc[i] < m_minTransfer[i]) 
	{
            m_minTransfer[i] = minTransferLoc[i];
        }
        if(maxTransferLoc[i] > m_maxTransfer[i])
	{
	    m_maxTransfer[i] = maxTransferLoc[i];
	}
      }
    }    
}
  std::cout << "Done" << std::endl;  
}

int Graph::saveGraph(std::string outputGraphFile)
{
  std::ofstream os(outputGraphFile, std::ios::out | std::ios::binary); 
  if(os.is_open())
  {
    int nbSite(m_site.size());
    os.write((char*) &nbSite, sizeof(int));
    
    for(int i(0) ; i < nbSite ; i++)
    {
      // Write energies
      os.write((char*) &m_nbrLevels2Read, sizeof(int));
      for(int j(0) ; j < m_nbrLevels2Read ; j++)
      {
	os.write((char*) &m_site[i].m_energyLevels[j], sizeof(double));
      }
      
      // Write transfer integrals
      int nbLevelsCalculated(m_nbrLevels2Read*m_nbrLevels2Read);
      os.write((char*) &nbLevelsCalculated, sizeof(int));
      int nbNeighborTransfer(m_site[i].m_transferIntegral.size());
      os.write((char*) &nbNeighborTransfer, sizeof(int));
      for(int j(0) ; j < nbNeighborTransfer ; j++)
      {
	for(int k(0) ; k < nbLevelsCalculated ; k++)
	{
	  os.write((char*) &m_site[i].m_transferIntegral[j][k], sizeof(doubleType));
	}
      }
      
    }
    
    os.close();
  }
  else
  {
    std::cerr << "Can't open file " << outputGraphFile << std::endl;
    return -1;
  }
  
  return 0;
}

int Graph::loadGraph(std::string inputGraphFile)
{
  std::ifstream is(inputGraphFile, std::ios::in  | std::ios::binary); 
  if(is.is_open())
  {
    int nbSite(0);
    is.read((char*) &nbSite, sizeof(int));
    if(nbSite != m_site.size() ) { std::cerr << "Incorrect number of sites in the graph file!" << std::endl; is.close(); return -1;}
    
    for(int i(0) ; i < nbSite ; i++)
    {
      // Read energies
      int nbrLevels2Read(0);
      is.read((char*) &nbrLevels2Read, sizeof(int));
      if(i==0) {std::cout << "The graph contains " << nbrLevels2Read << " orbitals per site" << std::endl;}
      m_nbrLevels2Read = nbrLevels2Read;
      for(int j(0) ; j < nbrLevels2Read ; j++)
      {
	double energy(0);
	is.read((char*) &energy, sizeof(double));
	m_site[i].m_energyLevels.push_back(energy);
      }
      
      // Read transfer integrals
      int nbLevelsCalculated(0);
      is.read((char*) &nbLevelsCalculated, sizeof(int));
      if(nbLevelsCalculated != m_nbrLevels2Read*m_nbrLevels2Read) { std::cerr << "Incorrect number of calculated transfer in the graph file!" << std::endl; is.close(); return -1;}
      int nbNeighborTransfer(0);
      is.read((char*) &nbNeighborTransfer, sizeof(int));
      for(int j(0) ; j < nbNeighborTransfer ; j++)
      {
	std::vector<doubleType> transferVect;
	for(int k(0) ; k < nbLevelsCalculated ; k++)
	{
	  doubleType transfer(0);
	  is.read((char*) &transfer, sizeof(doubleType));
	  //std::cout << i << "	" << j << "	" << k << "	" << transfer << std::endl;
	  transferVect.push_back(transfer);
	}
	m_site[i].m_transferIntegral.push_back(transferVect);
      }
      
    }
    
    is.close();
  }
  else
  {
    std::cerr << "Can't open file " << inputGraphFile << std::endl;
    return -1;
  }
  
  return 0;
}


int Graph::writeEnergy(std::string energyFile)
{
    std::cout << "Writing density of states to text file " << energyFile << "... ";
    
    std::ofstream osEnergy(energyFile);
    if(osEnergy.is_open())
    {
      int nbSite(m_site.size() / m_nbrLevels / repetition);
      for(int i(0) ; i < nbSite ; i++)
      {
	for(int j(0) ; j < m_nbrLevels ; j++)
	{
	  int idxSite( j * ( nbSite * repetition ) + i );
	  osEnergy << "Site" << i << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << j << "	" << m_site[idxSite].m_energy << std::endl;
	}
      }
      osEnergy.close();
    }
    else
    {
      std::cerr << "Can't open file " << energyFile << std::endl;
      return -1;
    }
    
    std::cout << "Done!" << std::endl; 
    return 0;
}

int Graph::writeTransfer(std::string transferFile)
{
    std::cout << "Writing transfer integral distribution to text file " << transferFile << "... ";
  
    std::ofstream osTransfer(transferFile);
    if(osTransfer.is_open())
    {
      int nbSite(m_site.size());
      int nbSiteNoLevel(m_site.size()/m_nbrLevels);
      for(int i(0) ; i < nbSite ; i++)
      {
	for(int j(0) ; j < m_site[i].m_transferIntegral.size() ; j++)
	{
	  int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );
	  int idxNeighbor( m_site[i].m_neighborIdx[j] );
	  int idOrb2( (idxNeighbor - (idxNeighbor % nbSiteNoLevel)) / nbSiteNoLevel );
	  osTransfer << "Site" << i << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb1 << "	";
	  osTransfer << "Site" << idxNeighbor << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb2 << "	";
	  osTransfer << m_site[i].m_transferIntegral[j][0] << std::endl;	  
	  
	  /*for(int k(0) ; k < m_site[i].m_transferIntegral[j].size() ; k++)
	  {
	    int idOrb2(k%m_nbrLevels2Read);
	    int idOrb1( (k - idOrb2) / m_nbrLevels2Read );
	    osTransfer << "Site" << i << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb1 << "	";
	    osTransfer << "Site" << m_site[i].m_neighborIdx[j] << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb2 << "	";
	    osTransfer << m_site[i].m_transferIntegral[j][k] << std::endl;
	  }*/
	}
      }
      osTransfer.close();
    }
    else
    {
      std::cerr << "Can't open file " << transferFile << std::endl;
      return -1;
    }
    
    std::cout << "Done!" << std::endl; 
    return 0;
}

int Graph::writeRate(std::string rateFile)
{
    std::cout << "Writing hopping rate distribution to text file " << rateFile << "... ";
  
    std::ofstream osRate(rateFile);
    if(osRate.is_open())
    {
      int nbSite(m_site.size());
      int nbSiteNoLevel(m_site.size()/m_nbrLevels);
      for(int i(0) ; i < nbSite ; i++)
      {
	for(int j(0) ; j < m_site[i].m_hopTo.size() ; j++)
	{
	  int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );
	  int idxNeighbor( m_site[i].m_neighborIdx[j] );
	  int idOrb2( (idxNeighbor - (idxNeighbor % nbSiteNoLevel)) / nbSiteNoLevel );
	  osRate << "Site" << i << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb1 << "	";
	  osRate << "Site" << idxNeighbor << "	" << m_levelType << (m_levelType=="HOMO" ? "-" : "+") << idOrb2 << "	";
	  osRate << m_site[i].m_hopTo[j] << std::endl;	  
	}
      }
      osRate.close();
    }
    else
    {
      std::cerr << "Can't open file " << rateFile << std::endl;
      return -1;
    }
    
    std::cout << "Done!" << std::endl; 
    return 0;
}

int Graph::completeGraph(int seed)
{
    /// Add molecular levels to the graph ///
    orderSites(0);
    addMolLevels2Site();
    
    /// Generate DOS from truncated Gaussian distribution ///
    if(m_dosBroadening >= 0)
    {
	gdm(m_dosBroadening,0, seed, false);
    }
    
    /// NEED TO ENSURE TRANSFER SYMMETRY ///
    averageTransfer();
    
    /// Compute rate so that the percolation curve can be calculated ///
    computeRate(maTransfer);
    
    return 0;
	  
}


int Graph::writeGraphVmd(std::string fileNameGraphVmd, int mode)
{
  // Show periodic links ?
  bool showPeriodic(true);

  double minTransferLog(1e15);
  double maxTransferLog(-1e15);
  int nbSiteNoLevel(m_site.size()/m_nbrLevels);
  for(int i(0) ; i < m_site.size() ; i++)
  {
    int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );
    for( int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++)
    {
      if( (!m_site[i].m_periodic[j] && !showPeriodic) || showPeriodic )
      {
        int idN(m_site[i].m_neighborIdx[j]);
        int idOrb2( (idN - (idN % nbSiteNoLevel)) / nbSiteNoLevel );
        double transferLog(0.0);
        
        if(mode == 0)
        {
          transferLog = (log(fabs(m_site[i].m_transferIntegral[j][0])));
        }
        if(mode == 1)
        {
          transferLog = (log(fabs(m_site[i].m_hopTo[j])));
        }
        if(mode == 2)
        {
          transferLog = m_site[i].m_occup * ( 1 - m_site[idN].m_occup ) * m_site[i].m_hopTo[j] * m_site[i].m_displacementField[j] - m_site[idN].m_occup * ( 1 - m_site[i].m_occup ) * m_site[i].m_hopFrom[j] * m_site[i].m_displacementField[j];
	  //transferLog = m_site[i].m_occup * ( 1 - m_site[idN].m_occup ) * m_site[i].m_hopTo[j] - m_site[idN].m_occup * ( 1 - m_site[i].m_occup ) * m_site[i].m_hopFrom[j];
        }
        
        
        if(transferLog < minTransferLog && transferLog != (log(fabs(m_H_level))))
        {
            minTransferLog = transferLog;
        }        
        if(transferLog > maxTransferLog && transferLog != (log(fabs(m_H_level))))
        {
            maxTransferLog = transferLog;
        }
      }
    }
  }

    std::ofstream of(fileNameGraphVmd.c_str());
    
    //Small vector to separate sites on the same conjugated units (i.e. levels)
    Vect3 sep(0.5,0.5,0.5);
    double siteRadius(0.5);
    
    //Color map for site coloring
    std::vector<std::string> colorSite = {"red", "green", "blue", "yellow", "silver", "cyan", "magenta"};
    
    // Link properties
    double cylinderRadius(0.1);
    
    // Get transfer integral min/max
    
    
    if(of)
    {
        for(int i(0) ; i < m_site.size() ; i++)
        {
        
                int idOrb1( (i - (i % nbSiteNoLevel)) / nbSiteNoLevel );

                Vect3 v1(m_site[i].m_pos);
                
                // Write site position as a sphere
                of << "draw color " << colorSite[idOrb1] << '\n';
                of << "set sphere [draw sphere {" << v1.m_x + idOrb1*sep.m_x << " " << v1.m_y + idOrb1*sep.m_y << " " << v1.m_z + idOrb1*sep.m_z << "} radius " << siteRadius << " resolution " << 12 << "]" << '\n';
                
                // Now write transfer integrals to neighbors as cylinder
                for( int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++)
                {
                
                  double cylinderRadiusTemp(cylinderRadius);
                
                  int idN(m_site[i].m_neighborIdx[j]);
                  int idOrb2( (idN - (idN % nbSiteNoLevel)) / nbSiteNoLevel );
                  Vect3 v2(m_site[idN].m_pos);
                
                  if(  (!m_site[i].m_periodic[j] && !showPeriodic) || showPeriodic )
                  {
                    double transferLog(0.0);
                    if(mode == 0)
                    {
                      transferLog = (log(fabs(m_site[i].m_transferIntegral[j][0])));
                      if(transferLog == (log(fabs(m_H_level))))
                      {
                        cylinderRadiusTemp = 0;
                      }
                    }
                    if(mode == 1)
                    {
                      transferLog = (log(fabs(m_site[i].m_hopTo[j])));
                    }
                    if(mode == 2)
                    {
                      transferLog = m_site[i].m_occup * ( 1 - m_site[idN].m_occup ) * m_site[i].m_hopTo[j] * m_site[i].m_displacementField[j] - m_site[idN].m_occup * ( 1 - m_site[i].m_occup ) * m_site[i].m_hopFrom[j] * m_site[i].m_displacementField[j];
                      //transferLog = m_site[i].m_occup * ( 1 - m_site[idN].m_occup ) * m_site[i].m_hopTo[j] - m_site[idN].m_occup * ( 1 - m_site[i].m_occup ) * m_site[i].m_hopFrom[j];
		      if(transferLog == 0)
                      {
                        cylinderRadiusTemp = 0;
                      }
                      else
		      {
			std::cout << transferLog << std::endl;
		      }
                    }

                  
                    
                    int colorNb = 33 + floor( (transferLog - minTransferLog) / (maxTransferLog - minTransferLog) * (1056.0 - 33.0) ); 
                    if(cylinderRadiusTemp == 0)
                    {
                      colorNb = 0;
                    }
                  
                  
                    of << "graphics top color "  << colorNb << '\n';
                    of << "graphics top cylinder {" << v1.m_x + idOrb1*sep.m_x << " " << v1.m_y + idOrb1*sep.m_y << " " << v1.m_z + idOrb1*sep.m_z << "} {" <<
                    v1.m_x + m_site[i].m_neighborDisplacement[j].m_x + idOrb2*sep.m_x << " " << v1.m_y + m_site[i].m_neighborDisplacement[j].m_y + idOrb2*sep.m_y << " " << v1.m_z + m_site[i].m_neighborDisplacement[j].m_z + idOrb2*sep.m_z << "} radius " << cylinderRadiusTemp << " resolution " << 12 << " filled yes" << '\n';
                  }
                
                }
                


        }


        of.close();
        return 0;

    }
    else
    {
        std::cerr << "Can't open script file for writing graph to vmd" << '\n';
        return -1;
    }
}

void Graph::averageTransfer()
{
  // Create symmetric matrix of transfer integral
  int sizeX(m_site.size());
  int sizeY(m_site.size());
  doubleType **ary = new doubleType*[sizeX];
  for(int i = 0; i < sizeX; ++i) 
  {
      ary[i] = new doubleType[sizeY];
      for(int j(0) ; j < sizeY ; j++)
      {
	ary[i][j] = 0.0;
      }
  }

  // Fill in matrix
  for(int i(0) ; i < sizeX ; i++)
  {
    for(int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++)
    {
        int idN(m_site[i].m_neighborIdx[j]);
	ary[i][idN] = ary[i][idN] + m_site[i].m_transferIntegral[j][0];
    }
  }

  // Average transfer and put back in sites
  for(int i(0) ; i < sizeX ; i++)
  {
    for(int j(0) ; j < m_site[i].m_neighborIdx.size() ; j++)
    {
        int idN(m_site[i].m_neighborIdx[j]);
	m_site[i].m_transferIntegral[j][0] = (ary[i][idN] + ary[idN][i]) / 2.0;
    }
  }

  // Delete transfer matrix
  for(int i = 0; i < sizeX; ++i) 
  {
      delete [] ary[i];
  }
  delete [] ary;  
 
}
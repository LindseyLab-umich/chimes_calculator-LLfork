/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

#include<iostream>
#include<vector>
#include<string>
#include<cmath>
#include<algorithm>
#include<fstream>

using namespace std;

#include "serial_jit_interface.h"

// Simple linear algebra functions (for arbitrary triclinic cell support)

inline double mag_a    (const vector<double> & a)
{
    double mag = 0;
    
    for(int i=0; i<a.size(); i++)
        mag += a[i]*a[i];
    
    mag = sqrt(mag);
    
    return mag;
}

inline double a_dot_b  (const vector<double> & a, const vector<double> & b)
{
    double dot = 0;
    
    if (a.size() != b.size())
    {
        cout << "ERROR in a_dot_b: Vectors of different length!" << endl;
        exit(0);
    }
    
    for(int i=0; i<a.size(); i++)
        dot += a[i]*b[i];
    
    return dot;
}

inline void   a_cross_b(const vector<double> & a, const vector<double> & b, vector<double> & cross)
{
    if( a.size() != b.size())
    {
        cout << "ERROR in a_cross_b: Vectors of different length!" << endl;
        exit(0);
    }
    if(a.size() != 3)
    {
        cout << "ERROR in a_cross_b: Vectors should be of length 3!" << endl;
        exit(0);
    }
    
    cross.resize(3);
    cross[0] =    (a[1]*b[2] - a[2]*b[1]);
    cross[1] = -1*(a[0]*b[2] - a[2]*b[0]);
    cross[2] =    (a[0]*b[1] - a[1]*b[0]);
    return;
}    

void set_hmat(const vector<double> & cell_a,const vector<double> & cell_b, const vector<double> & cell_c, vector<double> & hmat, vector<double> & invr_hmat, vector<int> replicates)
{
    // Define the h-matrix (stores the cell vectors locally)

    hmat[0] = cell_a[0]*(replicates[0]+1);  hmat[3] = cell_a[1]*(replicates[0]+1); hmat[6] = cell_a[2]*(replicates[0]+1);
    hmat[1] = cell_b[0]*(replicates[1]+1);  hmat[4] = cell_b[1]*(replicates[1]+1); hmat[7] = cell_b[2]*(replicates[1]+1);
    hmat[2] = cell_c[0]*(replicates[2]+1);  hmat[5] = cell_c[1]*(replicates[2]+1); hmat[8] = cell_c[2]*(replicates[2]+1);

    // Determine the h-matrix inverse
    
    double hmat_det = hmat[0] * (hmat[4]*hmat[8] - hmat[5]*hmat[7])
                    - hmat[1] * (hmat[3]*hmat[8] - hmat[5]*hmat[6])
                    + hmat[2] * (hmat[3]*hmat[7] - hmat[4]*hmat[6]);

    vector<double> tmp_vec(9);
    
    tmp_vec[0] =      (hmat[4]*hmat[8] - hmat[5]*hmat[7]); tmp_vec[3] = -1 * (hmat[1]*hmat[8] - hmat[2]*hmat[7]); tmp_vec[6] =      (hmat[1]*hmat[5] - hmat[2]*hmat[4]);
    tmp_vec[1] = -1 * (hmat[3]*hmat[8] - hmat[5]*hmat[6]); tmp_vec[4] =      (hmat[0]*hmat[8] - hmat[2]*hmat[6]); tmp_vec[7] = -1 * (hmat[0]*hmat[5] - hmat[2]*hmat[3]);
    tmp_vec[2] =      (hmat[3]*hmat[7] - hmat[4]*hmat[6]); tmp_vec[5] = -1 * (hmat[0]*hmat[7] - hmat[1]*hmat[6]); tmp_vec[8] =      (hmat[0]*hmat[4] - hmat[1]*hmat[3]);
    
    invr_hmat[0] = tmp_vec[0]; invr_hmat[3] = tmp_vec[1]; invr_hmat[6] = tmp_vec[2];
    invr_hmat[1] = tmp_vec[3]; invr_hmat[4] = tmp_vec[4]; invr_hmat[7] = tmp_vec[5];
    invr_hmat[2] = tmp_vec[6]; invr_hmat[5] = tmp_vec[7]; invr_hmat[8] = tmp_vec[8];

    invr_hmat[0] /= hmat_det; invr_hmat[3] /= hmat_det; invr_hmat[6] /= hmat_det;
    invr_hmat[1] /= hmat_det; invr_hmat[4] /= hmat_det; invr_hmat[7] /= hmat_det;
    invr_hmat[2] /= hmat_det; invr_hmat[5] /= hmat_det; invr_hmat[8] /= hmat_det;
    
}
    
    
// simulation_system member functions

simulation_system::simulation_system()
{
    hmat     .resize(9);
    invr_hmat.resize(9);
    n_layers.resize(3) ;
    n_replicates.resize(3) ;
    face_dist.resize(3) ;
    max_replicates = 0 ;
}

simulation_system::~simulation_system()
{}

void simulation_system::set_atomtyp_indices(vector<string> & type_list)
{
    sys_atmtyp_indices.resize(n_ghost);

    for(int i=0; i<n_ghost; i++)
    {
        sys_atmtyp_indices[i] = -1;
        
        for (int j=0; j<type_list.size(); j++)
        {
            if ( sys_atmtyps[i] == type_list[j])
            {
                sys_atmtyp_indices[i] = j;
                break;
            }
        }
        
        if (sys_atmtyp_indices[i] == -1)
        {
            cout << "ERROR: Couldn't assign an atom type index for (index/type) " << sys_parent[i] << " " << sys_atmtyps[i] << endl;
            exit(0);
        }
    }
}

void simulation_system::init(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, double max_2b_cut, bool small)
{
    allow_replication = small;
    max_cut = max_2b_cut;
    static bool called_before = false;
    
    //////////////////////////////////////////
    // STEP 1: Copy the system
    //////////////////////////////////////////
    
    n_atoms = x_in.size();
    
    // Sanity checks
    
    if (n_atoms != y_in.size())
    {
        cout << "ERROR: x and y coordinate vector lengths do not match!" << endl;
        exit(0);
    }
    if (n_atoms != z_in.size())
    {
        cout << "ERROR: x and z coordinate vector lengths do not match!" << endl;
        exit(0);
    }
    
    // Copy over the system
    
    n_ghost = n_atoms;
    n_non_repl = n_atoms ;
    
    sys_atmtyp_indices.resize(0);
    
    sys_x.resize(0);
    sys_y.resize(0);
    sys_z.resize(0);
    
    for (int a=0; a<n_atoms; a++)
    {
        sys_atmtyps.push_back(atmtyps[a]);    

        sys_x.push_back( x_in[a] );
        sys_y.push_back( y_in[a] );
        sys_z.push_back( z_in[a] );
        
        sys_parent.push_back(a); // for ghost
    }  
    
    // Determine if system is large enough  

    std::fill(n_replicates.begin(), n_replicates.end(), 0) ;
    std::fill(n_layers.begin(), n_layers.end(), 0) ;

    cell_face_distances(face_dist, n_layers, cella_in, cellb_in, cellc_in) ;
    
    max_replicates = 0 ;
    if (allow_replication)
    {
        for ( int i = 0 ; i < 3 ; i++ )
        {
            n_replicates[i] = ceil(max_cut/face_dist[i])-1;
            if ( n_replicates[i] > max_replicates ) 
                max_replicates = n_replicates[i] ;
        }
    }


    if ( max_replicates > 0)
    {
            if (!called_before)
            {
                called_before = true;

                // Avoid printing the number of replicates on each MD step (LEF).
                cout << "SerialchimesJIT: " << "Replicating the system prior to generating ghost atoms" << endl;
                cout << "SerialchimesJIT: " << "Maximum interaction cutoff = " << max_cut << endl ;
                cout << "SerialchimesJIT: " << "Number of replicates in a,b,c directions = " << n_replicates[0]
                     << " " << n_replicates[1] << " " << n_replicates[2] << endl ;
    
                cout << "SerialchimesJIT: " << "\t" << "Warning: At least one cell length is smaller than the ChIMES outer cutoff." << endl;
                cout << "SerialchimesJIT: " << "\t" << "System will be replicated prior to ghost atom generation." << endl;
                cout << "SerialchimesJIT: " << "\t" << "Results will only be correct for perfectly crystalline cells." << endl;
                cout << "SerialchimesJIT: " << "\t" << "For any other case, system size should be increased." << endl;
            }
    }

    set_hmat(cella_in, cellb_in, cellc_in, hmat, invr_hmat, {0,0,0} );
    
    // Build the replicates - required for correct counting of self-interactions.

    double tmp_x, tmp_y, tmp_z;

    int n_repl  = n_atoms;
    for (int i=0; i <= n_replicates[0] ; i++) // cella_in
    {
        for (int j=0; j <= n_replicates[1] ; j++) // cellb_in
        {
            for (int k=0; k <= n_replicates[2] ; k++) // cellc_in
            {                
                if ( (i==0) && (j==0) && (k==0) )
                   continue;
                
                for (int a=0; a<n_atoms; a++)
                {
                    n_ghost++;
                    n_repl++;    
                    
                    atmtyps.push_back(atmtyps[a]);
                    sys_atmtyps.push_back(atmtyps[a]);    
                    
                    sys_x.push_back(0.0); // Holder    
                    sys_y.push_back(0.0);
                    sys_z.push_back(0.0);
                    
                    
                    
                    tmp_x = sys_x[a] ;
                    tmp_y = sys_y[a] ;
                    tmp_z = sys_z[a] ;

                    // Translate by multiples of cell vectors.
                    
                    tmp_x += i * cella_in[0] + j * cellb_in[0] + k * cellc_in[0] ;
                    tmp_y += i * cella_in[1] + j * cellb_in[1] + k * cellc_in[1] ;
                    tmp_z += i * cella_in[2] + j * cellb_in[2] + k * cellc_in[2] ;                                        
                    
                    sys_x[n_ghost-1] = tmp_x ;
                    sys_y[n_ghost-1] = tmp_y ;
                    sys_z[n_ghost-1] = tmp_z ;
                    
                    sys_parent.push_back(n_repl-1); // As far as ghosts are concerned, these are real atoms
                }
            }
        }
    }

    // Note on replicates: We're essentially tricking the code into thinking the system is bigger than it is, 
    // before any ghost atoms are built 

    n_atoms = n_repl;

    set_hmat(cella_in, cellb_in, cellc_in, hmat, invr_hmat, n_replicates);

    //////////////////////////////////////////
    // STEP 2: Wrap atoms
    //////////////////////////////////////////
        
    // Wrap atoms by converting to inverse space (orthorhombic, cell lengths = unity)
    // Leave in fractional coordinates since we'll need to transform to the new basis later
        
    double tmp_ax, tmp_ay, tmp_az;

    for(int i=0; i<n_atoms; i++)
    {
        // Convert to fractional coordinates
        
        tmp_ax = invr_hmat[0]*sys_x[i] + invr_hmat[1]*sys_y[i] + invr_hmat[2]*sys_z[i];
        tmp_ay = invr_hmat[3]*sys_x[i] + invr_hmat[4]*sys_y[i] + invr_hmat[5]*sys_z[i];
        tmp_az = invr_hmat[6]*sys_x[i] + invr_hmat[7]*sys_y[i] + invr_hmat[8]*sys_z[i];

        // Wrap

        tmp_ax -= floor(tmp_ax);
        tmp_ay -= floor(tmp_ay);
        tmp_az -= floor(tmp_az);
        
        // Revert to absolute coordinates
        
        sys_x[i] = hmat[0]*tmp_ax + hmat[1]*tmp_ay + hmat[2]*tmp_az;    
        sys_y[i] = hmat[3]*tmp_ax + hmat[4]*tmp_ay + hmat[5]*tmp_az;
        sys_z[i] = hmat[6]*tmp_ax + hmat[7]*tmp_ay + hmat[8]*tmp_az;        
    }
    
    //////////////////////////////////////////
    // STEP 3: Determine cell volume after replication.
    //////////////////////////////////////////    

    vector<double> axb(3) ;

    vector<double> a_new({hmat[0], hmat[3], hmat[6]});
    vector<double> b_new({hmat[1], hmat[4], hmat[7]});
    vector<double> c_new({hmat[2], hmat[5], hmat[8]});

    a_cross_b(a_new, b_new, axb) ;
    vol = fabs( a_dot_b(axb, c_new) ) ;
    
}

void simulation_system::build_layered_system(vector<string> & atmtyps, vector<int> & poly_orders, double max_2b_cut, double max_3b_cut, double max_4b_cut)
{
    

    double eff_length = std::max({max_2b_cut,max_3b_cut,max_4b_cut}) ;    

    vector<double> cell_a({hmat[0], hmat[3], hmat[6]}) ;
    vector<double> cell_b({hmat[1], hmat[4], hmat[7]}) ;
    vector<double> cell_c({hmat[2], hmat[5], hmat[8]}) ;        

    find_cell_nlayers(n_layers, cell_a, cell_b, cell_c, eff_length) ;
    cell_face_distances(face_dist, n_layers, cell_a, cell_b, cell_c) ;

    double image_distance = std::min({face_dist[0], face_dist[1], face_dist[2]}) ;

    cout << "Calculated number of image layers = " << n_layers[0] << " " << n_layers[1] << " " << n_layers[2] << endl ;        
    cout << "Distances between image faces = " << face_dist[0] << " " << face_dist[1] << " " << face_dist[2] << endl ;

    
    if (poly_orders[1] >0)
    {
        if (max_3b_cut> image_distance) 
        {
            cout << "ERROR: Maximum 3b cutoff is greater than half at least one box length." << endl;
            cout << "       Increase requested n_layers." << endl;
            cout << "       Max 3b cutoff:            " << max_3b_cut << endl;
            cout << "       Nearest image distance = " << image_distance << endl ;            
            cout << "       nlayers:                  " << n_layers[0] 
                 << " " << n_layers[1] << " " << n_layers[2] << endl;
            exit(0);
        }
    }
    if (poly_orders[2] >0)
    {
        if (max_4b_cut> image_distance) 
        {
            cout << "ERROR: Maximum 4b cutoff is greater than half at least one box length." << endl;
            cout << "       Increase requested n_layers." << endl;
            cout << "       Nearest image distance = " << image_distance << endl ;                        
            cout << "       Max 4b cutoff:            " << max_4b_cut << endl;
            cout << "       nlayers:                  " << n_layers[0] 
                 << " " << n_layers[1] << " " << n_layers[2] << endl;            
            exit(0);
        }
    }

    // Build the layers

    double tmp_x, tmp_y, tmp_z;

    for (int i=-n_layers[0] ; i<=n_layers[0] ; i++) // cell_a
    {
        for (int j=-n_layers[1] ; j<=n_layers[1] ; j++) // cell_b
        {
            for (int k=-n_layers[2] ; k<=n_layers[2] ; k++) // cell_c
            {                
                if ((i==0)&&(j==0)&&(k==0))
                    continue;
                    
                for (int a=0; a<n_atoms; a++)
                {
                    sys_atmtyps.push_back(atmtyps[a]);    

                    // Copy atom position.
                    sys_x.push_back(sys_x[a]); 
                    sys_y.push_back(sys_y[a]);
                    sys_z.push_back(sys_z[a]);

                    // Increment ghost position by lattice vectors.
                    sys_x[n_ghost] += i * cell_a[0] + j * cell_b[0] + k * cell_c[0] ;
                    sys_y[n_ghost] += i * cell_a[1] + j * cell_b[1] + k * cell_c[1] ;
                    sys_z[n_ghost] += i * cell_a[2] + j * cell_b[2] + k * cell_c[2] ;                                        
                    n_ghost++;    

                    sys_parent.push_back(a);
                }
            }
        }
    }
}

void simulation_system::build_neigh_lists(vector<int> & poly_orders, vector<vector<int> > & neighlist_2b, vector<vector<int> > & neighlist_3b, vector<vector<int> > & neighlist_4b, double max_2b_cut, double max_3b_cut, double max_4b_cut)
{
    // Make the 2b neighbor lists
    
    neighlist_2b.resize(n_ghost);
    
    for (int i = 0; i < n_ghost; i++) 
        neighlist_2b[i].resize(0,0);
    
    neighlist_3b.resize(0);
    neighlist_4b.resize(0);

    // Determine search distances

    double search_dist = max_2b_cut;
    bool no_bins = false ;  // Do not use binned neighbor algorithm ?
     
    if (max_3b_cut > search_dist)
        search_dist = max_3b_cut;
    if (max_4b_cut > search_dist)
        search_dist = max_4b_cut;    
    
    // Prepare bins
    

    // Generate neighbor lists on basis of bins
    if ( no_bins )
    {
        // Simple order N^2 algorithm for testing.
        neighlist_no_bins(neighlist_2b, search_dist) ;

    }
    else
    {
        // Populate bins    
        neighlist_bins(neighlist_2b, search_dist) ;

    }

    if ((poly_orders[1] == 0)&&(poly_orders[2]==0))
        return;    
    
    // Make the 3- and 4-b neighbor lists

    bool valid_3mer;
    bool valid_4mer;
    
    vector<int> tmp_3mer(3);
    vector<int> tmp_4mer(4);
    
    int jj, kk, ll;
    
    double dist;
    double dist_ij, dist_ik, dist_il, dist_jk, dist_jl, dist_kl;
    
    for(int i=0; i<n_atoms; i++)
    {
        for(int j=0; j<neighlist_2b[i].size(); j++) // Neighbors of i
        {
            valid_3mer = true;
            valid_4mer = true;
            
            jj = neighlist_2b[i][j];
            
            dist_ij = get_dist(i,jj);

            if (dist_ij >= max_3b_cut)
                valid_3mer = false;
            if (dist_ij >= max_4b_cut)
                valid_4mer = false;
                
            if (!valid_3mer && !valid_4mer)
                continue;
            
            for(int k=0; k<neighlist_2b[i].size(); k++)
            {                                
                kk = neighlist_2b[i][k];
                
                if (jj == kk)
                    continue;
                if (sys_parent[jj] > sys_parent[kk])
                    continue;
                
                dist_ik = get_dist(i,kk); // Check i/k distance
            
                if (dist_ik >= max_3b_cut)
                    valid_3mer = false;
                if (dist_ik >= max_4b_cut)
                    valid_4mer = false;
                    
                dist_jk = get_dist(jj,kk); // Check j/k distance
        
                if (dist_jk >= max_3b_cut)
                    valid_3mer = false;
                if (dist_jk >= max_4b_cut)
                    valid_4mer = false;                                        
                
                if (!valid_3mer && !valid_4mer)
                {
                    if(dist_ij<max_3b_cut)
                        valid_3mer = true;
                    
                    if(dist_ij<max_4b_cut)
                        valid_4mer = true;
                
                    continue;
                }
                
                // If we're here then we have a valid 3-mer ... add it to the 3b neighbor list    
                
                tmp_3mer[0] = i;
                tmp_3mer[1] = jj;
                tmp_3mer[2] = kk;
                
                if (valid_3mer)
                    neighlist_3b.push_back(tmp_3mer);
                
                // Continue on to 4-body list
                
                if (poly_orders[2] == 0)
                    continue;
                
                if (!valid_4mer)
                {
                    if(dist_ij<max_4b_cut)
                        valid_4mer = true;
                    continue;
                }
                
                for(int l=0; l<neighlist_2b[i].size(); l++) 
                {                                
                    ll = neighlist_2b[i][l];
                
                    if (jj == ll)
                        continue;
                    if (kk == ll)
                        continue;
                    if (sys_parent[jj] > sys_parent[ll])
                        continue;                
                    if (sys_parent[kk] > sys_parent[ll])
                        continue;                

                    if (get_dist(i ,ll) >= max_4b_cut) // Check i/l distance
                        continue;                

                    if (get_dist(jj,ll) >= max_4b_cut) // Check j/l distance
                        continue;    

                    if (get_dist(kk,ll) >= max_4b_cut) // Check k/l distance
                        continue;        
                    
                    // If we're here then we have a valid 4-mer ... add it to the 4b neighbor list    
                
                    tmp_4mer[0] = i;
                    tmp_4mer[1] = jj;
                    tmp_4mer[2] = kk;
                    tmp_4mer[3] = ll;
                
                    if (valid_4mer)
                        neighlist_4b.push_back(tmp_4mer);        
                }                                                                        
            }
        }
    }
    /*
    cout << "2B neighbor list is of length:" << neighlist_2b.size() << endl;
    cout << "3B neighbor list is of length:" << neighlist_3b.size() << endl;
    cout << "4B neighbor list is of length:" << neighlist_4b.size() << endl;
    */

}

void simulation_system::run_checks(const vector<double>& max_cuts, vector<int>&poly_orders)
{
    // Sanity check 1: Are the cell vectors long enough?
    const double eps = 1.0e-12 ;
    
    for(int i=0;i<max_cuts.size(); i++)
    {
        if (
            // The non-ghost simulation cell is replicated so that it is at least 1 cutoff wide in all directions.
            // That non-ghost cell needs to be layered at least once in all directions, so you get 3 * cutoff as the limit.
            3.0 * max_cuts[i] > face_dist[0] * (1.0 + eps) || 
            3.0 * max_cuts[i] > face_dist[1] * (1.0 + eps) ||
            3.0 * max_cuts[i] > face_dist[2] * (1.0 + eps) )
        {
            cout << "ERROR: Layered system is smaller than 3x the model " << i+2 <<"-body maximum outer cutoff." << endl;
            cout << "Please report this error to the developers." << endl;
            cout << "Model maximum cutoff: " << max_cuts[i] << endl;
            cout << "Layered system face distance (a): " << face_dist[0] << endl;
            cout << "Layered system face distance (b): " << face_dist[1] << endl;
            cout << "Layered system face distance (c): " << face_dist[2] << endl;  
            exit(0);
            
        }
    }

    
    // Sanity check 2: Does the system have enough atoms?
    
    int bodiedness = 2;
    if (poly_orders[1] > 0)
        bodiedness++;
    if (poly_orders[2] > 0)
        bodiedness++;
    
    if (bodiedness > sys_x.size())
    {
        cout << "ERROR: Layered system contains too few atoms." << endl;
        cout << "   Model bodiedness:            " << bodiedness << endl;
        cout << "   No. atoms in layered system: " << sys_x.size() << endl;
        exit(0);
    }
}
    
// serial_chimes_interface member functions

serial_jit_interface::serial_jit_interface(bool small)
{
    // For small systems, allow explicit replication prior to ghost atom construction
    // This should ONLY be done for perfectly crystalline systems
    
    allow_replication = small;
    
    // Initialize Pointers, etc for chimes calculator interfacing (2-body only for now)
    // To set up for many body calculations, see the LAMMPS implementation

    
    max_2b_cut = 0.0;
    max_3b_cut = 0.0;
    max_4b_cut = 0.0;
    
}
serial_jit_interface::~serial_jit_interface()
{}

void serial_jit_interface::init_chimesJIT(int rank)
{
    // Initialize the chimesJIT object, read parameters from compiled-in string.
    
    init(rank);
    read_parameters() ;
    set_atomtypes(type_list);
    build_pair_int_trip_map() ; 
    build_pair_int_quad_map() ;
}

void serial_jit_interface::build_neigh_lists(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in)
{
    sys.build_neigh_lists(poly_orders, neighlist_2b, neighlist_3b, neighlist_4b, max_cutoff_2B(true), max_cutoff_3B(true), max_cutoff_4B(true));
}

void serial_jit_interface::calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, vector<string> & atmtyps, double & energy, vector<vector<double> > & force, vector<double> & stress)
{
    // Pointers, etc for chimes calculator interfacing (2-body only for now)
    // To set up for many body calculations, see the LAMMPS implementation
        
    double                     dist;
    vector        <double>     dist_3b(3) ;
    vector        <double>     dist_4b(6);
        
    vector        <double>     dr(3) ;
    vector        <double>     dr_3b(3*3) ;
    vector        <double>     dr_4b(6*3) ;
        
    vector<int>                typ_idxs_2b(2) ;
    vector<int>                typ_idxs_3b(3) ;
    vector<int>                typ_idxs_4b(4) ;

    // Read system, set up lattice constants/hmats

    // Determine the max outer cutoff (MUST be 2-body, based on ChIMES logic)

    // Initialize private members (LEF)
    max_2b_cut = max_cutoff_2B(true) ;
    max_3b_cut = max_cutoff_3B(true) ;
    max_4b_cut = max_cutoff_4B(true) ;

    vector<double> stress_chimes(6,0.0) ; // Switch Chimes to a packed stressed tensor.
    
    sys.init(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in, max_2b_cut, allow_replication);   
    
    sys.build_layered_system(atmtyps,poly_orders, max_2b_cut, max_3b_cut, max_4b_cut);

    sys.set_atomtyp_indices(type_list);
    
    sys.run_checks({max_2b_cut,max_3b_cut,max_4b_cut},poly_orders);

    build_neigh_lists(atmtyps, x_in, y_in, z_in, cella_in, cellb_in, cellc_in);

    
    // Setup vars
    
    int ii, jj, kk, ll;
    
    vector<double> force_4b(4*CHDIM) ;
    vector<double> force_3b(3*CHDIM) ;
    vector<double> force_2b(2*CHDIM) ;

    vector<double> force_all(sys.n_atoms*CHDIM, 0.0) ;  
    
    chimes2BTmp chimes_2btmp(poly_orders[0]) ;
    chimes3BTmp chimes_3btmp(poly_orders[1]) ;
    chimes4BTmp chimes_4btmp(poly_orders[2]) ;      
    
    ////////////////////////
    // interate over 1- and 2b's 
    ////////////////////////

    for(int i=0; i<sys.n_atoms; i++)
    {
        compute_1B(sys.sys_atmtyp_indices[i], energy);

        for(int j=0; j<neighlist_2b[i].size(); j++) // Neighbors of i
        {
            jj = neighlist_2b[i][j];

            dist = sys.get_dist(i,jj,dr); // Populates dr, which is passed by ref (overloaded)
            
            typ_idxs_2b[0] = sys.sys_atmtyp_indices[i ];
            typ_idxs_2b[1] = sys.sys_atmtyp_indices[jj];
            
            for (int idx=0; idx<2*CHDIM; idx++) {
                force_2b[idx] = 0.0 ;
            }
            
            compute_2B(dist, dr, typ_idxs_2b, force_2b, stress_chimes, energy, chimes_2btmp);

            for (int idx=0; idx<3; idx++)
            {
                force_all[i*CHDIM+idx]               += force_2b[0*CHDIM+idx] ;
                force_all[sys.sys_parent[jj]*CHDIM+idx]  += force_2b[1*CHDIM+idx] ;     
            }
        }
    }
    
    ////////////////////////
    // interate over 3b's 
    ////////////////////////
    
    if (poly_orders[1] > 0 )
    {
        for(int i=0; i<neighlist_3b.size(); i++)
        {
            ii = neighlist_3b[i][0];
            jj = neighlist_3b[i][1];
            kk = neighlist_3b[i][2];
        
            dist_3b[0] = sys.get_dist(ii,jj,&dr_3b[0]); 
            dist_3b[1] = sys.get_dist(ii,kk,&dr_3b[3]); 
            dist_3b[2] = sys.get_dist(jj,kk,&dr_3b[6]); 
        
            typ_idxs_3b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_3b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_3b[2] = sys.sys_atmtyp_indices[kk];
        
            for (int idx=0; idx<3*CHDIM; idx++)
            {
                force_3b[idx] = 0.0 ;               
            }    
        
            compute_3B(dist_3b, dr_3b, typ_idxs_3b, force_3b, stress_chimes, energy, chimes_3btmp);

            for (int idx=0; idx<3; idx++) {
                force_all[sys.sys_parent[ii]*CHDIM+idx] += force_3b[0*CHDIM+idx] ;
                force_all[sys.sys_parent[jj]*CHDIM+idx] += force_3b[1*CHDIM+idx] ;
                force_all[sys.sys_parent[kk]*CHDIM+idx] += force_3b[2*CHDIM+idx] ;
            }
        }
    }

    ////////////////////////
    // interate over 4b's 
    ////////////////////////

    if (poly_orders[2] > 0 )
    {
        for(int i=0; i<neighlist_4b.size(); i++)
        {
            ii = neighlist_4b[i][0];
            jj = neighlist_4b[i][1];
            kk = neighlist_4b[i][2];
            ll = neighlist_4b[i][3];
        
            dist_4b[0] = sys.get_dist(ii,jj,&dr_4b[0*CHDIM]); 
            dist_4b[1] = sys.get_dist(ii,kk,&dr_4b[1*CHDIM]); 
            dist_4b[2] = sys.get_dist(ii,ll,&dr_4b[2*CHDIM]); 
            dist_4b[3] = sys.get_dist(jj,kk,&dr_4b[3*CHDIM]); 
            dist_4b[4] = sys.get_dist(jj,ll,&dr_4b[4*CHDIM]); 
            dist_4b[5] = sys.get_dist(kk,ll,&dr_4b[5*CHDIM]);         

            typ_idxs_4b[0] = sys.sys_atmtyp_indices[ii];
            typ_idxs_4b[1] = sys.sys_atmtyp_indices[jj];
            typ_idxs_4b[2] = sys.sys_atmtyp_indices[kk];
            typ_idxs_4b[3] = sys.sys_atmtyp_indices[ll];        
        
            for (int idx=0; idx<4*CHDIM; idx++)
            {
                force_4b[idx] = 0.0 ;
            }    
        
            compute_4B(dist_4b, dr_4b, typ_idxs_4b, force_4b, stress_chimes, energy, chimes_4btmp);

            for (int idx=0; idx<3; idx++)
            {
                force_all[sys.sys_parent[ii]*CHDIM+idx] += force_4b[0*CHDIM+idx] ;
                force_all[sys.sys_parent[jj]*CHDIM+idx] += force_4b[1*CHDIM+idx] ;
                force_all[sys.sys_parent[kk]*CHDIM+idx] += force_4b[2*CHDIM+idx] ;
                force_all[sys.sys_parent[ll]*CHDIM+idx] += force_4b[3*CHDIM+idx] ;
            }    
        }    
    }

    // Correct for use of replicates, if applicable

    if ( sys.max_replicates > 0 )
    {
        double rep_scale = 1.0 ;
        for ( int i = 0 ; i < 3 ; i++ )
        {
            rep_scale *= sys.n_replicates[i] + 1.0 ;
        }
        energy /= rep_scale ;
    }

    for ( int ii = 0 ; ii < sys.n_non_repl ; ii++)
        for (int idx = 0 ; idx < CHDIM ; idx++)            
            force[ii][idx] = force_all[CHDIM * ii + idx] ;

    
    ////////////////////////
    // Finish pressure calculation
    ////////////////////////

    // Chimes calculates only unique elements of stress tensor.
    stress[0] = stress_chimes[0] ; // xx
    stress[1] = stress_chimes[1] ; // xy
    stress[2] = stress_chimes[2] ; // xz
    stress[3] = stress_chimes[1] ; // yx
    stress[4] = stress_chimes[3] ; // yy
    stress[5] = stress_chimes[4] ; // yz
    stress[6] = stress_chimes[2] ; // zx
    stress[7] = stress_chimes[4] ; // zy
    stress[8] = stress_chimes[5] ; // zz
    
    for (int idx=0; idx<9; idx++)
        stress[idx] /= sys.vol;  
}

double simulation_system::nearest_image(int layer)
// Returns the distance to the nearest image of an atom (does not depend on atom position in cell)
// belonging to a particular ghost particle layer number.
{
    vector<double> cell_a(3), cell_b(3), cell_c(3) ;
    vector<double> displacement(3) ;
    const double eps = 1.0e-12 ;
    
    cell_a[0] = hmat[0] ; cell_a[1] = hmat[3] ; cell_a[2] = hmat[6] ;
    cell_b[0] = hmat[1] ; cell_b[1] = hmat[4] ; cell_b[2] = hmat[7] ;
    cell_c[0] = hmat[2] ; cell_c[1] = hmat[5] ; cell_c[2] = hmat[8] ;

    if ( fabs(a_dot_b(cell_a, cell_b)) < eps &&
         fabs(a_dot_b(cell_a, cell_c)) < eps &&
         fabs(a_dot_b(cell_b, cell_c)) < eps ) {
        // Orthorhombic cell.  Use mag_a to allow for arbitrary orientation of cell vectors.
        return( layer * std::min({mag_a(cell_a), mag_a(cell_b), mag_a(cell_c)} ) ) ;
    }

    // Non-orthorhombic cell 
    double min_d = 1.0e100 ;
    
    for ( int i = -layer ; i <= layer ; i++ ) {
        for ( int j = -layer ; j <= layer ; j++ ) {
            for ( int k = -layer ; k <= layer ; k++ ) {
                // This restricts us to a specific layer, ignoring the interior layers.
                if ( i != layer && j != layer && k != layer )
                    continue ;

                double distance = 0.0 ;
                for ( int l = 0 ; l < 3 ; l++ ) {
                    displacement[l] = i * cell_a[l] + j * cell_b[l] + k * cell_c[l] ;
                    distance += displacement[l] * displacement[l] ;
                }
                if ( distance < min_d ) 
                    min_d = distance ;
            }
        }
    }
    return sqrt(min_d) ;
}

void simulation_system::neighlist_no_bins(vector<vector<int>> &neighlist_2b, double search_dist)
// Simple order N^2 algorithm for neighbor list generation.  Use for testing purposes.
{
            
    for(int ai=0; ai<n_atoms; ai++) {
        for ( int aj = 0 ; aj < n_ghost ; aj++ ) {
            if ( ai != aj && ai <= sys_parent[aj] && get_dist(ai, aj) < search_dist )
                neighlist_2b[ai].push_back(aj) ;
        }
    }
}

void simulation_system::neighlist_bins(vector<vector<int>> &neighlist_2b, double search_dist)
// Fast order N algorithm for neighbor list generation.  This uses a cartesian binning
// grid that does not necessarily align with cell boundaries.  There is no need for
// alignment when using the explicit image convention.
{

    vector<double> maxpos(3, -1.0e100) ;
    vector<double> minpos(3, +1.0e100) ;


    // Determine limits on position of all particles.
    for(int i=0; i<n_ghost; i++)
    {
        if ( sys_x[i] > maxpos[0] ) 
            maxpos[0] = sys_x[i] ;
        if ( sys_y[i] > maxpos[1] ) 
            maxpos[1] = sys_y[i] ;
        if ( sys_z[i] > maxpos[2] )
            maxpos[2] = sys_z[i] ;

        if ( sys_x[i] < minpos[0] ) 
            minpos[0] = sys_x[i] ;
        if ( sys_y[i] < minpos[1] ) 
            minpos[1] = sys_y[i] ;
        if ( sys_z[i] < minpos[2] )
            minpos[2] = sys_z[i] ;
        
    }

    // Determine the appropriate number of bins.
    vector<int> nbins(3) ;

    for ( int i = 0 ; i < 3 ; i++ )
    {
        nbins[i] = ceil((maxpos[i]-minpos[i]) / search_dist) ;
        if ( nbins[i] < 3 )
        {
            cout << "Error: require at least 3 neighbor bins in all directions.\n" ;
            cout << "The number of layers is not correct\n" ;
        }
    }
    
    int total_bins = nbins[0] * nbins[1] * nbins[2] ;
    vector<vector<int> > bin(total_bins);

    vector<int> bin_idx(3) ;

    // Put every atom into a bin.
    for ( int i = 0 ; i < n_ghost ; i++ )
    {
        
        bin_idx[0] = floor( (sys_x[i] - minpos[0]) / search_dist ) ;
        bin_idx[1] = floor( (sys_y[i] - minpos[1]) / search_dist ) ;
        bin_idx[2] = floor( (sys_z[i] - minpos[2]) / search_dist ) ;
    
        for ( int j = 0 ; j < 3 ; j++ )
        {
            if ( bin_idx[j] < 0 ) 
            {
                cout << "ERROR: Negative bin computed = " << j << " " << bin_idx[j] << endl; // Check .xyz box lengths  and atom coords
                exit(0);
            }
            if ( bin_idx[j] >= nbins[j] )
            {
                cout << "ERROR: Bin index overflow" << j << " " << bin_idx[j] << endl; // Check .xyz box lengths  and atom coords
                exit(0);
            }   
        }
        
        // Calculate bin BIN_IDX of the atom.
        int ibin = bin_idx[0] + bin_idx[1] * nbins[0] + bin_idx[2] * nbins[0]* nbins[1];
        
        if ( ibin >= total_bins || ibin < 0 ) 
        {
            cout << "Error: ibin out of range\n";
            cout << ibin << " " << total_bins << endl;
            exit(1);
        }

        // Push the atom into the bin 
        bin[ibin].push_back(i);
        
    }

    // Calculate neighbors of non-ghost atoms only.  Non-ghost atoms can have ghost neighbors.
    for(int ai=0; ai<n_atoms; ai++)
    {
        bin_idx[0] = floor( (sys_x[ai] - minpos[0]) / search_dist ) ;
        bin_idx[1] = floor( (sys_y[ai] - minpos[1]) / search_dist ) ;
        bin_idx[2] = floor( (sys_z[ai] - minpos[2]) / search_dist ) ;        

        for ( int j = 0 ; j < 3 ; j++ )
        {
            if ( bin_idx[j] < 1 && bin_idx[j] >= nbins[j] - 1 ) 
            {
                cout << "Error: Bad bin for non-ghost atom computed" << endl; // Check .xyz box lengths  and atom coords
                cout << bin_idx[0] << " " << bin_idx[1] << " " << bin_idx[2] << endl;
                cout << sys_x[ai] << " " << sys_y[ai] << " " << sys_z[ai] << endl;
                exit(0);
            }
        }
        

        // Loop over relevant bins only to find neighbors, not all atoms.
        int ibin, aj, ajj, ajend;
        
        for (int i=bin_idx[0]-1; i<= bin_idx[0]+1; i++)    
        {
            for (int j=bin_idx[1]-1; j<=bin_idx[1]+1; j++ ) 
            {
                for (int k=bin_idx[2]-1; k<=bin_idx[2]+1; k++ ) 
                {
                    ibin = i + j * nbins[0] + k * nbins[0] * nbins[1] ;

                    if (ibin >= total_bins || ibin < 0 ) 
                    {
                        cout << "Error: atom bin out of range\n";
                        cout << "ibin = " << ibin << endl ;
                        cout << "bin_idx[0] = " << i << "bin_idx[1] = " << j << "bin_idx[2] = " << k << endl;
                        exit(1);
                    }

                    ajend = bin[ibin].size();
                    
                    for (int aj=0; aj<ajend; aj++) 
                    {
                        ajj = bin[ibin][aj];

                        if ( ajj == ai ) 
                            continue;

                        if ( ai <= sys_parent[ajj]) 
                            if (get_dist(ai,ajj) < search_dist )
                                neighlist_2b[ai].push_back(ajj);        
                    }
                }
            }
        }
    }
}

void simulation_system::find_cell_nlayers(vector<int> & layers, const vector<double> &cell_a, const vector<double> &cell_b,
                                          const vector<double> &cell_c, double cutoff_length)
// Find the required number of layers based on the distances between the cell faces.
{
    // layers = 0 for primitive simulation cell.
    for ( int i = 0 ; i < 3 ; i++ ) {
        layers[i] = 0 ;
    }
    vector<double> w(3) ;

    cell_face_distances(w, layers, cell_a, cell_b, cell_c) ;

    // 3 * cutoff < w[j] * (1 + 2 * nlayers[j])
    for ( int j = 0 ; j < 3 ; j++ ) 
        layers[j] = std::max(ceil( (1.5 * cutoff_length / w[j]) - 0.5 ),1.0) ;
}


void simulation_system::cell_face_distances(vector<double> &w, const vector<int> layers,
                                            const vector<double> &a, const vector<double> &b, const vector<double> &c)
// Calculate the perpendicular distances between cell faces, which are stored in w, given that the
// cell vectors a, b, c are replicated by the given number of layers.
{
    vector<double> acell(3) ;
    vector<double> bcell(3) ;
    vector<double> ccell(3) ;

    // Vector cross products.
    vector<double> axb(3) ;
    vector<double> cxa(3) ;
    vector<double> bxc(3) ;
    
    const double eps = 1.0e-12 ;

    // Get simulation cell vectors.
    for ( int j = 0 ; j < 3 ; j++ )
    {
        acell[j] = a[j] * (2.0*layers[0] + 1.0) ;
        bcell[j] = b[j] * (2.0*layers[1] + 1.0) ;
        ccell[j] = c[j] * (2.0*layers[2] + 1.0) ;
    }

    a_cross_b(acell, bcell, axb) ;
    a_cross_b(ccell, acell, cxa) ;
    a_cross_b(bcell, ccell, bxc) ;    

    if ( mag_a(bxc) < eps || mag_a(cxa) < eps || mag_a(bxc) < eps ) {
        cout << "Error: cell vectors were degenerate\n" ;
        exit(0) ;
    }
    
    w[0] = fabs( a_dot_b(acell, bxc) / mag_a(bxc) )  ;
    w[1] = fabs( a_dot_b(bcell, cxa) / mag_a(cxa) ) ;
    w[2] = fabs( a_dot_b(ccell, axb) / mag_a(axb) ) ;
}



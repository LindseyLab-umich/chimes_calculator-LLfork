/* 
    ChIMES Calculator
    Copyright (C) 2020 Rebecca K. Lindsey, Nir Goldman, and Laurence E. Fried
    Contributing Author:  Rebecca K. Lindsey (2020) 
*/

/* ----------------------------------------------------------------------

This class demonstrates how chimesJIT{h,cpp} can be used to obtain the 
stress tensor, energy, and per-atom forces for a given system. See 
main.cpp for a usage example.

Notes: This class has been written for readability rather than speed.
Optimization and parallel distribution is recommended prior to use with 
large systems.

---------------------------------------------------------------------- */

#ifndef _serial_chimes_interface_h
#define _serial_chimes_interface_h

#include<vector>
#include<string>

using namespace std;

#include "chimesJIT.h"    

class simulation_system
{
    public: 
        
        simulation_system();
        ~simulation_system();
        
        inline double get_dist(int i,int j, vector<double> & rij);
        inline double get_dist(int i,int j, double* rij);   
        inline double get_dist(int i,int j);
        
        void init(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in, vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in, double max_2b_cut, bool small = false);
        void set_atomtyp_indices(vector<string> & type_list);
        void copy(simulation_system & to);
        void reorient();
        void build_layered_system(vector<string> & atmtyps, vector<int> & poly_orders, double max_2b_cut, double max_3b_cut, double max_4b_cut);
        void build_neigh_lists(vector<int> & poly_orders, vector<vector<int> > & neighlist_2b, vector<vector<int> > & neighlist_3b, vector<vector<int> > & neighlist_4b, double max_2b_cut, double max_3b_cut, double max_4b_cut);
        void run_checks(const vector<double>& max_cuts, vector<int>&poly_orders);
        double nearest_image(int layer) ;
    
        bool allow_replication;        // If true, replicates coordinates prior to calculation

        vector<int> n_replicates ;     // number of "real" replicate layers to make in the a,b,c directions.
        int max_replicates ;           // Maximum number of replicates.
        vector<int> n_layers;          // number of ghost layers to make
        int n_atoms;                   // number of real atoms
        int n_ghost;                   // number of real+ghost atoms
        int n_non_repl ;               // Number of non-replicated atoms.
        double max_cut;

        vector<double> face_dist ;         // Distance between faces of the simulation cell.
        vector<int>       sys_atmtyp_indices;   // Atom type indices for all (real+ghost) atoms        
        vector<string>    sys_atmtyps;          // Chemical symbols for all (real+ghost) atoms      
    
        vector<double> sys_x;          // System (i.e. ghost+real) x-coordinates
        vector<double> sys_y;          // System (i.e. ghost+real) y-coordinates
        vector<double> sys_z;          // System (i.e. ghost+real) z-coordinates
        vector<int>    sys_parent;     // Index of atom i's (real) parent
       
        double         vol;         // System volume    
        
    private: 

        void neighlist_bins(vector<vector<int>> &neighlist_2b, double search_dist) ;
        void neighlist_no_bins(vector<vector<int>> &neighlist_2b, double search_dist) ;
        void cell_face_distances(vector<double> &w, const vector<int> layers,
                                 const vector<double> &a, const vector<double> &b, const vector<double> &c) ;
        void find_cell_nlayers(vector<int> & layers, const vector<double> &cell_a, const vector<double> &cell_b,
                               const vector<double> &cell_c, double cutoff_length) ;
    
        vector<double>    hmat;        // System h-matrix
        vector<double>    invr_hmat;   // Inverse h-matrix
        
};


class serial_jit_interface : public chimesJIT
{
    public:
            
        serial_jit_interface(bool small = true);
        ~serial_jit_interface();
        
        bool allow_replication; // If true, replicates coordinates prior to calculation
           
        void    init_chimesJIT(int rank);
        void    calculate(vector<double> & x_in, vector<double> & y_in, vector<double> & z_in,
                          vector<double> & cella_in, vector<double> & cellb_in, vector<double> & cellc_in,
                          vector<string> & atmtyps, double & energy, vector<vector<double> > & force, vector<double> & stress);

    private:

        void build_neigh_lists(vector<string> & atmtyps, vector<double> & x_in, vector<double> & y_in,
                               vector<double> & z_in, vector<double> & cella_in, vector<double> & cellb_in,
                               vector<double> & cellc_in);
    
        simulation_system sys;      // Input system
        vector<string>    type_list;   // A list of possible unique atom types and thier corresponding numerical index, per the parameter file
        
        double max_2b_cut;    // Maximum 2-body outer cutoff
        double max_3b_cut;    // Maximum 3-body outer cutoff
        double max_4b_cut;    // Maximum 4-body outer cutoff
        
        vector<vector<int> > neighlist_2b;    // [real atom index][list of real/ghost atom neighbors]
        vector<vector<int> > neighlist_3b;    // [interaction set index][list of 3 atoms within interaction range] -- currently unused
        vector<vector<int> > neighlist_4b;    // [interaction set index][list of 4 atoms within interaction range] -- currently unused    
        
    

};

inline double simulation_system::get_dist(int i,int j, vector<double> & rij)
{
    rij[0] = sys_x[j] - sys_x[i];
    rij[1] = sys_y[j] - sys_y[i];
    rij[2] = sys_z[j] - sys_z[i];

    return sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);    
}


inline double simulation_system::get_dist(int i,int j, double *rij)
{
    /* Removed code using the h matrix.  Distance does not depend on cell vectors types unless the minimum 
       image convention is used */

    rij[0] = sys_x[j] - sys_x[i];
    rij[1] = sys_y[j] - sys_y[i];
    rij[2] = sys_z[j] - sys_z[i];

    return sqrt(rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2]);
}

inline double simulation_system::get_dist(int i,int j)
{
    vector<double> rij(3);
    
    return get_dist(i,j,rij);
}

#endif































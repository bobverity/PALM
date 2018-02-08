
#include <Rcpp.h>
#include <chrono>
#include <algorithm>
#include "misc.h"
#include "probability.h"

using namespace std;

// find path
vector< vector<int> > findPath(vector<double> &distance, vector<int> &visited, vector<int> &traceback, vector<double> &friction, vector< vector<int> > &neighbors, int source, vector<int> &dest, vector<double> &s);

vector<int> findPath(vector<double> &distance, vector<int> &visited, vector<int> &traceback, vector<double> &friction, vector< vector<int> > &neighbors, int source, int dest, double &s);

// propose hexs by searching from focal hex with given depth
void prop_pointSearch(vector<int> &h_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth);

// propose hexs by searching in random path from focal hex
void prop_squiggle(vector<int> &h_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth);

//------------------------------------------------
// Run model
// [[Rcpp::export]]
Rcpp::List fitModel_cpp(Rcpp::List args_data, Rcpp::List args_model) {
    
    // extract data
    vector<int> x = Rcpp_to_vector_int(args_data["hex_data"]);
    vector< vector<int> > hex_neighbors = Rcpp_to_mat_int(args_data["hex_neighbors"]);
    vector< vector<double> > stat = Rcpp_to_mat_double(args_data["stat"]);
    int n = x.size();
    int n2 = 0.5*n*(n-1);
    int nhex = hex_neighbors.size();
    
    // extract model parameters
    int reps = Rcpp_to_int(args_model["reps"]);
    vector<double> friction = Rcpp_to_vector_double(args_model["friction"]);
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // create objects for finding shortest path
    vector<double> distance(nhex);
    vector<int> visited(nhex);
    vector<int> traceback(nhex);
    
    // create vectors for conveniently moving between linear and matrix indexing
    vector<int> index_i(n2);
    vector<int> index_j(n2);
    vector< vector<int> > index_mat(n, vector<int>(n));
    int ind = 0;
    for (int i=0; i<(n-1); i++) {
        for (int j=(i+1); j<n; j++) {
            index_i[ind] = i;
            index_j[ind] = j;
            index_mat[i][j] = ind;
            ind++;
        }
    }
    
    // find initial shortest paths
    vector< vector<int> > path_mat(n2);
    vector<double> s(n2);
    vector<double> y(n2);
    for (int ind=0; ind<n2; ind++) {
        int i = index_i[ind];
        int j = index_j[ind];
        y[ind] = stat[i][j];
        path_mat[ind] = findPath(distance, visited, traceback, friction, hex_neighbors, x[i], x[j], s[ind]);
    }
    
    // store initial path lengths
    vector<double> s0 = s;
    vector<double> friction_best = friction;
    
    // store paths transecting each hex
    vector< vector<int> > hex_path(nhex);
    vector<int> hex_npath(nhex);
    for (int ind=0; ind<n2; ind++) {
        for (int k=0; k<int(path_mat[ind].size()); k++) {
            int thisHex = path_mat[ind][k];
            hex_path[thisHex].push_back(ind);
            hex_npath[thisHex] ++;
        }
    }
    
    // calculate regression coefficients
    double Sx=0, Sy=0, Sxx=0, Syy=0, Sxy=0;
    for (int ind=0; ind<n2; ind++) {
        Sx += s[ind];
        Sy += y[ind];
        Sxx += s[ind]*s[ind];
        Syy += y[ind]*y[ind];
        Sxy += s[ind]*y[ind];
    }
    double m = (n2*Sxy - Sx*Sy)/(n2*Sxx - Sx*Sx);
    double c = (Sy - m*Sx)/double(n2);
    double SS_res = Syy - 2*m*Sxy - 2*c*Sy + m*m*Sxx + 2*m*c*Sx + n2*c*c;
    double SS_res_global = SS_res;
    
    // create Guided Local Search (GLS) objects
    vector<int> GLSpenalty(nhex);
    
    // objects for storing results
    vector<double> SS_res_store(reps);
    
    // begin main model loop
    for (int rep=0; rep<reps; rep++) {
        //print(rep);
        
        // report progress
        if ( ((rep+1) % 1000)==0 ) {
            print("      iteration",rep+1);
        }
        
        // find best improvement to any hex
        int best_hex = -1;
        double best_delta = 0;
        double best_SS_res = SS_res;
        for (int i=0; i<nhex; i++) {
            
            // skip if no paths through this hex
            if (hex_npath[i]==0) {
                continue;
            }
            
            // calculate sum of x and y in this hex
            double Rx = 0, Ry = 0;
            for (int j=0; j<hex_npath[i]; j++) {
                Rx += s[hex_path[i][j]];
                Ry += y[hex_path[i][j]];
            }
            int k = hex_npath[i];
            
            // find optimal delta
            for (int j=0; j<101; j++) {
                
                // update coefficients
                double delta = 0.1*(-50+j);
                double Sx_new = Sx + k*delta;
                double Sxx_new = Sxx + 2*delta*Rx + k*delta*delta;
                double Sxy_new = Sxy + delta*Ry;
                
                double m_new = (n2*Sxy_new - Sx_new*Sy)/(n2*Sxx_new - Sx_new*Sx_new);
                double c_new = (Sy - m_new*Sx_new)/double(n2);
                double SS_res_new = Syy - 2*m_new*Sxy_new - 2*c_new*Sy + m_new*m_new*Sxx_new + 2*m_new*c_new*Sx_new + n2*c_new*c_new;
                
                // replace if this beats current best hex
                if ((SS_res_new + GLSpenalty[i]) < best_SS_res) {
                    best_hex = i;
                    best_delta = delta;
                    best_SS_res = SS_res_new;
                }
            }
        }
        
        // if no improvement then reached local minimum
        if (best_hex<0 || ((rep+1) % 100)==0 ) {
        //if (best_hex<0) {
            //print("local minumum at iteration ", rep+1);
            
            // find highest GLS utility
            //double best_hex = 0;
            double best_utility = 0;
            for (int i=0; i<nhex; i++) {
                /*
                // skip if no paths through this hex or zero friction
                if (hex_npath[i]==0 || friction[i]==0) {
                    continue;
                }
                
                // calculate sum of x and y in this hex
                double Rx = 0, Ry = 0;
                for (int j=0; j<hex_npath[i]; j++) {
                    Rx += s[hex_path[i][j]];
                    Ry += y[hex_path[i][j]];
                }
                int k = hex_npath[i];
                
                // update coefficients
                double delta = -friction[i];
                double Sx_new = Sx + k*delta;
                double Sxx_new = Sxx + 2*delta*Rx + k*delta*delta;
                double Sxy_new = Sxy + delta*Ry;
                
                double m_new = (n2*Sxy_new - Sx_new*Sy)/(n2*Sxx_new - Sx_new*Sx_new);
                double c_new = (Sy - m_new*Sx_new)/double(n2);
                double SS_res_new = Syy - 2*m_new*Sxy_new - 2*c_new*Sy + m_new*m_new*Sxx_new + 2*m_new*c_new*Sx_new + n2*c_new*c_new;
                
                // calculate cost and utility
                double GLScost = 1/(1 + SS_res_new - SS_res);
                double GLSutility = GLScost/double(1+GLSpenalty[i]);
                */
                
                if (friction[i]==0) {
                    continue;
                }
                double GLScost = 10.0/abs(friction[i]);
                double GLSutility = GLScost/double(1+GLSpenalty[i]);
                
                // replace if this beats current best hex
                if (GLSutility > best_utility) {
                    best_hex = i;
                    best_utility = GLSutility;
                }
                
            }
            
            // increment penalty
            GLSpenalty[best_hex] ++;
            friction[best_hex] = 0;
            
        } else {
            
            // implement improvement to best hex
            friction[best_hex] += best_delta;
            
        }
        
        // recalculate paths
        int best_npath = hex_npath[best_hex];
        vector< vector<int> > store_path(best_npath);
        vector<int> ind_vec(best_npath);
        for (int i=0; i<best_npath; i++) {
            int ind = hex_path[best_hex][i];
            ind_vec[i] = ind;
            store_path[i] = findPath(distance, visited, traceback, friction, hex_neighbors, x[index_i[ind]], x[index_j[ind]], s[ind]);
        }
        
        // remove old paths from hexs
        for (int i=0; i<best_npath; i++) {
            int ind = ind_vec[i];
            for (int j=0; j<int(path_mat[ind].size()); j++) {
                int h = path_mat[ind][j];
                hex_path[h].erase(remove(hex_path[h].begin(), hex_path[h].end(), ind), hex_path[h].end());
                hex_npath[h] --;
            }
        }
        
        // add new paths to hexs
        for (int i=0; i<best_npath; i++) {
            int ind = ind_vec[i];
            for (int j=0; j<int(store_path[i].size()); j++) {
                int h = store_path[i][j];
                hex_path[h].push_back(ind);
                hex_npath[h] ++;
            }
        }
        
        // replace path in path_mat
        for (int i=0; i<best_npath; i++) {
            int ind = ind_vec[i];
            path_mat[ind] = store_path[i];
        }
        
        // recalculate regression coefficients
        Sx=0; Sxx=0; Sxy=0;
        for (int ind=0; ind<n2; ind++) {
            Sx += s[ind];
            Sxx += s[ind]*s[ind];
            Sxy += s[ind]*y[ind];
        }
        m = (n2*Sxy - Sx*Sy)/(n2*Sxx - Sx*Sx);
        c = (Sy - m*Sx)/double(n2);
        SS_res = Syy - 2*m*Sxy - 2*c*Sy + m*m*Sxx + 2*m*c*Sx + n2*c*c;
        
        // store fit
        SS_res_store[rep] = SS_res;
        if (SS_res < SS_res_global) {
            SS_res_global = SS_res;
            friction_best = friction;
        }
        
    }   // end MCMC
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   MCMC completed in", time_span.count(), "seconds");
    
    // create ret object
    Rcpp::List ret;
    ret.push_back(Rcpp::wrap( s0 ));
    ret.push_back(Rcpp::wrap( s ));
    ret.push_back(Rcpp::wrap( y ));
    ret.push_back(Rcpp::wrap( path_mat ));
    ret.push_back(Rcpp::wrap( friction ));
    ret.push_back(Rcpp::wrap( friction_best ));
    ret.push_back(Rcpp::wrap( SS_res_store ));
    ret.push_back(Rcpp::wrap( GLSpenalty ));
    
    
    Rcpp::StringVector ret_names;
    ret_names.push_back("s_initial");
    ret_names.push_back("s_final");
    ret_names.push_back("y");
    ret_names.push_back("path_mat");
    ret_names.push_back("friction");
    ret_names.push_back("friction_best");
    ret_names.push_back("SS_res");
    ret_names.push_back("GLSpenalty");
    
    ret.names() = ret_names;
    return ret;
}

//------------------------------------------------
// propose hexs by searching from focal hex with given depth
void prop_pointSearch(vector<int> &h_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth) {
    
    // initialise
    int nhex = neighbors.size();
    h_list.clear();
    
    // keep searching until at least one hex included
    while (h_list.size()==0) {
        
        // propose focal hex
        int h0 = sample2(1,nhex)-1;
        h_list.push_back(h0);
        
        // add to target list by sequentially searching neighbors
        for (int d=0; d<searchDepth; d++) {
            int hsize = h_list.size();
            for (int j=0; j<hsize; j++) {
                for (int i=0; i<int(neighbors[h_list[j]].size()); i++) {
                    int h = neighbors[h_list[j]][i];
                    if (hex_npath[h]>0) {
                        h_list.push_back(h);
                    }
                }
            }
        }
    }
    
    // remove duplicates
    sort(h_list.begin(), h_list.end());
    h_list.erase(unique(h_list.begin(), h_list.end()), h_list.end());
}

//------------------------------------------------
// propose hexs by searching in random path from focal hex
void prop_squiggle(vector<int> &h_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth) {
    
    // initialise
    int nhex = neighbors.size();
    h_list.clear();
    
    // keep searching until at least one hex included
    while (h_list.size()==0) {
        
        // propose starting hex
        int h = sample2(1,nhex)-1;
        //if (hex_npath[h]>0) {
            h_list.push_back(h);
        //}
        
        // add to target list by randomly searching neighbors
        for (int d=0; d<searchDepth; d++) {
            h = neighbors[h][sample2(1,neighbors[h].size())-1];
            //if (hex_npath[h]>0) {
                h_list.push_back(h);
            //}
        }
    }
    
    // remove duplicates
    sort(h_list.begin(), h_list.end());
    h_list.erase(unique(h_list.begin(), h_list.end()), h_list.end());
}

//------------------------------------------------
// find path
vector< vector<int> > findPath(vector<double> &distance, vector<int> &visited, vector<int> &traceback, vector<double> &friction, vector< vector<int> > &neighbors, int source, vector<int> &dest, vector<double> &s) {
    
    // reset vectors as needed
    fill(distance.begin(), distance.end(), OVERFLO);
    distance[source] = friction[source];
    fill(visited.begin(), visited.end(), 0);
    fill(traceback.begin(), traceback.end(), 0);
    
    // initialise objects for search
    int nhex = neighbors.size();    // total number of hexs
    int ndest = dest.size();    // number of destinations
    int nfound = 0;             // number of destinations found
    //int focus = source;         // which hex is the current focus
    //double focus_dist = 0;         // distance to current focus
    vector<int> tempList(1,source); // list to store nodes that have been explored but not yet confirmed
    
    // initialise return object (list of paths)
    vector< vector<int> > ret(ndest);
    
    // search until run out of hexs
    for (int rep=0; rep<nhex; rep++) {
        
        // find smallest distance in tempList
        double minDist = distance[tempList[0]];
        int minIndex = 0;
        for (int i=1; i<int(tempList.size()); i++) {
            if (distance[tempList[i]] < minDist) {
                minDist = distance[tempList[i]];
                minIndex = i;
            }
        }
        
        // move focal node to smallest node
        int focus = tempList[minIndex];
        double focus_dist = distance[focus];
        
        // if focus is one of target destinations
        for (int i=0; i<ndest; i++) {
            if (focus==dest[i]) {
                nfound ++;
                
                // add distance to s
                s[i] = minDist;
                
                // add entire path to ret
                int tracer = dest[i];
                ret[i].push_back(tracer);
                while (tracer != source) {
                    tracer = traceback[tracer];
                    ret[i].push_back(tracer);
                }
            }
            
            // break if no more destinations
            if (nfound==ndest) {
                break;
            }
        }
        
        // mark focus as visited
        visited[focus] = 2;
        tempList.erase(tempList.begin()+minIndex);
        
        // expand search to neighbors
        int nn = neighbors[focus].size();
        for (int i=0; i<nn; i++) {
            int nextHex = neighbors[focus][i];
            
            // if not finalised
            if (visited[nextHex] != 2) {
                double nextDist = focus_dist + friction[nextHex];
                
                // replace distance if this is an improvement
                //if (distance[nextHex]==OVERFLO || nextDist<distance[nextHex]) {
                if (nextDist<distance[nextHex]) {
                    distance[nextHex] = nextDist;
                    traceback[nextHex] = focus;
                }
                
                // mark as visited
                if (visited[nextHex]==0) {
                    visited[nextHex] = 1;
                    tempList.push_back(nextHex);
                }
            }
        }
        
    }   // end loop over reps
    
    return ret;
}

vector<int> findPath(vector<double> &distance, vector<int> &visited, vector<int> &traceback, vector<double> &friction, vector< vector<int> > &neighbors, int source, int dest, double &s) {
    vector<int> dest_vec(1,dest);
    vector<double> s_vec(1,s);
    vector< vector<int> > ret = findPath(distance, visited, traceback, friction, neighbors, source, dest_vec, s_vec);
    s = s_vec[0];
    return(ret[0]);
}


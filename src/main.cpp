
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

// TODO
void prop_squiggle2(vector<int> &h_list, vector<double> &delta_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth);

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
    
    // objects for storing results
    vector<double> SS_res_store(reps);
    vector<double> rSquared_store(reps);
    
    // begin MCMC
    double rSquared_old = 0;
    double SS_res_old = OVERFLO;
    vector<double> s_new = s;
    vector<int> h_list;
    vector<double> delta_list;
    double frictionDelta = 1;
    int nHits = 0;
    for (int rep=0; rep<reps; rep++) {
        //print(rep);
        
        // report progress
        if ( ((rep+1) % 1000)==0 ) {
            print("      iteration",rep+1);
        }
        
        // sharpen surface
        //if ( (nHits % 100)==0 && true) {
        if ( ((rep+1) % 1000)==0 && false) {
            print("   sharpening");
            nHits ++;
            
            // sharpen surface
            vector<double> friction2 = friction;
            for (int i=0; i<nhex; i++) {
                for (int j=0; j<int(hex_neighbors[i].size()); j++) {
                    friction2[i] += friction[hex_neighbors[i][j]];
                }
                friction2[i] /= double(hex_neighbors[i].size()+1);
            }
            friction = friction2;
            sort(friction2.begin(), friction2.end());
            double frictionCutoff = friction2[nhex*0.9];
            for (int i=0; i<nhex; i++) {
                if (friction[i] < frictionCutoff) {
                    friction[i] = 1;
                } else {
                    //friction[i] = frictionCutoff;
                }
            }
            
            // recalculate shortest paths
            for (int ind=0; ind<n2; ind++) {
                int i = index_i[ind];
                int j = index_j[ind];
                path_mat[ind] = findPath(distance, visited, traceback, friction, hex_neighbors, x[i], x[j], s[ind]);
            }
            
            // recalculate paths transecting each hex
            hex_path = vector< vector<int> >(nhex);
            hex_npath = vector<int>(nhex);
            for (int ind=0; ind<n2; ind++) {
                for (int k=0; k<int(path_mat[ind].size()); k++) {
                    int thisHex = path_mat[ind][k];
                    hex_path[thisHex].push_back(ind);
                    hex_npath[thisHex] ++;
                }
            }
            //break;
        }
        
        // propose by increasing friction around focal node
        //prop_pointSearch(h_list, hex_neighbors, hex_npath, 0);
        //prop_pointSearch(h_list, hex_neighbors, hex_npath, sample2(0,2));
        //prop_squiggle(h_list, hex_neighbors, hex_npath, sample2(0,20));
        //prop_squiggle(h_list, hex_neighbors, hex_npath, sample2(0,nhex/2));
        prop_squiggle2(h_list, delta_list, hex_neighbors, hex_npath, sample2(0,nhex/2));
        
        
        // change friction surface
        frictionDelta = 1 - 2*sample2(0,1);
        //frictionDelta = 1;
        for (int i=0; i<int(h_list.size()); i++) {
            int h = h_list[i];
            //friction[h] += delta_list[i];
            friction[h] += frictionDelta;
        }
        
        // get unique list of path indices that pass through target hexs
        vector<int> u;
        for (int i=0; i<int(h_list.size()); i++) {
            int h = h_list[i];
            if (hex_npath[h]>0) {
                push_back_multiple(u, hex_path[h]);
            }
        }
        sort(u.begin(), u.end());
        u.erase(unique(u.begin(), u.end()), u.end());
        
        // recalculate paths
        vector< vector<int> > store_path(u.size());
        for (int i=0; i<int(u.size()); i++) {
            int ind = u[i];
            store_path[i] = findPath(distance, visited, traceback, friction, hex_neighbors, x[index_i[ind]], x[index_j[ind]], s_new[ind]);
        }
        
        // calculate regression coefficients
        double Sx=0, Sy=0, Sxx=0, Syy=0, Sxy=0;
        for (int ind=0; ind<n2; ind++) {
            Sx += s_new[ind];
            Sy += y[ind];
            Sxx += s_new[ind]*s_new[ind];
            Syy += y[ind]*y[ind];
            Sxy += s_new[ind]*y[ind];
        }
        
        // calculate model fit
        double m = (0.5*n*(n-1)*Sxy - Sx*Sy)/(0.5*n*(n-1)*Sxx - Sx*Sx);
        double c = (Sy - m*Sx)/(0.5*n*(n-1));
        
        double SS_tot = Syy - Sy*Sy/(0.5*n*(n-1));
        double SS_res = Syy - 2*m*Sxy - 2*c*Sy + m*m*Sxx + 2*m*c*Sx + 0.5*n*(n-1)*c*c;
        
        double rSquared = 1 - SS_res/SS_tot;
        
        
        // if accept proposed change
        //if (rSquared > rSquared_old) {
        if (SS_res < SS_res_old) {
        //if (true) {
            
            // update fit
            nHits ++;
            rSquared_old = rSquared;
            SS_res_old = SS_res;
            s = s_new;
            
            // remove old paths from hexs
            for (int i=0; i<int(u.size()); i++) {
                for (int k=0; k<int(path_mat[u[i]].size()); k++) {
                    int h = path_mat[u[i]][k];
                    hex_path[h].erase(remove(hex_path[h].begin(), hex_path[h].end(), u[i]), hex_path[h].end());
                }
            }
            
            // add new paths to hexs
            for (int i=0; i<int(u.size()); i++) {
                for (int k=0; k<int(store_path[i].size()); k++) {
                    int h = store_path[i][k];
                    hex_path[h].push_back(u[i]);
                }
            }
            
            // replace path in path_mat
            for (int i=0; i<int(u.size()); i++) {
                path_mat[u[i]] = store_path[i];
            }
            
        }
        // reject proposed change
        else {
            
            // revert to previous
            s_new = s;
            
            // reverse changes to friction surface
            for (int i=0; i<int(h_list.size()); i++) {
                int h = h_list[i];
                //friction[h] -= delta_list[i];
                friction[h] -= frictionDelta;
            }
        }
        
        // store results of this iteration
        SS_res_store[rep] = SS_res_old;
        rSquared_store[rep] = rSquared_old;
        
    }   // end MCMC
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   MCMC completed in", time_span.count(), "seconds");
    
    // create ret object
    Rcpp::List ret;
    ret.push_back(Rcpp::wrap( s ));
    ret.push_back(Rcpp::wrap( path_mat ));
    ret.push_back(Rcpp::wrap( friction ));
    ret.push_back(Rcpp::wrap( SS_res_store ));
    ret.push_back(Rcpp::wrap( rSquared_store ));
    
    Rcpp::StringVector ret_names;
    ret_names.push_back("s");
    ret_names.push_back("path_mat");
    ret_names.push_back("friction");
    ret_names.push_back("SS_res");
    ret_names.push_back("rSquared");
    
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
// TODO
void prop_squiggle2(vector<int> &h_list, vector<double> &delta_list, vector< vector<int> > &neighbors, vector<int> &hex_npath, int searchDepth) {
    
    prop_squiggle(h_list, neighbors, hex_npath, searchDepth);
    vector<int> h_temp = h_list;
    delta_list = vector<double>(h_list.size(),1);
    
    prop_squiggle(h_list, neighbors, hex_npath, searchDepth);
    for (int i=0; i<int(h_list.size()); i++) {
        delta_list.push_back(-1);
    }
    push_back_multiple(h_temp, h_list);
    h_list = h_temp;
    
    // remove duplicates
    //sort(h_list.begin(), h_list.end());
    //h_list.erase(unique(h_list.begin(), h_list.end()), h_list.end());
}

//------------------------------------------------
// find path
vector< vector<int> > findPath(vector<double> &distance, vector<int> &visited, vector<int> &traceback, vector<double> &friction, vector< vector<int> > &neighbors, int source, vector<int> &dest, vector<double> &s) {
    
    // reset vectors as needed
    fill(distance.begin(), distance.end(), OVERFLO);
    distance[source] = 0;
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


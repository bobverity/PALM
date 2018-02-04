
#include <Rcpp.h>
#include <chrono>
#include "misc.h"
#include "probability.h"

using namespace std;

// find path
vector< vector<int> > findPath(vector<int> &distance, vector<int> &visited, vector<int> &traceback, vector<int> &friction, vector< vector<int> > &neighbors, int source, vector<int> &dest, vector<int> &s);

vector<int> findPath(vector<int> &distance, vector<int> &visited, vector<int> &traceback, vector<int> &friction, vector< vector<int> > &neighbors, int source, int dest, int &s);

//------------------------------------------------
// Run model
// [[Rcpp::export]]
Rcpp::List fitModel_cpp(Rcpp::List args_data, Rcpp::List args_model) {
    
    // extract data
    vector<int> x = Rcpp_to_vector_int(args_data["hex_data"]);
    vector< vector<int> > hex_neighbors = Rcpp_to_mat_int(args_data["hex_neighbors"]);
    vector< vector<double> > stat = Rcpp_to_mat_double(args_data["stat"]);
    int n = x.size();
    double n2 = 0.5*n*(n-1);
    int nhex = hex_neighbors.size();
    
    // extract model parameters
    int reps = Rcpp_to_int(args_model["reps"]);
    
    // start timer
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    
    // create objects for finding shortest path
    vector<int> friction(nhex,1);
    vector<int> distance(nhex);
    vector<int> visited(nhex);
    vector<int> traceback(nhex);
    
    // find initial shortest paths
    vector< vector< vector<int> > > path_mat(n, vector< vector<int> >(n));
    vector< vector<int> > s(n, vector<int>(n));
    for (int i=0; i<(n-1); i++) {
        path_mat[i] = findPath(distance, visited, traceback, friction, hex_neighbors, x[i], x, s[i]);
    }
    
    // store paths transecting each hex
    vector< vector< vector<int> > > hex_path(nhex);
    vector<int> hex_npath(nhex);
    for (int i=0; i<(n-1); i++) {
        for (int j=(i+1); j<n; j++) {
            for (int k=0; k<(s[i][j]+1); k++) {
                int thisHex = path_mat[i][j][k];
                hex_path[thisHex].push_back(vector<int>(2));
                hex_path[thisHex][hex_npath[thisHex]] = {i,j};
                hex_npath[thisHex] ++;
            }
        }
    }
    
    // objects for linear modelling
    double rSquared = 0;
    double logLike_old = -OVERFLO;
    double alpha = 0;
    double beta = 1;
    double sigma2 = 10;
    
    // objects for storing results
    vector<double> alpha_store(reps);
    vector<double> beta_store(reps);
    
    // begin MCMC
    vector< vector<int> > s_new = s;
    vector< vector< vector<int> > > store_path_list;
    for (int rep=0; rep<reps; rep++) {
        //print(rep);
        
        // propose hex to update
        int h0 = sample2(1,nhex)-1;
        while (hex_npath[h0]==0) {
            h0 = sample2(1,nhex)-1;
        }
        vector<int> h_list(1,h0);
        for (int i=0; i<int(hex_neighbors[h0].size()); i++) {
            //int h = hex_neighbors[h0][i];
            //if (hex_npath[h]>0) {
            //    h_list.push_back(h);
            //}
        }
        printVector(h_list);
        
        
        // change friction surface
        int friction_delta = 1;
        //int friction_delta = 1 - 2*sample2(0,1);
        //int friction_delta = sample2(1,10);
        //int friction_delta = sample2(-5,5);
        
        for (int hi=0; hi<int(h_list.size()); hi++) {
            int h = h_list[hi];
            friction[h] += friction_delta;
        }
        
        /*
        for (int i=0; i<int(hex_neighbors[h].size()); i++) {
            friction[hex_neighbors[h][i]] += friction_delta;
        }
        */
        
        // recalculate paths
        store_path_list.clear();
        for (int hi=0; hi<int(h_list.size()); hi++) {
            int h = h_list[hi];
            store_path_list.push_back(vector< vector<int> > (hex_npath[h]));
            for (int k=0; k<hex_npath[h]; k++) {
                int i = hex_path[h][k][0]; // source index
                int j = hex_path[h][k][1]; // destination index
                store_path_list[hi][k] = findPath(distance, visited, traceback, friction, hex_neighbors, x[i], x[j], s_new[i][j]);
            }
        }
        
        // calculate regression coefficients
        double Sx=0, Sy=0, Sxx=0, Syy=0, Sxy=0;
        for (int i=0; i<(n-1); i++) {
            for (int j=(i+1); j<n; j++) {
                Sx += s_new[i][j];
                Sy += stat[i][j];
                Sxx += s_new[i][j]*s_new[i][j];
                Syy += stat[i][j]*stat[i][j];
                Sxy += s_new[i][j]*stat[i][j];
            }
        }
        
        // calculate model fit
        double m = (0.5*n*(n-1)*Sxy - Sx*Sy)/(0.5*n*(n-1)*Sxx - Sx*Sx);
        double c = (Sy - m*Sx)/(0.5*n*(n-1));
        
        double SS_tot = Syy - Sy*Sy/(0.5*n*(n-1));
        double SS_res = Syy - 2*m*Sxy - 2*c*Sy + m*m*Sxx + 2*m*c*Sx + 0.5*n*(n-1)*c*c;
        
        double rSquared_new = 1 - SS_res/SS_tot;
        
        /*
        // calculate likelihood
        double logLike_new = -n2/2*log(sigma2);
        for (int i=0; i<(n-1); i++) {
            for (int j=(i+1); j<n; j++) {
                double xi = s_new[i][j];
                double yi = stat[i][j];
                logLike_new -= (yi-alpha-beta*xi)*(yi-alpha-beta*xi)/(2*sigma2);
            }
        }
        //logLike_new += -friction[h]*friction[h]/10 + (friction[h]-friction_delta)*(friction[h]-friction_delta)/10;
        */
        
        
        // if accept proposed change
        if (rSquared_new > rSquared) {
        //if ( log(runif_0_1()) < (logLike_new-logLike_old) ) {
        //if (logLike_new > logLike_old) {
            //foo();
            rSquared = rSquared_new;
            //logLike_old = logLike_new;
            s = s_new;
            
            for (int hi=0; hi<int(h_list.size()); hi++) {
                int h = h_list[hi];
                
                vector< vector<int> > hex_path_old = hex_path[h];
                for (int k=0; k<int(hex_path_old.size()); k++) {
                    int i = hex_path_old[k][0]; // source index
                    int j = hex_path_old[k][1]; // destination index
                    
                    // erase this path from all hexs
                    for (int k=0; k<int(path_mat[i][j].size()); k++) {  // (note re-use of k index)
                        int thisHex = path_mat[i][j][k];
                        for (int l=0; l<hex_npath[thisHex]; l++) {
                            if (hex_path[thisHex][l][0]==i && hex_path[thisHex][l][1]==j) {
                                hex_path[thisHex].erase(hex_path[thisHex].begin()+l);
                                break;
                            } else if (hex_path[thisHex][l][0]==j && hex_path[thisHex][l][1]==i) {
                                hex_path[thisHex].erase(hex_path[thisHex].begin()+l);
                                break;
                            }
                        }
                        hex_npath[thisHex] --;
                    }
                    
                    // add new path to all hexs
                    for (int l=0; l<int(store_path_list[hi][k].size()); l++) {
                        int thisHex = store_path_list[hi][k][l];
                        hex_path[thisHex].push_back(hex_path_old[k]);
                        hex_npath[thisHex] ++;
                    }
                    
                    // replace path in path_mat
                    path_mat[i][j] = store_path_list[hi][k];
                }
            }
            
        }
        // reject proposed change
        else {
            //bar();
            s_new = s;
            
            // reverse changes to friction surface
            for (int hi=0; hi<int(h_list.size()); hi++) {
                int h = h_list[hi];
                friction[h] -= friction_delta;
            }
            //for (int i=0; i<int(hex_neighbors[h].size()); i++) {
            //    friction[hex_neighbors[h][i]] -= friction_delta;
            //}
        }
        
        /*
        // Gibbs sample alpha
        double alpha_postMean = 0;
        double alpha_postVar = pow(sigma2/n2, 0.5);
        for (int i=0; i<(n-1); i++) {
            for (int j=(i+1); j<n; j++) {
                double xi = s[i][j];
                double yi = stat[i][j];
                alpha_postMean += yi - beta*xi;
            }
        }
        alpha_postMean /= n2;
        alpha = rnorm1(alpha_postMean, alpha_postVar);
        
        // Gibbs sample beta
        double beta_postMean = 0;
        double beta_postVar = 0;
        for (int i=0; i<(n-1); i++) {
            for (int j=(i+1); j<n; j++) {
                double xi = s[i][j];
                double yi = stat[i][j];
                beta_postMean += (yi-alpha)*xi;
                beta_postVar += xi*xi;
            }
        }
        beta_postMean /= beta_postVar;
        beta_postVar = pow(sigma2/beta_postVar, 0.2);
        beta = rnorm1(beta_postMean, beta_postVar);
        */
        // store results
        alpha_store[rep] = alpha;
        beta_store[rep] = beta;
        
    }   // end MCMC
    
    // end timer
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    chrono::duration<double> time_span = chrono::duration_cast< chrono::duration<double> >(t2-t1);
    print("   MCMC completed in", time_span.count(), "seconds");
    
    // create ret object
    Rcpp::List ret;
    ret.push_back(Rcpp::wrap( s ));
    ret.push_back(Rcpp::wrap( path_mat ));
    ret.push_back(Rcpp::wrap( hex_npath ));
    ret.push_back(Rcpp::wrap( friction ));
    ret.push_back(Rcpp::wrap( alpha_store ));
    ret.push_back(Rcpp::wrap( beta_store ));
    
    Rcpp::StringVector ret_names;
    ret_names.push_back("s");
    ret_names.push_back("path_mat");
    ret_names.push_back("hex_npath");
    ret_names.push_back("friction");
    ret_names.push_back("alpha");
    ret_names.push_back("beta");
    
    ret.names() = ret_names;
    return ret;
}

//------------------------------------------------
// find path
vector< vector<int> > findPath(vector<int> &distance, vector<int> &visited, vector<int> &traceback, vector<int> &friction, vector< vector<int> > &neighbors, int source, vector<int> &dest, vector<int> &s) {
    
    // reset vectors as needed
    fill(distance.begin(), distance.end(), 9999);
    fill(visited.begin(), visited.end(), 0);
    fill(traceback.begin(), traceback.end(), 0);
    
    // initialise objects for search
    int ndest = dest.size();
    int focus = source;
    int focus_dist = 0;
    vector<int> tempList(1,source);
    distance[source] = 0;
    int minIndex = 0;
    int minDist = 0;
    int nfound = 0;
    
    // initialise return object
    vector< vector<int> > ret(ndest);
    
    // search until run out of options
    for (int rep=0; rep<int(neighbors.size()); rep++) {
        
        // find smallest distance in tempList
        minDist = distance[tempList[0]];
        minIndex = 0;
        for (int i=1; i<int(tempList.size()); i++) {
            if (distance[tempList[i]] < minDist) {
                minDist = distance[tempList[i]];
                minIndex = i;
            }
        }
        
        // move focal node to smallest node
        focus = tempList[minIndex];
        focus_dist = distance[focus];
        
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
        }
        
        // break if no more destinations
        if (nfound==ndest) {
            break;
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
                int nextDist = focus_dist + friction[nextHex];
                
                // replace distance if this is an improvement
                if (distance[nextHex]==9999 || nextDist<distance[nextHex]) {
                    distance[nextHex] = nextDist;
                    traceback[nextHex] = focus; // store trace back
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

vector<int> findPath(vector<int> &distance, vector<int> &visited, vector<int> &traceback, vector<int> &friction, vector< vector<int> > &neighbors, int source, int dest, int &s) {
    vector<int> dest_vec(1,dest);
    vector<int> s_vec(1,s);
    vector< vector<int> > ret = findPath(distance, visited, traceback, friction, neighbors, source, dest_vec, s_vec);
    s = s_vec[0];
    return(ret[0]);
}

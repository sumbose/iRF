#ifndef RITmain
#define RITmain

#include <vector>
#include <algorithm>
#include <set>
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <math.h>
#include "RaggedArray.h"
#include "RITaux.h"

using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
set<vector<int>> RIT_basic(RaggedArray &x, NumericVector &weights, const int L,
        const double branch, const int depth, const int n_trees,
        const unsigned int min_inter_sz, const int n_cores, const int n) {

    // Set up parameters
    const int fl_branch = floor(branch);
    const int cl_branch = ceil(branch);
    const double branch_diff = branch - fl_branch;
    const int depthFinal = depth - 2;

    // Set up vector of seeds for RNG
    IntegerVector seeds(n_cores);
    seeds = INT_MAX * runif(n_cores);

    // Set up output objects
    // Union of candidate interactions for all trees
    set<vector<int>> total_candidate_interactions;

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

#pragma omp parallel
    {
        // Set up RNG for each thread
#ifdef _OPENMP
        //Use Mersenne Twister as RNG
        mt19937_64 mt(seeds[omp_get_thread_num()]);
#else
        //Use Mersenne Twister as RNG
        mt19937_64 mt(seeds[0]);
#endif

        discrete_distribution<int> r_obs(weights.begin(), weights.end());
        // Use for random number of branches
        uniform_real_distribution<> r_unif(0.0, 1.0);

#pragma omp for schedule(static) nowait
        for (int tree = 0; tree < n_trees; tree++) {

            vector<int> root;

            // First intersection computed by walking along arrays
            // as sets will be of similar size

            const int i1 = r_obs(mt);
            const int i2 = r_obs(mt);
            set_intersection(x.begin(i1), x.end(i1),
                    x.begin(i2), x.end(i2), back_inserter(root));

            if (root.size() < min_inter_sz) {
                // Interactions must have size at least min_inter_sz
                continue;
            }

            if ((root.size() <= min_inter_sz) || (depth <= 2)) {
                // depth >= 3
#pragma omp critical(update_total_candidate_interactions)
                {
                    total_candidate_interactions.insert(root);
                }
                continue;
            }

            // Only run this code when the initial intersection
            // produces an interaction of size greater than
            // min_inter_sz

            // Set of candidate interactions for each tree
            set<vector<int>> candidate_interactions;

            // Initialise parents
            vector<RaggedArray> parents(depthFinal);
            parents[0].push_back(root.begin(), root.end());
            for (int depth = 1; depth <= depthFinal; depth++) {
                for (int node = 0; node < parents[depth-1].nrow(); node++) {
                    int cur_branch;
                    /*
                     *  if(floor(branch) == branch) {
                     *      // If branch is an integer
                     *      cur_branch=branch;  
                     *  }
                     */
                    if (r_unif(mt) < branch_diff) {
                        // If a random number in (0, 1) is less than
                        // decimal part of branch
                        cur_branch = cl_branch;
                    } else {
                        // If a random number in (0, 1) is greater than
                        // decimal part of branch
                        cur_branch = fl_branch;
                    }

                    for (int k = 0; k < cur_branch; k++) {
                        const int i = r_obs(mt);
                        vector<int> temp_interaction = \
                            binary_intersect(x.begin(i), x.end(i),
                                    parents[depth-1].begin(node),
                                    parents[depth-1].end(node));

                        if (temp_interaction.size() < min_inter_sz) {
                            continue;
                        }

                        if ((temp_interaction.size() == min_inter_sz)
                                || (depth == depthFinal)) {
                            candidate_interactions.insert(temp_interaction);
                        } else {
                            parents[depth].push_back(temp_interaction.begin(),
                                    temp_interaction.end());
                        }
                    }
                }
            }

#pragma omp critical(update_total_candidate_interactions)
            {
                total_candidate_interactions.insert(
                        candidate_interactions.begin(),
                        candidate_interactions.end());
            }
        }
    }
    return total_candidate_interactions;
}


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
set<vector<int>> RIT_minhash(RaggedArray &x, const int L, const double branch,
        const int depth, const int n_trees, const double theta0,
        const double theta1, const unsigned int min_inter_sz,
        const int n_cores, const int n, int **H0t,
        const double n0_plus_1_over_n0, const double recip_n0_plus_1) {

    // Set up parameters
    const int fl_branch = floor(branch);
    const int cl_branch = ceil(branch);
    const double branch_diff = branch - fl_branch;
    const int depthFinal = depth - 2;

    // Set up vector of seeds for RNG
    IntegerVector seeds(n_cores);
    seeds = INT_MAX * runif(n_cores);

    // Set up output objects
    // Union of candidate interactions for all trees
    set<vector<int>> total_candidate_interactions;

#ifdef _OPENMP
    omp_set_num_threads(n_cores);
#endif

#pragma omp parallel
    {
        // Set up RNG for each thread
#ifdef _OPENMP
        // Use Mersenne Twister as RNG
        mt19937_64 mt(seeds[omp_get_thread_num()]);
#else
        // Use Mersenne Twister as RNG
        mt19937_64 mt(seeds[0]);
#endif

        uniform_int_distribution<int> r_obs(0, n-1);
        // Use for random number of branches
        uniform_real_distribution<> r_unif(0.0, 1.0);

#pragma omp for schedule(static) nowait
        for (int tree = 0; tree < n_trees; tree++) {

            vector<int> root;

            // First intersection computed by walking along arrays
            // as sets will be of similar size

            const int i1 = r_obs(mt);
            const int i2 = r_obs(mt);
            set_intersection(x.begin(i1), x.end(i1),
                    x.begin(i2), x.end(i2), back_inserter(root));

            if ((root.size() >= min_inter_sz)
                    && (PrevEst(root, H0t, L,
                            n0_plus_1_over_n0, recip_n0_plus_1) < theta0)) {
                // Class 0 prevalence must be low
                // interactions must have size at least min_inter_sz
                continue;
            }

            if ((root.size() <= min_inter_sz) || (depth <= 2)) {
                // depth >= 3
#pragma omp critical(update_total_candidate_interactions)
                {
                    total_candidate_interactions.insert(root);
                }
                continue;
            }

            // Only run this code when the initial intersection
            // produces an interaction of size greater than
            // min_inter_sz

            // Set of candidate interactions for each tree
            set<vector<int>> candidate_interactions;

            // Initialise parents
            vector<RaggedArray> parents(depthFinal);
            parents[0].push_back(root.begin(), root.end());
            for (int depth = 1; depth <= depthFinal; depth++) {
                for (int node = 0; node < parents[depth-1].nrow(); node++) {
                    int cur_branch;
                    if (r_unif(mt) < branch_diff) {
                        // If a random number in (0, 1) is less than
                        // decimal part of branch
                        cur_branch = cl_branch;
                    } else {
                        // If a random number in (0, 1) is greater than
                        // decimal part of branch
                        cur_branch = fl_branch;
                    }

                    for (int k = 0; k < cur_branch; k++) {
                        const int i = r_obs(mt);
                        vector<int> temp_interaction = \
                            binary_intersect(x.begin(i), x.end(i),
                                    parents[depth-1].begin(node),
                                    parents[depth-1].end(node));

                        if ((temp_interaction.size() < min_inter_sz)
                                || (PrevEst(temp_interaction, H0t, L,
                                        n0_plus_1_over_n0, recip_n0_plus_1)
                                    >= theta0)) {
                            continue;
                        }

                        if ((temp_interaction.size() == min_inter_sz)
                                || (depth == depthFinal)) {
                            candidate_interactions.insert(temp_interaction);
                        } else {
                            parents[depth].push_back(temp_interaction.begin(),
                                    temp_interaction.end());
                        }
                    }
                }
            }

#pragma omp critical(update_total_candidate_interactions)
            {
                total_candidate_interactions.insert(
                        candidate_interactions.begin(),
                        candidate_interactions.end());
            }
        }
    }
    return total_candidate_interactions;
}

#endif

#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat Q(arma::rowvec x_ik, arma::vec beta, double time) {

    arma::mat betaMat = arma::reshape(beta, 4, 2);
    betaMat = betaMat.t();
    
    arma::rowvec x_ik_full = {1, x_ik(0), x_ik(1), time};
    double q_x = arma::dot(x_ik_full, betaMat.row(0));
    double q_t = arma::dot(x_ik_full, betaMat.row(1));
    
    arma::mat qmat = {{0, exp(q_x),        0},
                      {0,        0, exp(q_t)},
                      {0,        0,        0}};
    arma::vec q_row_sums = arma::sum(qmat, 1);
    qmat.diag() = -q_row_sums;
    
    return qmat;
}

// [[Rcpp::export]]
double fn_log_post(arma::vec pars, arma::field<arma::vec> prior_par, 
                   arma::field<arma::uvec> par_index, arma::mat x, 
                   arma::vec y, arma::vec t, arma::vec id, arma::vec eids,
                   arma::vec disc_t) {
    
    arma::rowvec init = {1, 0, 0};
    arma::mat resp_fnc = arma::eye(3,3);
    
    arma::vec beta = pars.elem(par_index(0) - 1); 
    
    arma::vec log_total_val_vec(eids.n_elem, arma::fill::zeros);
    
    omp_set_num_threads(10);
    # pragma omp parallel for
    for(int i = 0; i < eids.n_elem; i++) {

        arma::rowvec val;
        double log_norm = 0;

        arma::uvec sub_ind = arma::find(id == eids(i));
        arma::vec y_i = y.elem(sub_ind);
        arma::mat x_i = x.rows(sub_ind);
        arma::vec t_i = t.elem(sub_ind);
        arma::vec disc_t_i = disc_t.elem(sub_ind);
        
        arma::vec diag_part = resp_fnc.col(y_i(0) - 1);
        arma::mat D = arma::diagmat(diag_part);
        arma::mat f_i = init * D;

        for(int k = 1; k < t_i.n_elem; k++) { 
            
            // Time Homogeneous Transition Matrix -----------------------
            double t_diff = t_i(k) - t_i(k-1);
            arma::mat Q_mat = Q(x_i.row(k-1), beta, disc_t_i(k-1));
            Q_mat = t_diff * Q_mat;
            arma::mat P = arma::expmat(Q_mat);
            
            // Checking if the row is observed or a censored row
            if(y_i(k) != 99) {
                arma::vec diag_part = resp_fnc.col(y_i(k) - 1);
                arma::mat D = arma::diagmat(diag_part);
                val = f_i * P * D;

                arma::mat diag_val = arma::diagmat(val);
                arma::rowvec val_2 = val * diag_val;
                double norm_val = sqrt(arma::accu(val_2));

                f_i = val / norm_val;
                log_norm = log_norm + log(norm_val);
            } else { // censor row (only happens when disc == TRUE)
                val = f_i * P;

                arma::mat diag_val = arma::diagmat(val);
                arma::rowvec val_2 = val * diag_val;
                double norm_val = sqrt(arma::accu(val_2));

                f_i = val / norm_val;
                log_norm = log_norm + log(norm_val);
            }
            
        }
        
        log_total_val_vec(i) = log(arma::accu(f_i)) + log_norm;
    }
    
    double log_total_val = arma::accu(log_total_val_vec);
    
    arma::vec prior_mean = prior_par(0);
    arma::vec prior_sd_diag = prior_par(1);
    arma::mat prior_sd = arma::diagmat(prior_sd_diag);

    arma::vec prior = dmvnorm(pars.t(), prior_mean, prior_sd, true);
    double log_prior_dens = arma::as_scalar(prior);

    return log_total_val + log_prior_dens;
}
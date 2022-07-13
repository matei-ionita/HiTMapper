#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List som(arma::mat& data, int nx, int ny, int C=10, double lr=0.2) {
	int n = data.n_rows, d = data.n_cols, m = nx*ny, i, j, k, bmu, in;
	double dist, dist_curr, dif;
	arma::vec mapping = arma::vec(n);
	arma::rowvec diff = arma::rowvec(), update = arma::rowvec();
	arma::mat node_pos = arma::mat(m,d);

	// random initialization
	IntegerVector ind = Rcpp::sample(n, m, false);
	for (j=0; j<m; j++) {
		// Rcout << j << " " << ind(j) << "\n";
		node_pos.row(j) = data.row(ind(j) - 1);
	}

	double* p_data = data.memptr();
	double* p_node = node_pos.memptr();

	// som learning
	ind = Rcpp::sample(n, C*n, true)-1;

	for (i=0; i<C*n; i++) {
		in = ind(i);
		// find bmu (best matching unit)
		dist = arma::datum::inf;

		for (j=0; j<m; j++) {

			dist_curr=0;
			for (k=0; k<d; k++) {
				dif = p_data[in+k*n] - p_node[j+k*m];
				dist_curr += dif*dif;
			}

			if (dist_curr < dist) {
				dist = dist_curr;
				bmu = j;
			}
		}

		// update bmu position
		for (k=0; k<d; k++) {
			dif = p_data[in+k*n] - p_node[bmu+k*m];
			p_node[bmu+k*m] += lr * dif;
		}
	}

	// map data to nodes
	for (i=0; i<n; i++) {
		// find bmu (best matching unit)
		dist = arma::datum::inf;

		for (j=0; j<m; j++) {
			dist_curr=0;
			for (k=0; k<d; k++) {
				dif = p_data[i+k*n] - p_node[j+k*m];
				dist_curr += dif*dif;
			}

			if (dist_curr < dist) {
				dist = dist_curr;
				bmu = j;
			}
		}

		// assign data point to bmu
		mapping(i) = bmu+1;
	}

	return Rcpp::List::create(Rcpp::Named("node_pos")=node_pos,
							  Rcpp::Named("mapping")=mapping);
}

#include "RcppArmadillo.h"
using namespace Rcpp;


// [[Rcpp::export]]
List assign_datapoints(arma::mat& data, arma::mat& centroids) {
	int n = data.n_rows, d = data.n_cols, m = centroids.n_rows;
	double dist_curr, dif;
	IntegerVector mapping = IntegerVector(n);
	arma::mat sim = arma::mat(m,m, arma::fill::zeros);
	double* p_data = data.memptr();
	double* p_node = centroids.memptr();

	NumericVector dist = NumericVector(m);

	for (int i=0; i<n; i++) {
		// find bmu (best matching unit)

		for (int j=0; j<m; j++) {
			dist_curr=0;
			for (int k=0; k<d; k++) {
				dif = p_data[i+k*n] - p_node[j+k*m];
				dist_curr += dif*dif;
			}
			dist(j) = dist_curr;
		}

		int bmu = which_min(dist);
		dist(bmu) = arma::datum::inf;
		int bmu2 = which_min(dist);
		sim(bmu, bmu2) += 1;
		sim(bmu2, bmu) += 1;

		// assign data point to bmu
		mapping(i) = bmu+1;
	}

	return List::create(Named("mapping") = mapping, _["sim"] = sim);
}


arma::mat som(arma::mat& data, arma::uvec& bin, int m, 
	int n_passes=10, double lr=0.2) {
	int n = bin.n_elem, n0 = data.n_rows, d = data.n_cols, i, j, k, bmu, in;
	double dist, dist_curr, dif;
	arma::rowvec diff = arma::rowvec(), update = arma::rowvec();
	arma::mat node_pos = arma::mat(m,d);

	// random initialization
	IntegerVector ind = Rcpp::sample(n, m, false);
	for (j=0; j<m; j++)
		node_pos.row(j) = data.row(bin(ind(j)-1));

	double* p_data = data.memptr();
	double* p_node = node_pos.memptr();

	// som learning

	for (i=0; i<n_passes*n; i++) {
		IntegerVector s = Rcpp::sample(n,1)-1;
		in = bin(s(0));

		// find bmu (best matching unit)
		dist = arma::datum::inf;

		for (j=0; j<m; j++) {

			dist_curr=0;
			for (k=0; k<d; k++) {
				dif = p_data[in+k*n0] - p_node[j+k*m];
				dist_curr += dif*dif;
			}

			if (dist_curr < dist) {
				dist = dist_curr;
				bmu = j;
			}
		}

		// update bmu position
		for (k=0; k<d; k++) {
			dif = p_data[in+k*n0] - p_node[bmu+k*m];
			p_node[bmu+k*m] += lr * dif;
		}
	}

	return(node_pos);
}


arma::rowvec get_mean_manual(arma::mat& data, arma::uvec& bin) {
	int n = bin.n_elem, d = data.n_cols;
	arma::rowvec mean_vec(d, arma::fill::zeros);

	for (int i=0; i<n; i++) {
		int ind = bin(i);
		for (int j=0; j<d; j++)
			mean_vec(j) += data(ind,j);
	}
	mean_vec /= n;
	return mean_vec;
}


arma::mat get_centroids(arma::mat& data, arma::uvec& mapping, 
	arma::vec& all_k, int n_passes) {

	int n_bins = all_k.n_elem, i, k, d=data.n_cols, m=sum(all_k), start=0;
	arma::uvec un = unique(mapping);
	arma::mat centroids = arma::mat(m, d);

	for (i=0; i<n_bins; i++) {
		k = all_k(i);
		arma::uvec bin = find(mapping==un(i));

		if (k==0)
			continue;

		if (k==1)
			centroids.row(start) = get_mean_manual(data, bin);
		else
			centroids.rows(start,start+k-1) = som(data, bin, k, n_passes);

		start = start + k;
	}

	return centroids;
}


arma::uvec binning(arma::mat& filter, arma::vec& b1, arma::vec& b2) {
	int n = filter.n_rows, m1 = b1.n_elem, m2 = b2.n_elem, i, j, k;
	double val;
	arma::uvec bin_mapping = arma::uvec(n);

	for (i=0; i<n; i++) {
		val = filter(i,0);
		j=0;
		while (val > b1(j) && j < m1-1)
			j++;

		val = filter(i,1);
		k=0;
		while (val > b2(k) && k < m2-1) {
			k++;
		}

		bin_mapping(i) = j + k*m1;
	}

	return bin_mapping;
}


// [[Rcpp::export]]
arma::mat compute_centroids(arma::mat& data, arma::ivec& mapping, int k) {
	int d = data.n_cols, i;
	arma::mat centroids = arma::mat(k,d);

	for (i=0; i<k; i++) {
		arma::uvec idx = find(mapping == i+1);
		centroids.row(i) = get_mean_manual(data, idx);
	}

	return centroids;
}


arma::mat get_PCA(arma::mat& data, arma::mat& cova) {
	int d = data.n_cols;
	arma::mat means = mean(data), eigvec;

	arma::vec eigval;
	eig_sym(eigval, eigvec, cova);

	arma::mat scores = data * eigvec.cols(d-2,d-1);
	arma::mat diff = means * eigvec.cols(d-2,d-1);
	scores.each_row() -= diff;

	return scores;
}


arma::vec get_boundaries(arma::mat& filter, arma::uvec& grid_size, int ind) {
	int d = grid_size(ind), i;
	double m=min(filter.col(ind)), M=max(filter.col(ind));
	double width=(M-m)/d;
	arma::vec b(d);

	for (i=0; i<d; i++) {
		b(i) = m + (i+1)*width;
	}
	return b;
}


arma::vec allocate(int total_nodes, arma::vec max_bin, 
	arma::vec weight, arma::vec alloc, int depth) {

	if (sum(alloc)>= total_nodes || depth > 10)
		return(round(alloc));

	double to_alloc = total_nodes - sum(alloc);
	arma::uvec not_at_max = find(alloc != max_bin);
	double renorm = sum(weight(not_at_max));
	arma::vec alloc_new = min(weight * to_alloc/renorm, max_bin-alloc);

	return allocate(total_nodes, max_bin, weight, alloc+alloc_new, depth+1);
}

arma::vec allocate_nodes(arma::mat&data, arma::uvec& mapping, 
	int min_node_size, int total_nodes, double frac_sum_sq) {

	arma::uvec un = unique(mapping);
	arma::uvec bin;
	int n_bins = un.n_elem, i, l, j, k, d=data.n_cols;
	double v, s, sq, val;

	arma::vec bin_length = arma::vec(n_bins);
	arma::vec bin_var = arma::vec(n_bins);
	arma::vec log_size = arma::vec(n_bins);

	for (i=0; i<n_bins; i++) {
		bin = find(mapping==un(i));
		l = bin.n_elem;
		bin_length(i) = l;

		if (l < 2) {
			bin_var(i) = 0;
			log_size(i)=0;
		}
		else {
			v=0;
			for (k=0;k<d;k++) {
				sq=0;
				s=0;
				for (j=0;j<l;j++) {
					val=data(bin(j),k);
					s+=val;
					sq+=val*val;
				}
				s/=l;
				sq/=l;
				v+= (sq-s*s) * l/(l-1);
			}
			bin_var(i) = sqrt(v);
			log_size(i) = log2(l);
		}
	}

	arma::vec max_bin = ceil(bin_length / min_node_size);
	arma::vec weight = frac_sum_sq * (bin_var / sum(bin_var)) +
						(1-frac_sum_sq) * (log_size / sum(log_size));
	arma::vec alloc(n_bins, arma::fill::zeros);
	return allocate(total_nodes, max_bin, weight, alloc, 1);
}


arma::uvec get_bin_mapping(arma::mat& data, arma::mat& cova, arma::uvec& grid_size) {
	arma::mat filter = get_PCA(data, cova);
	arma::vec b1 = get_boundaries(filter, grid_size, 1);
	arma::vec b2 = get_boundaries(filter, grid_size, 0);
	arma::uvec bin_mapping = binning(filter, b1, b2);
	return(bin_mapping);
}


// [[Rcpp::export]]
arma::mat clustering_main(arma::mat& data, arma::mat& cova, arma::uvec& grid_size, 
	int total_nodes, int min_node_size, int n_passes) {

	arma::uvec bin_mapping = get_bin_mapping(data, cova, grid_size);
	arma::vec all_k = allocate_nodes(data, bin_mapping, min_node_size, total_nodes, 1);
	arma::mat centroids = get_centroids(data, bin_mapping, all_k, n_passes);
	return centroids;
}





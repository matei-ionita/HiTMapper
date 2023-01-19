#include "RcppArmadillo.h"
using namespace Rcpp;


// [[Rcpp::export]]
IntegerVector predict_datapoints(arma::mat& data, arma::mat& centroids) {
	int n = data.n_rows, d = data.n_cols, m = centroids.n_rows, i, j, k, bmu;
	double dist, dist_curr, dif;
	IntegerVector mapping = IntegerVector(n);
	double* p_data = data.memptr();
	double* p_node = centroids.memptr();

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

	return mapping;
}



arma::mat get_weights(arma::mat& sim, IntegerVector mapping) {
	int n = mapping.length(), m = sim.n_rows;
	arma::uvec tab = arma::uvec(m, arma::fill::zeros);
	arma::vec norms = arma::vec(m), col = arma::vec(m);
	arma::mat weights = arma::mat(m,m);

	for (int i=0; i<n; i++)
		tab(mapping(i)-1)++;

	for (int i=0; i<m; i++)
		for (int j=0; j<m; j++) {
			sim(i,j) /= tab(i)+tab(j);
			if (i==j)
				sim(i,j)=1;
		}

	for (int i=0; i<m; i++) {
		col = sim.col(i);
		norms(i) = sqrt(sum(col % col));
	}

	for (int i=0; i<m; i++)
		for (int j=i; j<m; j++) {
			double weight = sum(sim.col(i) % sim.col(j));
			weight /= norms(i) * norms(j);

			if (weight < 0.05 && sim(i,j)==0)
				weight = 0;

			weights(i,j) = weight;
			weights(j,i) = weight;
		}

	return weights;
}


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
		double d1 = dist(bmu);
		dist(bmu) = arma::datum::inf;
		int bmu2 = which_min(dist);
		double d2 = dist(bmu2);
		sim(bmu, bmu2) += exp(-(d2-d1)*(d2-d1));
		sim(bmu2, bmu) += exp(-(d2-d1)*(d2-d1));

		// assign data point to bmu
		mapping(i) = bmu+1;
	}

	arma::mat weights = get_weights(sim, mapping);

	return List::create(Named("mapping") = mapping, _["weights"] = weights);
}




arma::mat km(arma::mat& data, arma::uvec& bin, int m, int n_passes, double tol=0.01) {
	int n = bin.n_elem, n0 = data.n_rows, d = data.n_cols, i, j, k, in, curr, bmd;
	double dist, dist_curr, dif;
	arma::mat node_pos = arma::mat(m,d), node_pos_old = arma::mat(m,d);
	arma::uvec bmu = arma::uvec(n), occu = arma::uvec(n), outliers = arma::uvec();
	arma::vec bd = arma::vec(n);
	arma::uvec count = arma::uvec(m, arma::fill::zeros);
	arma::vec quant = {0.25,0.75};

	// random initialization
	IntegerVector ind = Rcpp::sample(n, m, false);
	for (j=0; j<m; j++)
		node_pos.row(j) = data.row(bin(ind(j)-1));

	double* p_data = data.memptr();
	double* p_node = node_pos.memptr();

	for (int iter=0; iter<n_passes; iter++) {

		// assign each datapoint to best matching centroid
		for (i=0; i<n; i++) {
			in = bin(i);
			dist = arma::datum::inf;

			for (j=0; j<m; j++) {

				dist_curr=0;
				for (k=0; k<d; k++) {
					dif = p_data[in+k*n0] - p_node[j+k*m];
					dist_curr += dif*dif;
				}

				if (dist_curr < dist) {
					dist = dist_curr;
					bmu(i) = j;
					bd(i) = dist;
				}
			}
		}

		arma::vec q = quantile(bd, quant);
		double upper = q(1) + 1.5*(q(1)-q(0));

		count *= 0;
		for (i=0; i<n; i++) {
			if (bd(i) <= upper)
				count(bmu(i))++;
		}

		// update centroid position
		node_pos *= 0;
		for (i=0; i<n; i++) {
			curr = bmu(i);
			in = bin(i);

			if (bd(i) > upper)
				continue;

			for (k=0; k<d; k++) {
				p_node[curr+k*m] += p_data[in+k*n0] / count(curr);
			}
		}

		occu = occu * 0;

		for (j=0; j<m; j++) {
			dist = arma::datum::inf;
			for (i=0; i<n; i++) {
				in = bin(i);
				dist_curr=0;

				for (k=0; k<d; k++) {
					dif = p_data[in+k*n0] - p_node[j+k*m];
					dist_curr += dif*dif;
				}

				if (dist_curr < dist) {
					dist = dist_curr;
					bmd=i;
				}
			}

			while(occu(bmd)==1) {
				ind = Rcpp::sample(n, 1, false);
				bmd=ind(0)-1;
			}

			node_pos.row(j) = data.row(bin(bmd));
			occu(bmd)=1;
		}

		if (iter > 0 && accu(abs(node_pos - node_pos_old)) < tol)
			break;

		node_pos_old = node_pos;

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
	arma::uvec un = unique(mapping), bin;
	arma::mat centroids = arma::mat(m, d);

	for (i=0; i<n_bins; i++) {
		k = all_k(i);
		bin = find(mapping==un(i));

		if (k==0)
			continue;

		centroids.rows(start, start+k-1) = km(data, bin, k, n_passes);

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
	arma::uvec idx;

	for (i=0; i<k; i++) {
		idx = find(mapping == i+1);
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

arma::vec allocate_nodes(arma::mat& data, arma::uvec& mapping, 
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





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
		// find bmu (best matching unit
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


// [[Rcpp::export]]
arma::mat get_edgelist(arma::mat& idx) {
	int n = idx.n_rows, m = idx.n_cols, cnt=0;
	arma::mat edgelist = arma::mat(n*m,3);
	arma::rowvec node_nbrs, nbr_nbrs, inter;

	for (int node=0; node<n; node++) {
		node_nbrs = idx.row(node);
		for (int i=0; i<m; i++) {
			int nbr = idx(node,i)-1;
			nbr_nbrs = idx.row(nbr);
			inter = intersect(node_nbrs, nbr_nbrs);

			if (node <= nbr || !any(nbr_nbrs==(node+1))) {
				double n_inter = inter.n_elem;
				edgelist(cnt,0) = node+1;
				edgelist(cnt,1) = nbr+1;
				edgelist(cnt,2) = n_inter / (2*m-n_inter);
				cnt++;
			}
		}
	}

	edgelist = edgelist.rows(0,cnt-1);

	return edgelist;
}


// [[Rcpp::export]]
arma::ivec sample_cells(arma::ivec& mapping, arma::ivec& uniq, int m) {
 
    int start = 0, k = uniq.n_elem, i, j;
    arma::ivec chosen = arma::ivec(k*m), tmp;
    IntegerVector ind;
    arma::uvec idx;
  
    for (j=0; j<k; j++) {
      	idx = find(mapping==uniq(j));
      	int n_idx = idx.n_elem;
      	int n_sam = std::min(n_idx,m);
      	ind = Rcpp::sample(n_idx, n_sam, false)-1;
    
      	for (i=0; i<n_sam; i++)
        	chosen(i+start) = idx(ind(i))+1;
        	start += n_sam;
        }
  
    chosen = chosen.subvec(0, start-1);
    return(chosen);
}


arma::mat som(arma::mat& data, arma::uvec& bin, int m, int n_passes=10, double lr=0.05) {
  int n = bin.n_elem, n0 = data.n_rows, d = data.n_cols, i, j, k, bmu, in, s=100;
  double dist, dist_curr, dif;
  arma::rowvec diff = arma::rowvec(), update = arma::rowvec();
  arma::mat node_pos = arma::mat(m,d);
  
  // random initialization
  IntegerVector ind = Rcpp::sample(n, m, false);
  for (j=0; j<m; j++) {
    node_pos.row(j) = data.row(bin(ind(j) - 1));
  }
  
  double* p_data = data.memptr();
  double* p_node = node_pos.memptr();
  
  // som learning
  ind = IntegerVector(s);
  
  for (i=0; i<n_passes*n; i++) {
  	if (i % s == 0)
  		ind = Rcpp::sample(n, s, true)-1;

  	lr = 0.05 - 0.04 * (i+0.0001)/(n_passes*n);

    in = bin(ind(i % s));

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


arma::rowvec get_mean_inplace(arma::mat& data, arma::uvec& bin) {
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


// [[Rcpp::export]]
arma::mat compute_centroids(arma::mat& data, arma::ivec& mapping, int k) {
	int d = data.n_cols, i;
	arma::mat centroids = arma::mat(k,d);
	arma::uvec idx;

	for (i=0; i<k; i++) {
		idx = find(mapping == i+1);
		centroids.row(i) = get_mean_inplace(data, idx);
	}

	return centroids;
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

		centroids.rows(start, start+k-1) = som(data, bin, k, n_passes);

		start = start + k;
	}

	return centroids;
}


arma::uvec binning(arma::vec& filter1, arma::vec& filter2, arma::vec& b1, arma::vec& b2) {
	int n = filter1.n_elem, m1 = b1.n_elem, m2 = b2.n_elem, i, j, k;
	double val;
	arma::uvec bin_mapping = arma::uvec(n);

	for (i=0; i<n; i++) {
		val = filter1(i);
		j=0;
		while (val > b1(j) && j < m1-1)
			j++;

		val = filter2(i);
		k=0;
		while (val > b2(k) && k < m2-1) {
			k++;
		}

		bin_mapping(i) = j + k*m1;
	}

	return bin_mapping;
}


arma::vec get_PCA(arma::mat& data, arma::mat& means, arma::mat& cova, int ind) {
	arma::mat eigvec;
	arma::vec eigval;
	eig_sym(eigval, eigvec, cova);

	arma::vec filter = data * eigvec.col(ind);
	arma::vec diff = means * eigvec.col(ind);
	filter -= diff(0);
	return filter;
}


arma::vec get_boundaries(arma::vec& filter, int siz) {
	arma::vec tmp = {0.02, 0.98};
	arma::vec q = quantile(filter,tmp);
	arma::vec b = arma::linspace(q(0), q(1), siz+1);
	b = b(arma::span(1,siz));
	return b;
}


arma::uvec get_bin_mapping(arma::mat& data, arma::mat& cova, arma::uvec& grid_size) {
	int d = cova.n_cols;
	arma::mat means = mean(data);
	arma::vec filter1 = get_PCA(data, means, cova, d-1);
	arma::vec filter2 = get_PCA(data, means, cova, d-2);
	arma::vec b1 = get_boundaries(filter1, grid_size(0));
	arma::vec b2 = get_boundaries(filter2, grid_size(1));
	arma::uvec bin_mapping = binning(filter1, filter2, b1, b2);
	return(bin_mapping);
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
	int min_node_size, int total_nodes) {

	arma::uvec un = unique(mapping), bin;
	int n_bins = un.n_elem, i, l;

	arma::vec bin_length = arma::vec(n_bins);

	for (i=0; i<n_bins; i++) {
		bin = find(mapping==un(i));
		l = bin.n_elem;
		bin_length(i) = l;
	}

	arma::vec max_bin = ceil(bin_length / min_node_size);
	arma::vec weight = arma::vec(n_bins, arma::fill::ones) / n_bins; 
	arma::vec alloc(n_bins, arma::fill::zeros);
	return allocate(total_nodes, max_bin, weight, alloc, 1);
}


// [[Rcpp::export]]
arma::mat clustering_main(arma::mat& data, arma::mat& cova, arma::uvec& grid_size, 
	int total_nodes, int min_node_size, int n_passes) {

	arma::uvec bin_mapping = get_bin_mapping(data, cova, grid_size);
	arma::vec all_k = allocate_nodes(data, bin_mapping, min_node_size, total_nodes);
	arma::mat centroids = get_centroids(data, bin_mapping, all_k, n_passes);
	return centroids;
}




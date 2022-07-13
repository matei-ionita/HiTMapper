#include "RcppArmadillo.h"
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List intersections(List nodes) {
	int n = nodes.length();
	double l1, l2, l12;
	IntegerVector node1, node2;

	std::vector<int> source, target;
	std::vector<double> iou;

	for (int i=0; i<n-1; i++) {
		node1=nodes[i];
		l1=node1.length();

		for (int j=i+1; j<n; j++) {
			node2=nodes[j];
			l2=node2.length();

			std::vector<int> intersect;
			std::set_intersection(node1.begin(), node1.end(), 
				node2.begin(), node2.end(),
				std::back_inserter(intersect));
			l12=intersect.size();

			if (l12>0) {
				source.push_back(i+1);
				target.push_back(j+1);
				iou.push_back(l12/(l1+l2-l12));
			}
		}
	}

	return Rcpp::List::create(Rcpp::Named("source")=source,
							  Rcpp::Named("target")=target,
							  Rcpp::Named("iou")=iou);
}
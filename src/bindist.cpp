#include <Rcpp.h>
using namespace Rcpp;

//' @export

// [[Rcpp::export]]
Rcpp::List bindist(NumericMatrix binmat, NumericMatrix data) {

  Rcpp::List nodes(binmat.nrow()); // create list object to receive the binned datapoints
  Rcpp::List nodedist(binmat.nrow());

  for(int i = 0; i < data.nrow(); i++)
    {
      NumericVector dists(binmat.nrow());
      std::vector<int> winners;
      std::vector<double> windist;
      for(int n = 0; n < binmat.nrow(); n++)
        {
          dists[n] = sqrt(sum(pow(binmat.row(n)-data.row(i), 2.0))); // get dist of datapoint i to each node n in binmat (why use euclidean distance here???)
        }
      // put the data index in the correct node bin
      if(nodes[which_min(dists)] == R_NilValue)
        {
          int winners = i + 1;
          nodes[which_min(dists)] = winners;
        }
      else
        {
          winners = as<std::vector<int> >(nodes[which_min(dists)]); // get current inhabitants of winning node
          winners.push_back(i + 1); // append winner index to end of vector. i + 1 to correct for 0-indexing in cpp
          nodes[which_min(dists)] = winners; // replace winning node vector with updated vector
        }
      if(nodedist[which_min(dists)] == R_NilValue)
        {
          double windist = min(dists);
          nodedist[which_min(dists)] = windist;
        }
      else
        {
          windist = as<std::vector<double> >(nodedist[which_min(dists)]); // repeat process to store the calculated distance
          windist.push_back(min(dists));
          nodedist[which_min(dists)] = windist;
        }
    }
  Rcpp::List out = Rcpp::List::create(Rcpp::Named("node.pops") = nodes,
                                      Rcpp::Named("node.dists") = nodedist);
  return out;
}


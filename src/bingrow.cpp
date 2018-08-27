#include <Rcpp.h>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/count.hpp>
#include <boost/accumulators/statistics/sum.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/skewness.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <random>
#include <numeric>

using namespace Rcpp;
using namespace boost::accumulators;
using namespace boost::math;

//
// dist.thresh <- 4
// markers <- markers.norm
//  growfact <- 2
// minpts <- 5
// data <- expdata.norm.sub

// [[Rcpp::plugins(cpp17)]]
std::vector<double> markdist(std::vector<double> markdat, int growfact)
{

  // create boost accumulator to calc variables used to model the distribution of marker j expression in node.
  accumulator_set< double, features<tag::mean, tag::count, tag::sum, tag::variance, tag::immediate_mean,
                                    tag::skewness, tag::moment<2>, tag::moment<3> > > nacc;
  nacc = std::for_each(markdat.begin(), markdat.end(), nacc); // add the data to the accumulator
  // create convenient extractors of necessary data.
  boost::accumulators::extractor<boost::accumulators::tag::mean> mean_;
  boost::accumulators::extractor<boost::accumulators::tag::variance> stdev_;
  boost::accumulators::extractor<boost::accumulators::tag::skewness> skew_;

  // Setup generators
  std::random_device rand;
  std::default_random_engine noise;

  // Sample from a uniform distribution i.e. [0,1)
  std::uniform_real_distribution<double> uniform_dist(0,1.0);

  auto markdist = boost::math::skew_normal_distribution<> (mean_(nacc), std::sqrt(stdev_(nacc)), skew_(nacc));

  std::vector<double> nmarks(growfact);
  for(int x = 0; x < growfact; x++)
  {
    noise.seed(rand());
    auto prob = uniform_dist(noise);
    nmarks[x] = boost::math::quantile(markdist, prob);
  }

  return nmarks;
}

std::vector<double> resamp(std::vector<double> markdat, int growfact)
{
  std::default_random_engine generator;
  std::uniform_int_distribution<int> distribution(0,markdat.size());

  std::vector<double> nmarks(growfact);
  for(int x = 0; x < growfact; x++)
  {
    nmarks[x] = markdat[distribution(generator)];
  }
  return(nmarks);
}

// [[Rcpp::export]]
NumericMatrix bingrow(double dthresh, int minpts, NumericMatrix data, int lmarkers, NumericMatrix oldbins, Rcpp::List ndist, Rcpp::List nlist, int growfact) {

  NumericVector binmat;

  for(int i = 0; i < ndist.length(); i++)
  {
    std::vector<double> node;
    node = ndist[i];

    std::vector<int> members;
    members = nlist[i];

    double mu;
    mu = std::accumulate(node.begin(), node.end(), 0.0)/node.size();

    NumericMatrix ndata; // get data for node members to model
    for(int j = 0; j < members.size(); j++)
    {
      ndata.row(j) = data.row(members[j] - 1); // subtract 1 from member value b/c cpp is 0-indexed UPDATE helper function bindist to return 0-indexed values.
    }

    NumericVector bin;

    if(mu > dthresh) // if mean dist within node is greater than threshold, redefine the node.
    {
      NumericMatrix ndef(growfact, ndata.ncol());
      for(int j = 0; j < data.ncol(); j++)
      {
        std::vector<double> markdat;
        NumericVector ncol = ndata(_,j);
        markdat = std::for_each(ncol.begin(), ncol.end(), markdat);
        if(node.size() > minpts) // if node membership is high enough, model the distribution and redefine using probability weighted sample
        {
          std::vector<double> output;
          output = markdist(markdat, growfact);
          for(int n = 0; n < output.size(); n++)
          {
            ndef(n,j) = output[n];
          }
        }
        else // if membership is low sample directly from data w/o modelling
        {
          std::vector<double> output;
          output = resamp(markdat, growfact);
          for(int n = 0; n < output.size(); n++)
          {
            ndef(n,j) = output[n];
          }
        }
      }
      for(int n = 0; n < ndef.nrow(); n++)
        {
          bin.push_back(ndef(n,_));
        }
    }

    else // if dist within node is safe, keep it.
    {
      if(node.size() > minpts) // if node membership is high enough, try to improve node def by modelling
      {
        std::vector<double> ndef;
        for(int j = 0; j < ndata.ncol(); j++)
        {
          std::vector<double> markdat;
          NumericVector ncol = ndata(_,j);
          markdat = std::for_each(ncol.begin(), ncol.end(), markdat);
          std::vector<double> output;
          output = markdist(markdat, 1);

          ndef.push_back(output);
        }

        bin.push_back(ndef);
      }
      else // or just keep the same node def.
      {
        bin.push_back(oldbins(i,_));
      }
    }

    // attach new node def to new bin.mat result
    binmat.push_back(bin);

  }

  NumericMatrix binmatrix(binmat.size()/data.ncol(),data.ncol());
  int n = binmatrix.ncol() * binmatrix.nrow();
  for(int j = 0; j < n; j++)
  {
    binmatrix[j] = binmat[j];
  }

  return(binmatrix);
}


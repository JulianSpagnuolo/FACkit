#include <Rcpp.h>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/median.hpp>
#include <boost/accumulators/statistics/p_square_cumul_dist.hpp>
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

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector snlocation(std::vector<double> markdat)
{
  
  // create boost accumulator to calc variables used to model the distribution of marker j expression in node.
  accumulator_set< double, features<tag::median(with_p_square_cumulative_distribution),
                                    tag::count, tag::sum, tag::variance, tag::mean(immediate),
                                    tag::skewness, tag::moment<2>, tag::moment<3> > > nacc( p_square_cumulative_distribution_num_cells = 500 );

  
// use this if you want to get location directly using the median extractor: (simpler more direct method without using distributions).  
//  accumulator_set< double, features<tag::median(with_p_square_cumulative_distribution),
//                                    tag::count > > nacc( p_square_cumulative_distribution_num_cells = 500 );
  
  
  nacc = std::for_each(markdat.begin(), markdat.end(), nacc); // add the data to the accumulator
  // create convenient extractors of necessary data.
  boost::accumulators::extractor<boost::accumulators::tag::median> median_;
  boost::accumulators::extractor<boost::accumulators::tag::variance> stdev_;
  boost::accumulators::extractor<boost::accumulators::tag::skewness> skew_;
  
  auto markdist = boost::math::skew_normal_distribution<> (median_(nacc), std::sqrt(stdev_(nacc)), skew_(nacc));
  
  NumericVector loc;
  
  loc = markdist.location();

  // use this if you want to get location directly using the median extractor: (simpler more direct method without using distributions).  
  //loc = median_(nacc);
  
  return loc;
}



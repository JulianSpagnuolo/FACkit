#include <Rcpp.h>
using namespace Rcpp;

//' @export

// [[Rcpp::export]]
std::vector< Rcpp::CharacterVector > binclust2(IntegerMatrix binmat, CharacterVector rowids) {

  CharacterVector binned;
  std::vector< Rcpp::CharacterVector > binlist;


  for(int i = 0; i < binmat.nrow(); i++)
    {
    IntegerVector row_x = binmat.row(i);
    std::vector<int> x = as<std::vector<int> >(row_x);

    // if cell has already been binned, skip, else create new bin and find matching cells
    if(std::find(binned.begin(), binned.end(), rowids(i)) != binned.end()){
      // figure out better way to do this ... ie. opposite of this if statement such that if true, run the rest, else nothing.
    }else{
      binned.push_back(rowids(i)); // add the cell to the list of cells already binned

      CharacterVector current_bin; // create a new bin for matching cells
      current_bin.push_back(rowids(i)); // add the first cell at the start

      // find matching cells
      for(int n = 0; n < binmat.nrow(); n++)
      {
        if(n != i) // skip matching self
        {
          IntegerVector row_y = binmat.row(n);
          std::vector<int> y = as<std::vector<int> >(row_y);
          if(x == y) {
            // add the matching cell to the current bin and the already binned vector
            current_bin.push_back(rowids(n));
            binned.push_back(rowids(n));
          }
        }
      }
      // Add the completed bin to the binlist
      binlist.push_back(current_bin);
    }
    }
  return(binlist);
}

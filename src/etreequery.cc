/*
The Evolving Tree
Copyright (C) 2004 Jussi Pakkanen, Petri Turkulainen

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/

/*! \file
  \brief This file contains the query program of the Evolving Tree.
*/

#include<iostream>
#include"etree.hh"

#define OUTPUT_SEPARATOR "-----\n"

//! Gets classification result for a given tree and data.
/*!
 Usage 'etreequery \<training data file> \<testing data file> \<tree file>'

 Data files must be in etree format.

 Output on success:

 first test vector data label (usually a class id)<br>
 vector labels that were mapped to the BMU (may be empty), one per line<br>
 -----<br>
 second test vector data label<br>
 vector labels that were mapped to the BMU (may be empty), one per line<br>
 -----<br>
 ...<br>

 In case of an error, some kind of an error message is printed.
*/

int main(int argc, char **argv) {
  ifstream treefile;
  vector<int> result;

  if(argc != 4) {
    cout << "Usage:\n" << argv[0] << " <training data file> <testing data file> <tree file>\n";
    exit(-1);
  }

  DataCache train, test;
  if(!train.ReadFromDisk(argv[1])) {
    cout << "Error: could not read training data file " << argv[1] << "\n";
    return -1;
  }

  if(!test.ReadFromDisk(argv[2])) {
    cout << "Error: could not read test data file " << argv[2] << "\n";
    return -1;
  }

  Etree etree(&train);

  treefile.open(argv[3]);
  if(treefile.bad()) {
    cout << "Error: could not open tree file " << argv[3] << "\n";
    return -1;
  }
  treefile >> etree;
  treefile.close();

  // Create the classification index and use it for classification.
  etree.CreateIndex();

  // First the query label.
  for(int k=0; k<test.GetSize(); k++) {
    cout << test.GetLabel(k) << "\n";
    
    // Then the best matches.
    result = etree.GetClassification(test.GetVectorNumber(k));
    for(unsigned int j=0; j<result.size(); j++)
      cout << train.GetLabel(result[j]) << "\n";

    // Print the separator
    cout << OUTPUT_SEPARATOR;
  }

  return 0;
}














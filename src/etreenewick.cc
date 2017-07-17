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
  \brief This file prints the tree and the data it indexes in Newick format.
*/

#include<iostream>
#include<iomanip>
#include"etree.hh"

//! Print Newick tree for a given tree and data.
/*!
 Usage 'etreenewick \<data file> \<tree file> \<output file>'

 Input data files must be in etree format.

 The output format can be found here: 
 http://evolution.genetics.washington.edu/phylip/newicktree.html

 In case of an error, an unspecified error message is printed.
*/

int main(int argc, char **argv) {
  ifstream treefile;
  ofstream outfile;
  vector<int> result;

  if(argc != 4) {
    cout << "Usage:\n" << argv[0] << " <data file> <tree file> <output file>\n";
    exit(-1);
  }

  DataCache data;
  if(!data.ReadFromDisk(argv[1])) {
    cout << "Error: could not read data file " << argv[1] << "\n";
    return -1;
  }

  Etree etree(&data);

  treefile.open(argv[2]);
  if(treefile.bad()) {
    cout << "Error: could not open tree file " << argv[2] << "\n";
    return -1;
  }
  treefile >> etree;
  treefile.close();

  outfile.open(argv[3], ios::trunc | ios::binary);
  if(outfile.bad()) {
    cout << "Error: could not open output file " << argv[2] << "\n";
    return -1;
  }
  outfile << setprecision(3);

  // Prepare the tree for printing.
  etree.CreateIndex();

  // Print out the result.
  etree.PrintNewickTree(outfile);

  outfile.close();
  return 0;
}














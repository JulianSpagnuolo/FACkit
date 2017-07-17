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
  \brief This file contains the training program of the Evolving Tree.
*/

#include<iostream>
#include"etree.hh"

//! Trains an evolving tree with the given parameters and saves it to disk.
/*!
 Usage 'etreetrain \<training data file> \<output file> \<parameter file>'

 Training data must be in etree format.

 Parameter file has on line with the following structure:

 <i>divisionthreshold divisionfactor eta0 sigma0 tau1 tau2 bmudecay k-means rounds</i>

 Output:

 If everything went ok, print "Success".

 Otherwise some kind of an error message is printed.
*/



void etreetrain(DataCache *data, double eta0, double sigma0, double tau1, double tau2, int divisionThreshold, int divisionFactor, int kmeanscount, double bdecay) {
  ofstream treefile;
//  ifstream params;

  srandom(time(NULL));

/*
  if(argc != 4) {
    cout << "Usage:\n" << argv[0] << " <training data file> "
	 <<"<output file> <parameter file>\n";
    exit(-1);
  }

  DataCache data;
  if(!data.ReadFromDisk(argv[1])) {
    cout << "Error: Unable to read data file " << argv[1] << "\n";
    return -1;
  }
*/

  Etree etree(&data);

/*
  params.open(argv[3]);
  if(params.bad()) {
    cout << "Error: Unable to read parameter file " << argv[3] << "\n";
    return -1;
  }

  params >> divisionThreshold >> divisionFactor >>
    eta0 >> sigma0 >> tau1 >> tau2 >> bdecay >> kmeanscount;
  params.close();
*/

  etree.SetParameters(divisionThreshold, divisionFactor, eta0,
		      sigma0, tau1, tau2, bdecay);

/*
  // Open tree output file before training.
  treefile.open(argv[2]);
  if(treefile.bad()) {
    cout << "Error: Unable to open tree output file " << argv[2] << "\n";
    return(-1);
  }
*/

  etree.Initialize();

  while (etree.HasTreeGrown()) {
    etree.SingleRound();
  }

  // Fine-tune with K-means.
  for(int i=0; i<kmeanscount; i++)
    etree.KmeansAdjust();

//  treefile << etree;
//  treefile.close();

  return etree;
}



















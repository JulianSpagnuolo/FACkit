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
  \brief  An implementation of DataCache.

  All the functions needed to have a working data cache.
*/

#include"datacache.hh"

//! Default constructor, sets everything to zero.

DataCache::DataCache() {
  dim = 0;
  shuffler = NULL;
  currentVector = 0;
}

//! Deallocates pointer data.
DataCache::~DataCache() {
  DeleteData();
}

//! Totally obliterates all traces of the data space and shuffler.

void DataCache::DeleteData() {
  unsigned int i;
  for(i=0; i<data.size(); i++)
    delete[] data[i];
  data.clear();
  delete []shuffler;
  shuffler = NULL;
  labels.clear();
}

//! Reads an etree data file from disk..
/*!  A data vector's label will be the last element on its line.
  \return true on success.
*/

bool DataCache::ReadFromDisk(string filename) {
  ifstream ifile;
  char line[1024];
  string label;
  streampos fsize, thisline;
  double *newvec;

  DeleteData();
  ifile.open(filename.c_str(), ios::in);
  if(ifile.bad())
    return false;

  ifile.seekg(0, ios::end);
  fsize = ifile.tellg();
  ifile.seekg(0, ios::beg);

  // Skip the header, meaning the lines that start with "#".
  do {
    thisline = ifile.tellg();
    ifile.getline(line, 1024);
  } while(line[0] == '#');
  ifile.seekg(thisline, ios::beg);

  // Now read in data.
  ifile >> dim;
  if(ifile.bad() || dim==0) {
    ifile.close();
    return false;
  }
  ifile.ignore(2, '\n');

  //  cout << "Dimensio: " << dim << ".\n";

  // If some cases tellg() may return value of -1. Be prepared for it.
  while(ifile.tellg() < fsize && ifile.tellg() >= 0) {
    newvec = new double[dim];
    for(int i=0; i<dim; i++)
      ifile >> newvec[i];
    if(ifile.bad()) {
      ifile.close();
      return false;
    }
    data.push_back(newvec); // Deleting memory is delegated.
    ifile >> label;
    labels.push_back(label);
    ifile.ignore(2, '\n');
  }

  ifile.close();
  CreateShuffler();
  return true;
}



//! Creates the shuffler data structure. 

/*! 
  Note that this shuffler returns the data in the original order. To actually randomize the
  data, call Shuffle().
*/

void DataCache::CreateShuffler() {
  if(shuffler)
    delete []shuffler;
  shuffler = new int[data.size()];
  for(unsigned int i=0; i<data.size(); i++)
    shuffler[i] = i;
}

//! Return a random vector from the data set.

const double* DataCache::GetRandomVector() {
  return data[random()%data.size()];
}

//! Randomizes the order of the vectors.

void DataCache::Shuffle() {
  unsigned int to, from;
  int temp;
  for(from=0; from<data.size(); from++) {
    to = random() % data.size();
    temp = shuffler[to];
    shuffler[to] = shuffler[from];
    shuffler[from] = temp;
  }
}

//! Reset the shuffler to its beginning state.

void DataCache::StartNewRound() {
  currentVector = 0;
}

//! Returns the next vector in the cache.
/*
  \param last set to true if the current data vector is the last one in the cache.
  \return the data vector, or NULL if out of bounds.
*/

const double* DataCache::GetNextVector(bool &last) {
  // First check for index out of bounds;
  if(currentVector < 0 || currentVector >= (int) data.size()){
    last = true;
    return NULL;
  }

  if(currentVector == (int)data.size()-1)
    last = true;
  else
    last = false;

  return data[shuffler[currentVector++]];
}

//! Get a vector with the given number.

const double* DataCache::GetVectorNumber(int number) const {
  return data[number];
}

//! Find the min and max values of every vector element.
/*!
 Returns an array of size dim*2. The first dim elements have the
 minimum values of each element and the second dim elements have the
 max values.

 It is the caller's responsibility to free[] the array.
*/

double* DataCache::GetMinMax() const {
  double *minmax = new double[2*dim];
  int i;
  unsigned int j;

  for(i=0; i<dim; i++) {
    minmax[i] = data[0][i];
    minmax[i+dim] = data[0][i];
  }

  for(j=1; j<data.size(); j++) 
    for(i=0; i<dim; i++) {
      if(data[j][i] < minmax[i])
	minmax[i] = data[j][i];
      if(data[j][i] > minmax[i+dim])
	minmax[i+dim] = data[j][i];
    }

  return minmax; 
}

//! Return a textual label of the given element.
/*! The label is usually a class label, but can be anything.

\param n Index to the desired vector.
\return The class label or "INDEXOUTOFBOUNDS", if, surprisingly, the index is out of bounds.
*/

string DataCache::GetLabel(int n) const {
  string foo;
  if(n < 0 || n >= (int)labels.size())
    return "INDEXOUTOFBOUNDS";
  foo = labels[n];
  return foo;
}

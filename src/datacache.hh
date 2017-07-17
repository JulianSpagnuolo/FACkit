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
  \brief  A container for data vectors.

  This file defines a data vector container and accessor expected by 
  the Evolving Tree. 
*/

#ifndef _datacache_hh_
#define _datacache_hh_

#include<string>
#include<vector>
#include<fstream>

#include"etreedef.hh"

//! A contained and accessor for data vectors used by the Evolving Tree.
/*! This particular implementation reads SomPAK-type data files, and
  keeps them in memory. For more demanding applications (e.g.
  accessing data stored in SQL) use subclassing.
*/

class DataCache {

private:

protected:

  vector<double*> data; //!< The data vectors. Remember that vector does not do autodisallocation.
  int dim;              //!< Data vectors' dimension.

  int *shuffler;     //!< Used by GetNextVector to pass through data randomly.
  int currentVector; //!< Which vector are we processing at the moment.

  void DeleteData();
  void CreateShuffler();

  vector<string> labels; //!< Class labels for the data vectors.

public:

  DataCache();
  ~DataCache();
  bool ReadFromDisk(string filename);

  // Accessor functions.
  const double* GetRandomVector();
  void Shuffle();
  void StartNewRound();
  const double* GetNextVector(bool &last);
  const double* GetVectorNumber(int number) const;
  int GetDimension() const { return dim; }
  int GetSize() const { return data.size(); }
  double *GetMinMax() const;
  string GetLabel(int n) const;
};

#endif




















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
    \brief This file defines the interface of a single node.

 This file defines a single non-intelligent data node. It knows only
 its own prototype, neighbors and how to update itself towards a given
 data vector. Higher-level decisions are made outside the node.
*/

#ifndef _node_hh_
#define _node_hh_

#include<iostream>
#include<fstream>
#include<vector>
#include<math.h>

#include"etreedef.hh"

/*! \brief A node used in the evolving tree.

    A simple node that forms a basis of the Evolving Tree. It only
    contains a prototype vector and knowledge of its parents and
    children.
*/

class Node {

private:

  double *proto; //!< The prototype vector.
  int dim;       //!< Dimension of the prototype vector.

  int mynumber;         //!< This node's number (or index) in the tree.
  int parent;           //!< The parent's number in the tree.
  vector<int> children; //!< Children's numbers, if any.

  double bestmatches; //!< How many times this node has been the BMU.

protected:

public:

  Node();
  Node(const Node &source);
  ~Node();
  double Distance(const double *sample) const;
  double Distance(const Node &other) const { return Distance(other.proto); }
  int GetParent() const { return parent; }
  void SetParent(const int newparent) { parent = newparent; }
  const vector<int>& GetChildren() const { return children; }
  void SetChildren(const vector<int> newchildren) { children = newchildren; }
  void AddChild(int newchild) { children.push_back(newchild); }
  void ClearChildren() { children.clear(); }
  int GetNumber() const { return mynumber; }
  void SetNumber(int newnumber) { mynumber = newnumber; }
  int GetDimension() { return dim; }

  const double* GetDataVector() const;
  void SetDataVector(int ndim, const double* newvector);
  void SetDataVector(const double* newvector);
  void UpdateTowards(const double* targetvector, const double weight);

  double IncrementBMU();
  double BestMatches() const { return bestmatches; }
  void ClearBestMatches() { bestmatches = 0; }
  void SetBestMatches(double bm) { bestmatches = bm; } 

  void CopyDataFrom(const Node &source);
  void PrintInTextForm(ostream &os);
  Node& operator=(const Node &source);
  friend ostream& operator<<(ostream &os, Node &outnode);
  friend istream& operator>>(istream &is, Node &innode);

  // Functions that print the state of the node in human interpretable
  // form.
  void PrintState(ostream &os) const;
  void PrintProto(ostream &os) const;
};

#endif






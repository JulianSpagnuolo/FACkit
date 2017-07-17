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


/*! 
  \file
  \brief The implementation of a node.

  This file contains the implementation of the nodes used in the Evolving Tree.
*/

#include"node.hh"

//! Basic constructor. 
/*!
  Sets everything to zero.
*/

Node::Node() {
  proto = NULL;
  dim = 0;
  mynumber = 0;
  parent = 0;
  bestmatches = 0;
}

//! Copy constructor.
/*! 
  \param source A node to be copied.
*/

Node::Node(const Node &source) {
  proto = NULL;
  CopyDataFrom(source);
}

//! Destructor.
/*!
  Very basic, just unallocates the data vector.
*/

Node::~Node() {
  if(proto)
    delete[] proto;
}


//! Euclidean distance between prototype and sample.
/*!
  Calculates the euclidian distance between the node and a given data vector.
  \param sample A data vector, must have the same dimension as the node. (not checked)
  \return A double containing the distance.
*/

double Node::Distance(const double *sample) const {
  double result = 0;

  for(int i=0; i<dim; i++)
    result += pow(sample[i] - proto[i], 2);

  return sqrt(result);
}

//! Get the data vector.
/*! 
  Returns a copy of the node's data vector.
  \return The data vector.
*/

const double* Node::GetDataVector() const {
  return proto;
}

//! Assign a new data vector.
/*!
  Replaces the node's data vector with a new one with the given dimension.
  \param ndim The new dimension.
  \param newvector The new data vector of size <i>ndim</i>.
*/

void Node::SetDataVector(int ndim, const double* newvector) {
  if(proto != NULL)
    delete[] proto;
  dim = ndim;

  proto = new double[dim];
  SetDataVector(newvector);
}

//! Assign a new data vector.
/*!
  Sets the node's data vector to a new one.
  \param newvector The new vector, assumed to have the same size as the current vector.
*/

void Node::SetDataVector(const double* newvector) {
  memcpy(proto, newvector, dim*sizeof(double));
}

//! Move towards the given vector.
/*! 
  Updates the node towards the given data vector by the amount
  specified by weight. Currently implemented with the Kohonen rule.
  
  \param targetvector The vector to update towards.
  \param weight How much to update, should be between 0 and 1.
*/

void Node::UpdateTowards(const double* targetvector, const double weight) {
    for(int j=0; j<dim; j++)
      proto[j] += weight*(targetvector[j] - proto[j]);
}

//! Increment BMU count.
/*!
  Call this when the node has been selected as the BMU
  \return BMU count thus far, including the current one.
*/

double Node::IncrementBMU() {
  return ++bestmatches;
}

//! Duplicate a state.
/*!
   Duplicate the state of some other node to this one. Also copies the
   count numbers, so be careful not to shoot yourself in the foot with
   those.
   \param source The source to copy data from.
*/

void Node::CopyDataFrom(const Node &source) {
  if(proto)
    delete[] proto;

  dim = source.dim;
  if(source.proto == NULL)
    proto = NULL;
  else {
    proto = new double[dim];
    memcpy(proto, source.proto, dim*sizeof(double));
  }

  mynumber = source.mynumber;
  parent = source.parent;
  children = source.children;

  bestmatches = source.bestmatches;
}

//! Assignment operator
/*
  Copies the state from a given node. See warnings in CopyDataFrom.
  \param source The source of the data.
*/

Node& Node::operator=(const Node &source) {
  CopyDataFrom(source);
  return *this;
}

//! Serializes the node.
/*! Writes the node to the given stream as a binary dump. The dump
  does not have any magic numbers or other identifiers.

  \param os The stream to write to (usually a file).
  \param outnode The node to serialize.
  \return The stream written to.
*/

ostream& operator<<(ostream &os, Node &outnode) {
  //  int i;
  unsigned int s = outnode.children.size();
  os.write((char*)(&(outnode.dim)), sizeof(int));
  os.write((char*)(outnode.proto), outnode.dim*sizeof(double));
  os.write((char*)(&(outnode.mynumber)), sizeof(int));
  os.write((char*)(&(outnode.parent)), sizeof(int));
  os.write((char*)(&s), sizeof(unsigned int));
  os.write((char*)(&(outnode.children.front())), s*sizeof(int));
  os.write((char*)(&(outnode.bestmatches)), sizeof(int));
  return os;
}

//! Creates a new node from binary stream.
/*! Reads the binary stream and uses it to build a new node. No error
  checking of any kind is done, so make sure that the stream points to
  a sensible location.

  \param is The stream to read (usually a file).
  \param innode The node to write the information to.
  \return The stream just read.
*/

istream& operator>>(istream &is, Node &innode) {
  int i;
  unsigned int children, j;
  if(innode.proto)
    delete[] innode.proto;
  innode.children.clear();

  is.read((char*)(&(innode.dim)), sizeof(int));
  // Allocate memory.
  if(innode.dim == 0)
    innode.proto = NULL;
  else
    innode.proto = new double[innode.dim];

  is.read((char*)(innode.proto), innode.dim*sizeof(double));
  is.read((char*)(&(innode.mynumber)), sizeof(int));
  is.read((char*)(&(innode.parent)), sizeof(int));
  is.read((char*)(&children), sizeof(unsigned int));

  for(j=0; j<children; j++) {
    is.read((char*)(&i), sizeof(int));
    innode.children.push_back(i);
  }
  is.read((char*)(&(innode.bestmatches)), sizeof(int));

  return is;
}


//! Prints the data vector.
/*! Prints the elements of the location vector to the stream. The
elements are separated by a space and the line is ended by a linefeed.
This function can be used, for example, exporting the nodes as Matlab
matrices.

\param os The stream to write to.
*/

void Node::PrintInTextForm(ostream &os) {
  for(int i=0; i<dim; i++)
    os << proto[i] << " ";
  os << "\n";
}

//! Print state info.
/*!
  Outputs the state of the node in plain text to the specified
  stream. Does not output prototype vector.

  \param os The stream to write to.
*/

void Node::PrintState(ostream &os) const {
  os << "I am node number " << mynumber << ".\n";
  if(parent == ROOT_LOOP) {
    os << "I have no parent node.\n";
  } else {
    os << "My parent is node number " << parent << ".\n";
  }
  os << "I have been the best matching unit " << bestmatches << " times.\n";
  if(children.size() == 0)
    os << "I have no children.\n";
  else {
    if(children.size() == 1) 
      os << "My child is node ";
    else 
      os << "My children are nodes ";
    for(unsigned int i=0; i<children.size()-1; i++)
      os << children[i] << ", ";
    os << children.back() << ".\n";
  }

  os << "\n"; // Leave an empty line for separation.
}

//! A verbose description of the data vector.
/*!
  Prints a free form description of the data vector.

  \param os The stream to write to.
*/

void Node::PrintProto(ostream &os) const {
  if(proto == NULL) {
    os << "The prototype vector of node " << mynumber << " is non-existant.\n\n";
    return;
  }
  // Here dimension is guaranteed to be at least 1.

  os << "The prototype vector of node " << mynumber << " has " << dim << " ";
  if(dim == 1) {
    os << "element, which is " << *proto << ".\n";
  } else {
    os << "elements, which are ";
    for(int i=0; i<dim-1; i++)
      os << proto[i] << ", ";
    os << proto[dim-1] << ".\n";
  }
}

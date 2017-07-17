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
  \brief The implementation of the Evolving Tree.

  This file contains all the functions that make the Evolving Tree tick.
*/

#include<cfloat>
#include<iostream>
#include<list>

#include"etree.hh"

// First we have some helper stuff that are required by the query modules.
static const double *cur_query = NULL;
static const vector<Node> *cur_nodes = NULL;
static bool distsorter(pair<int, int> x, pair<int, int> y);

//! A helper function for the old query method. Should not be used anymore.
static bool distsorter(pair<int, int> x, pair<int, int> y) {

  // Check tree distance first.
  if(x.first != y.first)
    return x.first < y.first;

  // Tree distance is the same, sort by data space distance
  return (*cur_nodes)[x.second].Distance(cur_query) < 
    (*cur_nodes)[y.second].Distance(cur_query);
}

//! A basic constructor.
/*!
  Sets all values to sensible guesstimates. For better performance, 
  tweak them by hand.

  \param dc A data cache holding the training data.
*/
Etree::Etree(DataCache *dc) : datavecs(dc) {
  round = 0;
  dim = datavecs->GetDimension();

  SetParameters(50, 7, 0.6, 0.7, 4.0, 4.0, 0.9);

  weightThreshold = 0.00001;
  growingThreshold = 0.05;
  
  prevTreeSize = -1;
}

//! Initialization
/*!
  Initializes the tree and its root node.
*/
void Etree::Initialize() {
  Node root;

  // Initialize the tree with one vector at the middle of the data cloud.
  double *minmax;
  double *center;
  
  minmax = datavecs->GetMinMax();
  center = new double[dim];
  for(int i=0; i<dim; i++)
    center[i] = (minmax[i] + minmax[i+dim])/2;
  root.SetDataVector(dim, center);
  root.SetParent(ROOT_LOOP);
  root.SetNumber(0);
  nodes.push_back(root);
  fringe.push_back(0);
    
  delete []center;
  delete []minmax;
}

//! Sets the training parameters to the given values.
/*!
*/
void Etree::SetParameters(int dT, int dF, double e0, double s0, double t1, 
			  double t2, double bdecay) {
  divisionThreshold= dT;
  divisionFactor = dF;
  eta0 = e0;
  sigma0 = s0;
  tau1 = t1;
  tau2 = t2;
  bmuDecay = bdecay;
}

//! Find the best matching unit.
/*!
 Starting at the root, progresses down along the path until a leaf
 node is found.

 Note that this assumes nodes[0] is the root. Being a leaf node is
 equivalent to having zero children.
*/
int Etree::FindBMU(const double *datavec) const {
  int curnode = 0;
  double curdist, bestdist;
  int bestind;

  // curnode points to the node currently examined. Check children,
  // find the one with the smallest distance, process that.
  while((nodes[curnode]).GetChildren().size() > 0) {
    const vector<int> &childs = nodes[curnode].GetChildren();
    unsigned int i;
    bestind = 0;

    bestdist = nodes[childs[bestind]].Distance(datavec);
    
    // Find the child with the smallest distance.
    for(i=1; i<childs.size(); i++) {
      curdist = nodes[childs[i]].Distance(datavec);
      if(curdist < bestdist) {
	bestdist = curdist;
	bestind = (int) i;
      }
    }
    curnode = childs[bestind];
  }

  return curnode;
}

//! Calculate the update weight.
/*! Figure out how much to update a node towards a data vector. Has
  exponential neighborhood and decay.

  \param distance The tree distance between the BMU and the node to be updated.
  \return An update factor, between 0 and 1.
*/
double Etree::CalculateUpdateFactor(const double distance) const {
  double sigma;
  double foo;
  sigma = sigma0 * exp(-(double)round/tau1);
  foo =  exp(-pow(distance, 2) / (2*pow(sigma, 2)));
  return foo;
}


static double weightbase;  // eta0*exp(-round/tau2)

//! Moves the fringe nodes towards the data vector.
/*!
  Updates the fringe node locations towards the data vector.

  \param BMU Index of the best matching unit
  \param datavec The data vector to update towards.
*/
void Etree::UpdateNodes(int BMU, const double *datavec) {
  weightbase = eta0*exp(-round/tau2);

  // update the BMU
  double weight = weightbase * CalculateUpdateFactor(0);
  nodes[BMU].UpdateTowards(datavec,weight);

  // update the neighbours
  if (nodes[BMU].GetParent()!=ROOT_LOOP)
    UpdateNode(nodes[BMU].GetParent(),BMU,BMU,datavec,1);
}

//! Moves the fringe nodes towards the data vector.
/*!
  Recursively goes through the tree and updates weights of
  the leaf nodes.

  \param curnode A current node being examined. 
  \param prevnode A previously examined node
  \param BMU Index of the best matching unit
  \param datavec The data vector to update towards.
  \param distance A number of hops between the current node and the BMU
*/
void Etree::UpdateNode(int curnode, int prevnode, int BMU, 
		       const double *datavec, double distance) {
  unsigned int i;

  double treedist = distance - 1;
  double weight = weightbase *
    CalculateUpdateFactor(treedist);
  
  if (weight < weightThreshold) return;

  // update leaf nodes only
  if (nodes[curnode].GetChildren().size()==0) {
    nodes[curnode].UpdateTowards(datavec,weight);
  } else {
    // update the children of the node
    const vector<int> &childs = nodes[curnode].GetChildren();
    for (i=0;i<childs.size();i++) {	
      if (childs[i]!=prevnode)
	UpdateNode(childs[i],curnode,BMU,datavec,distance+1);
    }

    // update the siblings of the node
    int parent = nodes[curnode].GetParent();
    if (parent!=ROOT_LOOP && parent!=prevnode)
      UpdateNode(parent,curnode,BMU,datavec,distance+1);	
  }	
  
}

//! Calculate the tree distance between nodes.
/*! Figuring out the tree distance is quite simple. First we find out
  the paths from the nodes to root. Then we find the place where they
  differ. Then we just add the lengths of the remaining paths and
  subtract one.

  \param node1 Index to the first node.
  \param node2 Index to the second node.
  \return The tree distance.
*/
int Etree::TreeDistance(const int node1, const int node2) const {
  // This is the only case the algorithm below does not work with. So
  // deal with it first.

  if(node1 == node2)
    return 0;

  vector<int> path1 = FindTreePath(node1);
  vector<int> path2 = FindTreePath(node2);
  unsigned int i = 0;

  while(i < path1.size() && i < path2.size() && path1[i] == path2[i])
    i++;
  
  // Make i point to the last common element of both paths.
  if(i!=0)
    i--;
  
  // Now calculate the distance. That is the sum of the depths of the
  // left and the right subtrees minus one.
  return path1.size()-i-1 + path2.size()-i-1 - 1;
}

//! Find path from root to the specified node.
/*!
 Returns a vector of integers that tell the path from the root node
 to the node specified along the tree. The first element is the root
 and the last element is the node in question.
 
 \param nodenum Index to the node to be calculated.
 \return The path.
*/
vector<int> Etree::FindTreePath(int nodenum) const {
  vector<int> path;
  int newparent;
  newparent = nodes[nodenum].GetParent();
  while(newparent != ROOT_LOOP) {
    path.push_back(nodenum);
    nodenum = newparent;
    newparent = nodes[newparent].GetParent();
  }
  path.push_back(nodenum);
  
  // The path is currently backwards. Flip it.
  return vector<int>(path.rbegin(), path.rend());
}

//! Trains the tree with one data vector.
/*!
 This takes a datavector and does one iteration of the training
 algorithm. It also splits the BMU should that be necessary. It does
 not update the round variable, because one round consists of many
 iterations.

 \param datavec The data vector to use in training.
*/
void Etree::SingleIteration(const double *datavec) {
  int bmu = FindBMU(datavec);
  int i;

  // With the BMU known, we can now update the fringe nodes.
  UpdateNodes(bmu, datavec);

  // Do we need to split the node?
  if(nodes[bmu].IncrementBMU() >= divisionThreshold) {
    for(i=0; i<divisionFactor; i++) {
      Node node = nodes[bmu];
      
      nodes.push_back(node); // The new nodes are duplicates
      nodes.back().SetParent(bmu); // of the old except for links.
      nodes.back().SetNumber(nodes.size()-1);
      nodes.back().ClearChildren();
      nodes.back().ClearBestMatches();
      
      // Now add the referencing housekeeping info.
      fringe.push_back(nodes.size()-1);
      nodes[bmu].AddChild(nodes.size()-1);
    }

    for(vector<int>::iterator it = fringe.begin();
	it != fringe.end(); it++) {
      if(*it == bmu) {
	fringe.erase(it);
      }
    }
  }
}

//! Get tree distances from the given node.
/*!
  This function calculates the tree distances between the given node
  and all fringe nodes. The result is sorted according to distance.

  The first element in pair<int, int> is the distance and the second
  one is the index of the node. This is so that STL's sort can be
  used directly.

  \param nodenum Index to the node
  \return The list of distances.
*/
vector<pair<int, int> > Etree::GetSortedFringeDistances(const int nodenum) 
  const {
  vector<pair<int, int> > distances;
  pair<int, int> adder;
  unsigned int i;

  for(i=0; i<fringe.size(); i++) {
    adder.second = fringe[i];
    adder.first = TreeDistance(nodenum, adder.second);
    distances.push_back(adder);
  }

  sort(distances.begin(), distances.end());
  return distances;
}

//! Sets a new data cache.
/*!
  \return true on success and false if tree is non-empty and data 
  vector sizes differ.
*/
bool Etree::SetDataCache(DataCache &newcache) {
  if(nodes.size() > 0 && newcache.GetDimension() != dim)
    return false;

  datavecs = &newcache; // We are not responsible for deleting the cache.
  return true;
}

//! Check for leafnodedness.
/*!
  \return true if the given node is a leaf node, false otherwise.
*/
bool Etree::IsLeafNode(unsigned int i) const {
  return nodes[i].GetChildren().size() == 0;
}

//! This function calculates one round of the training algorithm.
/*! In a nutshell: shuffle data vectors, run SingleIteration on them. Stop.

  \param exportDate If true, log every step to a file in Matlab form.
*/
void Etree::SingleRound(bool exportData) {
  bool last = false;
  ofstream ofile;
  if(datavecs->GetSize() == 0)
    return; // Training with no data is tricky.
  datavecs->Shuffle();
  datavecs->StartNewRound();

  if(exportData)
    ofile.open("flatvecs.ddt", ios::out | ios::app);

  // Go for it.
  while(!last) {
    const double *dv = datavecs->GetNextVector(last);
    SingleIteration(dv);

    if(exportData) 
      ExportShapeAsMatlab(ofile);
  }

  if(exportData)
    ofile.close();

  DecayBMUcounts();

  round++;
}

//! Export the current tree as a Matlab matrix.
/*!
 The matrix is built from rows of node prototype vectors. 
 Only exports the nodes that are in the fringe.

 \param varname the name of the Matlab matrix to write.
*/
void Etree::ExportChildrenAsMatlab(string varname, ostream &os) {
  unsigned int i;

  os << varname << " = [\n";
  for(i=0; i<fringe.size(); i++)
    nodes[fringe[i]].PrintInTextForm(os);

  os << "];\n";
}

//! Prints every non-leaf node's relationship data as plain text.
/*!
  Uses the following format.

  numberofelements<br>
  numberofchildren<br>
  parent	  <br>
  child1	  <br>
  child2	  <br>
  child3	  <br>
  ...		  <br>
  parent	  <br>
  child1          <br>
  ...

  \param os The output stream.
*/
void Etree::ExportShapeAsMatlab(ostream &os) {
  unsigned int i, j;

  // How many elements in total.
  os << nodes.size() - fringe.size() << "\n";
  for(i=0; i<nodes.size(); i++) {
    if(nodes[i].GetChildren().size() == 0)
      continue;

    // First number of children and the parent node.
    os << nodes[i].GetChildren().size() << "\n";
    nodes[i].PrintInTextForm(os);
    for(j=0; j<nodes[i].GetChildren().size(); j++)
      nodes[nodes[i].GetChildren()[j]].PrintInTextForm(os);
  }
}

//! Serializes the tree to the stream.
/*! Note that datacache is not exported. Also no magic numbers or 
  checksums are used.
  
  The order is as follows. First the constants in the order they
  appear in the class declaration. Next comes the fringe, then the
  serialized nodes.

  \param os The output stream.
  \param outtree The tree to output
*/
ostream& operator<<(ostream &os, Etree &outtree) { 
  int out;
  unsigned int i;
  os.write((char*)(&(outtree.dim)), sizeof(int));
  os.write((char*)(&(outtree.round)), sizeof(int));
  os.write((char*)(&(outtree.eta0)), sizeof(double));
  os.write((char*)(&(outtree.sigma0)), sizeof(double));
  os.write((char*)(&(outtree.tau1)), sizeof(double));
  os.write((char*)(&(outtree.tau2)), sizeof(double));
  os.write((char*)(&(outtree.bmuDecay)), sizeof(double));

  i = outtree.fringe.size();
  os.write((char*)&i, sizeof(int));
  for(i=0; i<outtree.fringe.size(); i++) {
    out = outtree.fringe[i];
    os.write((char*)(&out), sizeof(int));
  }

  i = outtree.nodes.size();
  os.write((char*)&i, sizeof(int));
  for(i=0; i<outtree.nodes.size(); i++)
    os << outtree.nodes[i];
  return os;
}

//! Reads an Evolving Tree from a binary stream.
/*! No checking of any kind is done. 
  \param is The input stream.
  \param intree The tree to write the data to.
*/
istream& operator>>(istream &is, Etree &intree) {
  intree.Unserialize(is);
  return is;
}

//! Read tree data from given file.
/*!
  Datacache must be set to something proper beforehand.

  \return false on failure, meaning file not found or data dimensions disagree.
*/
bool Etree::ReadFromDisk(string fname) {
  ifstream ifile(fname.c_str());
  //  ifile.open(fname.c_str(), ios::in | ios::nocreate);
  if(!ifile)
    return false;
  Unserialize(ifile);
  if(nodes.size() > 0 && dim != datavecs->GetDimension())
    return false;
  return true;
}

//! Extract the tree structure from a binary mess.
void Etree::Unserialize(istream &is) {
  int in;
  unsigned int i, tot;
  Node node;
  nodes.clear();
  fringe.clear();
  labels.clear();
  is.read((char*)(&(dim)), sizeof(int));
  is.read((char*)(&(round)), sizeof(int));
  is.read((char*)(&(eta0)), sizeof(double));
  is.read((char*)(&(sigma0)), sizeof(double));
  is.read((char*)(&(tau1)), sizeof(double));
  is.read((char*)(&(tau2)), sizeof(double));
  is.read((char*)(&(bmuDecay)), sizeof(double));

  is.read((char*)&tot, sizeof(int));
  for(i=0; i<tot; i++) {
    is.read((char*)(&in), sizeof(int));
    fringe.push_back(in);
  }

  is.read((char*)&tot, sizeof(int));
  for(i=0; i<tot; i++) {
    is >> node;
    nodes.push_back(node);
  }
}

//! Returns the nodes that are at the specified depth in the tree.
/*!  Depth of 0 means root, 1 means the first level and so on.

  \return The node list.
*/
vector<int> Etree::NodesAtDepth(int depth) {
  vector<int> res;
  for(unsigned int i=0; i<nodes.size(); i++)
    if(FindTreePath(i).size() == (unsigned int) depth)
      res.push_back(i);
  return res;
}

//! Returns the index of the deepest node. Ties broken arbitrarily.
unsigned int Etree::MaxDepthNode() const {
  unsigned int maxind = 0;
  int maxdepth = FindTreePath(maxind).size();
  for(unsigned int i=1; i<nodes.size(); i++)
    if(FindTreePath(i).size() > (unsigned int) maxdepth)
      maxind = i;
  return maxind;
}

//! Return the path to the deepest node.
int Etree::GetMaxDepth() const {
  unsigned int deepest = MaxDepthNode();
  return FindTreePath(deepest).size();
}

//! Returns vector indicating how many nodes are at each level of the
//! tree.

vector<int> Etree::GetDepthDistribution() const {
  unsigned int maxdepthnode, maxdepth, i;
  vector<int> res;
  maxdepthnode = MaxDepthNode();
  maxdepth = FindTreePath(maxdepthnode).size();
  for(i=0; i<=maxdepth; i++)
    res.push_back(0);

  for(i=0; i<nodes.size(); i++)
    res[FindTreePath(i).size()]++;

  return res;
}

//! Prints some statistics in text form.
void Etree::PrintTreeInfo(ostream &os) const {
  vector<int> ddist = GetDepthDistribution();
  os << "Thise are the statistics for the current evolving tree.\n";
  os << "Data dimension is " << dim << ".\n";
  os << "Current round is number " << round << ".\n";
  os << "The training parameters are:\n";
  os << "Eta0 " << eta0 << ", sigma0 " << sigma0 << "\n";
  os << "Tau1 " << tau1 << ", tau2 " << tau2 << "\n";
  os << "DivisionThreshold " << divisionThreshold
     << ", divisionFactor " << divisionFactor << "\n";

  os << "Size of the tree is " << nodes.size() << ".\n";
  os << "Size of the fringe is " << fringe.size() << ".\n";

  if(ddist.size() != 0) {
    os << "The deepest node is at level " << GetMaxDepth() << ".\n";
    os << "The depth distribution is:\n";
    for(unsigned int i=0; i<ddist.size(); i++)
      os << i << ": " << ddist[i] << "\n";
  }
  os << "\n";
}

//! Print the relationship info of a given node in plain text.
void Etree::PrintNodeInfo(int nodenum, ostream &os) const {
  if(nodenum < 0 || nodenum >= (int)nodes.size())
    return;
  nodes[nodenum].PrintState(os);
}

//! Print the prototype vector of a given node as plaintext.
void Etree::PrintNodeProto(int nodenum, ostream &os) const {
  if(nodenum < 0 || nodenum >= (int)nodes.size())
    return;
  nodes[nodenum].PrintProto(os);
}

//! Builds the index information. 
/*! This should be called after training, but before querying.
 */
void Etree::CreateIndex() {
  int dcsize;
  int i;
  string label;

  ClearIndex();
  dcsize = datavecs->GetSize();

  // Create the string vectors to add.
  for(i=0; i<(int)nodes.size(); i++) {
    vector<int> drEvil;
    labels.push_back(drEvil);
  }

  for(i=0; i<dcsize; i++) {
    int ind = FindBMU(datavecs->GetVectorNumber(i));
    labels[ind].push_back(i);
  }
}

//! Seeks the best matches to the given query vector.
/*! Finds the BMU, and returns the labels of the data vectors mapped to it.
  \param queryvec The data vector containing the query, assumed to have the 
  correct dimension.
  \return Index list of training vectors that have been mapped to the BMU.
*/
vector<int> Etree::GetClassification(const double *queryvec) const {
  vector<int> result;
  int BMU;
  vector<pair<int, int> > distances;

  BMU = FindBMU(queryvec);
  result = GetIndexLabels(BMU);

  if(result.size() > 0)
    return result;

  // BMU was empty, so seek the next closest nonempty fringe node.

  // Sort results according to:
  // A) Node distance
  // B) Data space distance to query vector
  cur_query = queryvec;
  cur_nodes = &nodes;
  sort(distances.begin(), distances.end(), distsorter);
  cur_query = NULL;
  cur_nodes = NULL;

  for(unsigned int i=0; i<distances.size(); i++) {
    result = GetIndexLabels(distances[i].second);
    if(result.size() > 0)
      return result;
  }

  result.clear();
  return result;
}

//! Prints the distribution of mapped data vectors on leaf nodes.
/*! 
   Note that you must call CreateIndex before calling this. 
*/
void Etree::PrintIndexStats(ostream &out) const {
  unsigned int i, overflow;
  vector<unsigned int> counts;
  const unsigned int maxcount = 10;
  overflow = 0;

  for(i=0; i<maxcount; i++)
    counts.push_back(0);
      
  for(i=0; i<fringe.size(); i++) {
    int nodenum = fringe[i];
    unsigned int count = labels[nodenum].size();
    if(count > maxcount)
      overflow++;
    else
      counts[count]++;
  }

  // Now print the results.
  out << "Distribution of data vectors in leaf nodes.\n";
  for(i=0; i<counts.size(); i++)
    out << i << ": " << counts[i] << "\n";
  out << "n: " << overflow << ".\n";
}

//! Decrease BMU counters
/*! 
  Multiplies all the BMU counters of the leaves by regularization factor.
*/
void Etree::DecayBMUcounts() {
  unsigned int i;

  for(i=0; i<fringe.size(); i++) {
    nodes[fringe[i]].SetBestMatches(nodes[fringe[i]].BestMatches() *
				    bmuDecay);
  }
}

//! Fine-tune leaf locations.
/*! 
  Simple k-means optimization. Map all nodes to leafs, and then
  place the node to their center of mass.
*/

void Etree::KmeansAdjust() {
  unsigned int i, j, k;
  double *newloc = new double[dim];
  CreateIndex();

  for(i=0; i<fringe.size(); i++) {
    vector<int> cloud = GetIndexLabels(fringe[i]);
    if(cloud.size() == 0)
      continue;
    
    // Initialize new location to zero
    for(k=0; k<(unsigned int)dim; k++)
      newloc[k] = 0.0;

    // Add locations of child nodes.
    for(j=0; j<cloud.size(); j++) {
      const double *curvec = datavecs->GetVectorNumber(cloud[j]);

      for(k=0; k<(unsigned int)dim; k++)
	newloc[k] += curvec[k];
    }

    // Divide by # of nodes to get center of mass.
    for(k=0; k<(unsigned int)dim; k++)
      newloc[k] /= double(dim);

    // Set new location.
    nodes[fringe[i]].SetDataVector(newloc);
  }

  ClearIndex();
  delete[] newloc;
}

//! Optimizes locations of trunk nodes.
/*! A bottom-up recursion that places trunk nodes to the center of
  mass of their children.

  \param nodenum The index to the node to adjust. Should point to root
  when called from outside this function.
*/

void Etree::TrunkReshape(int nodenum) {
  const vector<int> ch = nodes[nodenum].GetChildren();
  unsigned int i;
  int j;
  double *newloc;

  if(ch.size() == 0)
    return; // End recursion when at leaf nodes.

  for(i=0; i<ch.size(); i++)
    TrunkReshape(ch[i]);

  newloc = new double[dim];
  for(i=0; i<(unsigned int) dim; i++)
    newloc[i] = 0;

  // Add up child vectors.
  for(i=0; i<ch.size(); i++) {
    const double *curvec = nodes[ch[i]].GetDataVector();
    for(j=0; j<dim; j++)
      newloc[j] += curvec[j];
  }
  
  // Divide to get center of mass.
  for(j=0; j<dim; j++)
    newloc[j] /= double(ch.size());

  nodes[nodenum].SetDataVector(newloc);

  delete[] newloc;
}

//! Stopping criteria 
/*! Check if the tree has grown between rounds. Can be used
  to stop the training.
  \return true if tree has grown enough otherwise false
*/
bool Etree::HasTreeGrown() {
  bool result;

  if (prevTreeSize==-1) 
    result = true;
  else if ((double(nodes.size()) / double(prevTreeSize)) < 
	   (1+growingThreshold))
    result = false;
  else
    result = true;

  prevTreeSize = nodes.size();
  return result; 
}

//! Returns the size of the subtree
/*! 
  Returns a number of nodes belonging to the subtree starting from 
  the <code>node<\code>.

  \param node the root of the subtree
  \return the size of the subtree
*/
int Etree::SubtreeSize(int node) {
  const vector<int> &childs = nodes[node].GetChildren();
  int size = childs.size();

  if (size == 0) return 0;

  for (unsigned int i=0; i<childs.size(); i++)
    size += SubtreeSize(childs[i]);
  
  return size;
}


//! Print the tree and data in Newick format.
/*! 

  Goes recursively through the tree and data and outputs the result to
  the given stream in Newick format. You probably want to map your
  data to leaf nodes with CreateIndex before calling this function.

  \param out The stream to output to.
  \param curnode The node we are currently processing.
*/
void Etree::PrintNewickTree(ostream &out, const int curnodenum) const {
  const Node &curnode = nodes[curnodenum];
  const int curparent = curnode.GetParent();
  const vector<int> &curchildren = curnode.GetChildren();

  // Start by printing wacky fun parentheses.
  out << "(";

  // Are we at a leaf node?
  if(curchildren.size() == 0) {
    const vector<int> &curdata = labels[curnode.GetNumber()];
    const unsigned int numchildren = curdata.size();
    for(unsigned int i=0; i<numchildren; i++) {
      const int curvec = curdata[i];
      if(i != 0)
	out << ",";
      out << datavecs->GetLabel(curvec) << ":"
	  << curnode.Distance(datavecs->GetVectorNumber(curvec));
	}
  } else { // No, print children recursively.
    const vector<int> &curchildren = curnode.GetChildren();
    const unsigned int numchildren = curchildren.size();
    for(unsigned int i=0; i<numchildren; i++) {
      const int curchild = curchildren[i];
      if(i != 0)
	out << ",";
      PrintNewickTree(out, curchild);
    }
  }

  // Wrap up by a closing parenthesis and the distance to parent node.
  out << ")";
  if(curparent != ROOT_LOOP)
    out << ":" << curnode.Distance(nodes[curparent].GetDataVector());

  // The tree ends with a semicolon.
  if(curnodenum == 0) 
    out << ";\n";
}

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

#ifndef _etree_hh_
#define _etree_hh_

#include<utility>
#include<algorithm>

#include"node.hh"
#include"datacache.hh"

/*! \brief This class defines the Evolving Tree.

  The Evolving Tree contains a group of nodes, the training parameters
  and the functions necessary to train and query the tree.
*/


class Etree {

private:

protected:

  vector<Node> nodes;  //!< A container for the nodes.
  vector<int> fringe;  //!< Indices to bottom layer nodes. 
                       //!< Only these may be BMUs.
  int dim;             //!< Dimension of data
  double round;        //!< How manyth iteration is going on.
  double eta0;         //!< Initial learning-rate decay value. 
  double sigma0;       //!< Initial neighborhood decay.
  double tau1;         //!< A time constant.
  double tau2;         //!< A time constant. See Haykin, pages 450-
  double bmuDecay;     //!< A regularization factor

  int divisionThreshold; //!< When a node is BMU this many times, it is split.
  int divisionFactor;    //!< When split, a node gets this many children.

  double weightThreshold; //!< A cutoff threshold for neighbor updating
  double growingThreshold; //!< A cutoff threshold for tree growing 
  int prevTreeSize; //!< The size of the tree in the previous round.

  DataCache *datavecs;  //!< A container for the training vectors.
  vector<vector<int> > labels; //!< A container for data vector indices that have been mapped to nodes.

  void Unserialize(istream &is);

  void UpdateNode(int curnode, int prevnode, int BMU,const double *datavec, 
		    double distance);
public:
  void Initialize();
  bool HasTreeGrown();

  friend class VisualizerData;

  Etree(DataCache *dc);

  void SetParameters(int dT, int dF, double e0, double s0, 
		     double t1, double t2, double bdecay);

  // Helper functions.
  int FindBMU(const double *datavec) const;
  double CalculateUpdateFactor(const double distance) const;
  void UpdateNodes(int BMU, const double *datavec);
  int TreeDistance(const int node1, const int node2) const;
  vector<int> FindTreePath(int nodenum) const;
  vector<pair<int, int> > GetSortedFringeDistances(const int nodenum) const;
  bool SetDataCache(DataCache &newcache);
  bool IsLeafNode(unsigned int i) const;

  //const vector<double *> MapData(int node);

  // Updating functions.
  void SingleIteration(const double *datavec);
  void SingleRound(bool exportData = false);
  void DecayBMUcounts();
  void KmeansAdjust(); 
  void TrunkReshape(int nodenum);

  void ExportChildrenAsMatlab(string varname, ostream &os);
  void ExportShapeAsMatlab(ostream &os);
  friend ostream& operator<<(ostream &os, Etree	&outtree);
  friend istream& operator>>(istream &is, Etree	&intree);
  bool ReadFromDisk(string fname);

  // Functions that print tree status info.
  int SubtreeSize(int node);
  int TotalNodes() { return nodes.size(); }
  int FringeSize() { return fringe.size(); }
  vector<int> NodesAtDepth(int depth);
  unsigned int MaxDepthNode() const;
  vector<int> GetDepthDistribution() const;
  int GetMaxDepth() const;
  void PrintNewickTree(ostream &out, const int curnodenum=0) const;

  // Statistics calculators and printers.
  void PrintTreeInfo(ostream &os) const;
  void PrintNodeInfo(int nodenum, ostream &os) const;
  void PrintNodeProto(int nodenum, ostream &os) const;

  // Functions dealing with the index.
  void CreateIndex();
  void ClearIndex() { labels.clear(); }
  vector<int> GetIndexLabels(int i) const { return labels[i]; }
  vector<int> GetClassification(const double *queryvec) const;
  void PrintIndexStats(ostream &out) const;

};

#endif














package Unused;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

import dataStructures.Node;
import dataStructures.Path;
import likelihood.PairwiseLikelihoodCoreStream2;

//This is the data structure that stores information about all of the individuals in the sample
//Essentially the nodes themselves encode the pedigree graph, and this class contains
//a sort of adjacency matrix, where instead of binary adjacencies, only the individuals are represented
//and the i,j^th entry represents the relationship between the i^th and j^th individuals
public class PedigreeOld {

	//pedigree variables
	public final int depth;
	public final int maxDepthForSamples;
	public final int numIndiv;
	public final double genTime;
	public final Node[] inds;
	public final Path[][] relationships;
	private final PairwiseLikelihoodCoreStream2 core;
	private final Random rGen;
	public Node[] nodes;
	private int nActiveNodes;
	public double logLikelihood; 
	

	
	//for age
	final double muGenTime = 29;
	final double varGenTime = 36;

	

	////// CONSTRUCTOR ///////
	//TODO handle known relationships
	public PedigreeOld(int depth, int maxDepthForSamples, Node[] inds, PairwiseLikelihoodCoreStream2 core, String marginalPath, String lkhdPath, Random rGen, int maxNumNodes, double genTime) throws IOException{
		
		this.numIndiv = inds.length;
		this.depth = depth;
		this.maxDepthForSamples = maxDepthForSamples;
		this.genTime = genTime;
		this.inds = inds;
		this.core = core;
		this.rGen = rGen;
		

		//initialize everyone to be unrelated
		this.relationships = new Path[numIndiv][numIndiv];
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				this.relationships[i][j] = new Path(0,0,0);
				this.relationships[j][i] = new Path(0,0,0);
			}
		}
		
	 
		//update nodes list
		this.nodes = new Node[maxNumNodes];;
			//add sampled nodes
			for(int i=0; i<inds.length; i++){
				nodes[i] = inds[i];
			}
			//add ghost nodes
			for(int i = inds.length; i<maxNumNodes; i++){
				nodes[i] = new Node(false, i);
			}
		
		
		this.nActiveNodes = inds.length;
		
		
		//initialize pairwise & marginal likelihoods
		core.setMarginals(marginalPath);
		core.setLikelihoods(lkhdPath);
		
		//compute current likelihood
		this.logLikelihood = 0d;
		for(Node i : inds){
			logLikelihood += core.getMarginal(i);
		}
		
		
		
	}
	
	
	////// SETTERS ////////
	public void setLogLikelihood(double lkhd){
		//this.logLikelihood = lkhd;
		logLikelihood = lkhd;
	}
	
	public void setNActiveNodes(int n){
		//this.nActiveNodes = n;
		this.nActiveNodes = n;
	}
	
	
	////// GETTERS ////////// 
	public double getLogLikelihood(){
		//return logLikelihood;
		return logLikelihood;
	}
	
	
	public Path[][] getRelationships(){
		return relationships;
	}
	
	public int getNActiveNodes(){
		return this.nActiveNodes;
	}
	
	public Node getIndividual(int index){
		return this.inds[index];
	}
	
	public Node getNode(int index){
		return this.nodes[index];
	}

	public Node getRandomNode(){ //works
		return this.nodes[rGen.nextInt(this.nActiveNodes)];
	}
	
	
	public Node[] getNRandomNodes(int sampleSize){//works

		Node[] toReturn = new Node[sampleSize];

		int nSampled = 0;
		int i = 0;
		int N = this.nActiveNodes;

		while(nSampled < sampleSize){
			
			if(rGen.nextDouble()*(N-i) < (sampleSize - nSampled)){
				toReturn[nSampled++] = this.nodes[i];
			}
			
			i++;
			
		}
		
		return toReturn;
	}
	


	
	public Node getRandomSampledNode(){
		return inds[rGen.nextInt(numIndiv)];
	}
	
	
	public Node getRandomPathAncestor(Node node, int targetDepth, int targetSex){
		
		if(node.getDepth() > targetDepth) 
			throw new RuntimeException("starting node has depth greater than target depth");

		Node currNode = node;
		
		//go up to target depth
		int sex;
		while(currNode.getDepth() < targetDepth){

			if(currNode.getDepth()==targetDepth-1){ //sex is deterined by targetSex
				sex = targetSex;
			}
			else{//randomly choose mom or dad
				sex = rGen.nextInt(2);
			}
			
			Node parent = currNode.getParentWithSex(sex);
			
			//if there was no parent with the given sex, make one
			if(parent==null){
				parent = makeNewNode(currNode.getDepth() + 1, sex);
				currNode.addParent(parent);
				parent.addChild(currNode);
			}
			
			currNode = parent;
					
		}

		if(currNode.getSex()!=targetSex) 
			throw new RuntimeException("node's sex at targetDepth does not match targetSex");
		
		return currNode;

		
	}
	
	
	
	////// CREATE/DELETE NODES ///////
	public Node makeNewNode(int depth, int sex){ //works

		//get new node from the end of array and increment pointer
		Node newNode = this.nodes[this.nActiveNodes++];
		newNode.setDepth(depth);
		newNode.setSex(sex);

		return newNode;
	}
	

	//return number of parent- and child-edges that were deleted
	public void deleteNode(Node nodeToDelete){ //works
	
		if(nodeToDelete.sampled)
			throw new RuntimeException("Trying to delete a sampled node");
		
		//clear outgoing edges and reset sex, depth
		nodeToDelete.reset();
		
		//update node list
		if(nodeToDelete.getIndex() != this.nActiveNodes-1){
			Node temp = this.nodes[this.nActiveNodes-1];
			this.nodes[nodeToDelete.getIndex()] = temp;
			this.nodes[temp.getIndex()] = nodeToDelete;
			temp.setIndex(nodeToDelete.getIndex());
			nodeToDelete.setIndex(this.nActiveNodes - 1);
		}
		

		
		this.nActiveNodes--;
		
	}
	
	
	
	
	////////// UPDATE PEDIGREE STRUCTURE////////////	
 	public void connect(Node parent, Node child){
 		
 		parent.addChild(child);
 		child.addParent(parent);
 		
 	}
	
	//delete unnecessary ghost nodes
	public Node clean(Node node){//works
		
		
		if(node==null) return null;
		
		
		//full sib case; here we don't care about the last existing nodes
		else if(!node.sampled && node.getNumEdges()==2 && node.getParents().size()==2){
			
			Node p1 = node.getParents().get(0);
			Node p2 = node.getParents().get(1);
			
			deleteNode(node);
			
			clean(p1);
			clean(p2);
			
			return null;
			
		}
		
		
		
		else if(node.sampled || node.getNumEdges() > 1){
			return node;
		}
		
		
		else{
			
			Node neighbor = null;
			
			if(node.getNumEdges()!=0){
				neighbor = node.getParents().size() > 0? node.getParents().get(0) : node.getChildren().get(0); 
			}
			
			deleteNode(node);
			
			//recurse on neighbor
			return clean(neighbor);
		}
		
		
	}
	
	
	
	//returns highest node touched
	public void cut(Node child, Node parent, boolean hasFullSib){ //works
		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeCut = parent.getConnectedNodes(new ArrayList<Node>());
		this.logLikelihood -= likelihoodLocalPedigree(nodesBeforeCut);
		
		
		// update graph structure
		child.removeParent(parent);
		parent.removeChild(child);
		
		
		for(Node ind : nodesBeforeCut){
			updateAdjMat(ind);
		}
		
		//update adjmat and likelihood
		if(hasFullSib){ //if child had a full sibling, the original pedigree is not split
			
			this.logLikelihood += likelihoodLocalPedigree(nodesBeforeCut);
			
		}

		else{ //split pedigrees

			clearVisit();
			List<Node> childPed = child.getConnectedNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> parentPed = parent.getConnectedNodes(new ArrayList<Node>());

			this.logLikelihood += likelihoodLocalPedigree(childPed);
			this.logLikelihood += likelihoodLocalPedigree(parentPed);
		}
		
	}
	
	
	//make a ghost copy and randomly assign children to the copy
	public void split(Node parent, Node splitParent, List<Node> splitChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit = parent.getConnectedNodes(new ArrayList<Node>());
		this.logLikelihood -= likelihoodLocalPedigree(nodesBeforeSplit);
		
		//make ghost parent
		splitParent.setChildren(splitChildren);
		
		//assign children to ghost parent
		for(Node i : splitChildren){
			i.removeParent(parent);
			i.addParent(splitParent);
			parent.removeChild(i);
		}
		
		
		for(Node ind : nodesBeforeSplit){
			updateAdjMat(ind);
		}
		
		
		//add new lkhd
		if(hasFullSib){//one pedigree
			
			this.logLikelihood += likelihoodLocalPedigree(nodesBeforeSplit);
			
		}
		else{
			
			clearVisit();
			List<Node> parentPed = parent.getConnectedNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> splitPed = splitParent.getConnectedNodes(new ArrayList<Node>());
			
			this.logLikelihood += likelihoodLocalPedigree(splitPed);
			this.logLikelihood += likelihoodLocalPedigree(parentPed);

		}
		       
		
		
	}
	

	
	//make a ghost copy and randomly assign children to the copy
	public void split2(Node parent, Node stayParent, List<Node> stayChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit =  parent.getConnectedNodes(new ArrayList<Node>());
		this.logLikelihood -= likelihoodLocalPedigree(nodesBeforeSplit);
		
		
		//grand parents
		for(Node i : parent.getParents()){
			i.addChild(stayParent);
			stayParent.addParent(i);
			i.removeChild(parent);
		}
		parent.getParents().clear();
		
		//children
		for(Node i : stayChildren){
			i.removeParent(parent);
			parent.removeChild(i);
			i.addParent(stayParent);
			stayParent.addChild(i);
		}
		
		for(Node ind : nodesBeforeSplit){
			updateAdjMat(ind);
		}
		
		//add new lkhd
		if(hasFullSib){//one pedigree
			
			this.logLikelihood += likelihoodLocalPedigree(nodesBeforeSplit);
			
		}
		else{
			clearVisit();
			List<Node> parentPed = parent.getConnectedNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> stayPed = stayParent.getConnectedNodes(new ArrayList<Node>());
			
			//setUnrelated(parentPed, stayPed);
			
			this.logLikelihood += likelihoodLocalPedigree(stayPed);
			this.logLikelihood += likelihoodLocalPedigree(parentPed);
			
		}
		
		
	}
	
	
	//recipient gets donor's children & parents
	//onePed==true if donor and recipient are already connected
	public void merge(Node donor, Node recipient, boolean onePed) {//works
						
		//subtract current subpedigree likelihoods
		clearVisit();
		List<Node> recipientPed = recipient.getConnectedNodes(new ArrayList<Node>());
		this.logLikelihood -= likelihoodLocalPedigree(recipientPed);
		
		
		//boolean onePed = recipientPed.contains(donor);
		if(!onePed){ //if two pedigrees, subtract donor pedigree as well
			clearVisit();
			List<Node> donorPed = donor.getConnectedNodes(new ArrayList<Node>());
			this.logLikelihood -= likelihoodLocalPedigree(donorPed);
		}

		
		//merge
		for(Node i : donor.getChildren()){
			i.removeParent(donor);
			i.addParent(recipient);
			recipient.addChild(i);
		}
		donor.getChildren().clear();
		
		for(Node i : donor.getParents()){
			i.removeChild(donor);
			i.addChild(recipient);
			recipient.addParent(i);
		}
		donor.getParents().clear();

		/*
		// update relationship matrix 
		for(Node i : inds){
			updateAdjMat(i);
		}
		*/
		
		
		//add new likelihood 
		if(onePed){ //if they were already in the same pedigree
			
			//update
			for(Node i : recipientPed){
				updateAdjMat(i);
			}
			
			this.logLikelihood += likelihoodLocalPedigree(recipientPed);
		}
		else{ //if two different pedigrees, get new merged pedigree
			
			clearVisit();
			List<Node> mergedPed = recipient.getConnectedNodes(new ArrayList<Node>());
			
			//update
			for(Node i : mergedPed){
				updateAdjMat(i);
			}
			
			this.logLikelihood += likelihoodLocalPedigree(mergedPed);
		}

		
		
	}

	
	//at least one of them has to be sampled
	
	public void swap(Node child, Node parent){//works
		
		//get cluster
		clearVisit();
		List<Node> ped = child.getConnectedNodes(new ArrayList<Node>());


		//subtract old terms
		this.logLikelihood -= likelihoodLocalPedigree(ped);		
		
		//switch nodes
		switchParentChild(parent, child);
	
		//update adj matrix 
		for(Node ind : ped){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood += likelihoodLocalPedigree(ped);
	
		
	}


	
	private void switchParentChild(Node parent, Node child){//works

		// remove child & parent from each other's list
		parent.removeChild(child);
		child.removeParent(parent);
		
		
		//shallow copy child's info
		List<Node> childParents = new ArrayList<Node>(child.getParents());
		List<Node> childChildren = new ArrayList<Node>(child.getChildren());
		
		
		//update child
		child.setDepth(parent.getDepth());
		child.setParents(parent.getParents());
		child.setChildren(parent.getChildren());
		
		for(Node p : child.getParents()){
			p.removeChild(parent);
			p.addChild(child);
		}
		for(Node c : child.getChildren()){
			c.removeParent(parent);
			c.addParent(child);
		}
		
		child.addChild(parent);
		
		
		//update parents
		parent.setDepth(parent.getDepth()-1);
		parent.setChildren(childChildren);
		parent.setParents(childParents);
		
		for(Node p : parent.getParents()){
			p.removeChild(child);
			p.addChild(parent);
		}
		for(Node c : parent.getChildren()){
			c.removeParent(child);
			c.addParent(parent);
		}
		
		parent.addParent(child);
		
		
	}
	
	
	public void stretch(Node parent){
		
		//get cluster
		clearVisit();
		List<Node> ped = parent.getConnectedNodes(new ArrayList<Node>());

		//subtract old terms
		this.logLikelihood -= likelihoodLocalPedigree(ped);	

		
		//stretch
		List<Node> ghostParents = new ArrayList<Node>(parent.getChildren().size());
	
		for(Node child : parent.getChildren()){
			Node ghostParent = makeNewNode(parent.getDepth(), parent.getSex());
			ghostParents.add(ghostParent);
			
			child.removeParent(parent);
			child.addParent(ghostParent);
			ghostParent.addChild(child);
			ghostParent.addParent(parent);
		}
		
		parent.setDepth(parent.getDepth()+1);
		parent.setChildren(ghostParents);
		
		
		//update adj matrix 
		for(Node ind : inds){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood += likelihoodLocalPedigree(ped);
	
		
	}
	
	
	public void compress(Node grandParent){
		
		//get cluster
		clearVisit();
		List<Node> ped = grandParent.getConnectedNodes(new ArrayList<Node>());

		//subtract old terms
		this.logLikelihood -= likelihoodLocalPedigree(ped);	
		
		
		//compress
		//get grand children
		List<Node> grandChildren = new ArrayList<Node>();
		List<Node> ghostParents = new ArrayList<Node>(grandParent.getChildren().size());
		for(Node ghostParent : grandParent.getChildren()){
			
			ghostParents.add(ghostParent);
			
			//update every grand child
			for(Node grandChild : ghostParent.getChildren()){
				
				grandChildren.add(grandChild);
			
			}
			
		}
		
		//delete ghost parents
		for(Node gc : ghostParents){
			deleteNode(gc);
		}

		
		//update edges
		for(Node gc : grandChildren){
			gc.addParent(grandParent);
			grandParent.addChild(gc);
		}
		grandParent.setDepth(grandParent.getDepth()-1);
		

		
		//update adj matrix 
		for(Node ind : inds){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood += likelihoodLocalPedigree(ped);
		
		
	}

	
	
	///////// UPDATE ADJACENCY MATRIX ////////
	public void updateAdjMat(Node node){
		
		if(!node.sampled) return;
		
		updatePathFromThisNode(node);
		updateAdjMatForThisNode(node, 0);
		
	}
	
	//update paths for this node
	public void updateAdjMatForThisNode(Node node, int offset){ //works
		
		if(!node.sampled) return;
		
		int i = node.getIndex();
		
		for(Node ind : inds){
		
			int j = ind.getIndex();
			
			if(i==j) //skip self
				continue; 
			
			if(ind.numVisit==0){ //unrelated
				relationships[i][j].updatePath(0,0,0);
				relationships[j][i].updatePath(0,0,0);
			}
			else{ //update relationship
				relationships[i][j].updatePath(ind.up + offset, ind.down, ind.numVisit);
				relationships[j][i].updatePath(ind.down, ind.up + offset, ind.numVisit);	
			}
			
		}
	
		
	}
 

	
	
	
	//////// PATH Update /////////
	//updates the relationship between the given node to everyone else
	public void updatePathFromThisNode(Node node){ //works

		  //record path from source node to every related node
		  clearVisit();
		  node.numVisit = 1;
		  updateUpDownPath(node, 1); //up down
		  updateDownPath(node, null, 0, 1); //down

	}


	//records the path from node to its relatives via its parents
	private void updateUpDownPath(Node node, int up){
					
		//for every parent
		for (Node parent :node.getParents()){
			
			//update path to parent
			if (parent.sampled){
				parent.recordPath(up, 0);
			}
				
			//update children of this parent
			updateDownPath(parent, node, up, 1);
				
			//recurse
			updateUpDownPath(parent, up+1);
			
		}

		
	}
	


	//records the path from node to its descendants, excluding excludeChild
	private void updateDownPath(Node node, Node excludeChild, int up, int down){

		//for every child
		for (Node child : node.getChildren()){
			
			if (child==excludeChild) continue;
			
			//update path to this child
			if (child.sampled){
				child.recordPath(up, down);
			}
			
			//recurse
			updateDownPath(child, null, up, down+1);
				
		}

	}

	
	
	///////// UPDATE LIKELIHOOD //////////
	//return the lkhd of the local pedigree
	public double likelihoodLocalPedigree(List<Node> connectedSamples){//works
			
		
		int n = connectedSamples.size();
		
		if(n==0) return 0d;
		
		double lkhd = 0d;
		
		//pairwise
		for(int i=0; i<connectedSamples.size(); i++){
			Node ind1 = connectedSamples.get(i);
			
			for(int j=i+1; j<connectedSamples.size(); j++){
				Node ind2 = connectedSamples.get(j);
				
				lkhd += core.getLikelihood(ind1, ind2, relationships[ind1.getIndex()][ind2.getIndex()]);
				
			}
		}
		
		//denom or founder
		double marginals = 0d;
		for(Node ind : connectedSamples){
			marginals += core.getMarginal(ind);
		}
		int coeff = n>1? -(n-2) : 1;
		lkhd += coeff * marginals;
		
		
		
		//TODO testing
		//lkhd += ageLikelihood(connectedSamples);
		

		return lkhd;
		
	}
	
	
	public double ageLikelihood(List<Node> connectedSamples){
		
		double toReturn = 0d;
		
		for(int i=0; i<connectedSamples.size(); i++){
			Node ind1 = connectedSamples.get(i);
			
			if(ind1.getAge()==-1) continue;
			
			for(int j=i+1; j<connectedSamples.size(); j++){
				Node ind2 = connectedSamples.get(j);
				
				if(ind2.getAge()==-1) continue;
				
				Path rel = relationships[ind1.getIndex()][ind2.getIndex()];
				if(rel.getNumVisit()==0) continue; //unrelated
				
				int up = rel.getUp();
				int down = rel.getDown();
				
				
				NormalDistribution normalDist = new NormalDistribution(ind1.getAge() + (up-down)*muGenTime, Math.sqrt((up+down)*varGenTime)); //TODO precompute
				
				toReturn += Math.log(normalDist.density(ind2.getAge()));
				
			}
		}
		
		
		return toReturn;
		
	}
	
	
	
	
	/////////////////// REVERSE MOVE /////////////////////
	public void reverseMove(){
		

	}

	
	public void copyCurrPedigree(){
		
	}
	
	
	
	//////////// MISC ////////////		
	public void clearVisit(){
		for (int i=0; i<this.nActiveNodes; i++){
			this.nodes[i].numVisit = 0;
		}
	}
	
	
	//marks reachable nodes as visited
	public void performDFS(Node node){
		
		clearVisit();
		dfs(node);
	}
	
	
	
	private void dfs(Node node){
		
		node.numVisit++;
		
		for(Node i : node.getParents()){
			if(i.numVisit > 0) continue;
			else dfs(i);
		}
		
		for(Node i : node.getChildren()){
			if(i.numVisit > 0) continue;
			else dfs(i);
		}
		
	}
	

	
	public boolean fullSibs(Node child1, Node child2){//works
		
		if(child1==child2) //skip self 
			return false;
		
		List<Node> parents1 = child1.getParents();
		List<Node> parents2 = child2.getParents();
		
		if(parents1.size()!=2 || parents2.size()!=2){
			return false;
		}
		

		if((parents1.get(0)==parents2.get(0) && parents1.get(1)==parents2.get(1)) || (parents1.get(0)==parents2.get(1) && parents1.get(1)==parents2.get(0))){
			return true;
		}
		
		return false;

		
	}
	
	
	
	
	public Node getHighestNode(Node node){//works
		
		//i.e. the nodes has at least 3 edges, including the child that is about to be cut
		if(node.getParents().size()==0 || node.sampled || node.getNumEdges() > 2){
			return node;
		}
		
		else{
			
			for(Node i : node.getParents()){
				return getHighestNode(i);
			}
			
		}
		
		throw new RuntimeException("highestnode not found!");
	}
	
	
	// returns the maximum age of the descendants of given node; excludes given node
	public Node getDescendantWithMaxAge(Node node){//TODO make this more efficient
		
		Node toReturn = null;
		double currAge = -1;
		
		for(Node i : node.getDescendants(new ArrayList<Node>())){
			
			if(i.sampled && i.getAge() > currAge){
				toReturn = i;
				currAge = i.getAge();
			}
			
		}
		
		
		return toReturn;
	
	}
	
	
	public Node getAncestorWithMinAge(Node node){//TODO make this more efficient
		
		Node toReturn = null;
		double currAge = -1;
		
		for(Node i : node.getAncestors(new ArrayList<Node>())){
			
			if(i.sampled && i.getAge() < currAge){
				toReturn = i;
				currAge = i.getAge();
			}
			
		}
		
		
		return toReturn;
	
	}
	
		

	//count the number of descendants at each level
	public int[] depthToCount(Node node, int[] depthToCount){ //works
		
		if(node==null) return depthToCount;
		
		//count
		depthToCount[node.getDepth()] += 1;
		
			
		//recurse
		for(Node i : node.getChildren()){
			depthToCount(i, depthToCount);
		}
			
		
		return depthToCount;
		
		
	}
	
	
 	public void printAdjMat(){
		
		for(int i=0; i<relationships.length; i++){
			for(int j=0; j<relationships.length; j++){
				
				if(!(i<j)){
					System.out.print(String.format("(%s, %s, %s) ","-","-","-"));	
				}
				else{
					Path path = relationships[i][j];
					System.out.print(String.format("(%d, %d, %d) ", path.getUp(), path.getDown(), path.getNumVisit()));
				}
			}
			System.out.println();
		}
		System.out.println();
	}




 
 	

}
	 
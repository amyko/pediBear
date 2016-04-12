package dataStructures;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.NormalDistribution;

import utility.ArrayUtility;
import likelihood.PairwiseLikelihoodCoreStream2;

//This is the data structure that stores information about all of the individuals in the sample
//Essentially the nodes themselves encode the pedigree graph, and this class contains
//a sort of adjacency matrix, where instead of binary adjacencies, only the individuals are represented
//and the i,j^th entry represents the relationship between the i^th and j^th individuals
public class Pedigree {

	//pedigree variables
	public final int maxDepth;
	public final int maxDepthForSamples;
	public final int numIndiv;
	private final PairwiseLikelihoodCoreStream2 core;
	public final Random rGen;

	
	//for reverse move
	private final Path[][][] relationships; 
	private final List<ArrayList<Node>> nodes = new ArrayList<ArrayList<Node>>(2);
	private int[] nActiveNodes = new int[2];
	private double[] logLikelihood = new double[2];
	private int curr;
	private int copy;

	
	//for age
	public final double genTime;
	final double muGenTime = 29;
	final double varGenTime = 36;

	

	////// CONSTRUCTOR ///////
	//TODO handle known relationships
	public Pedigree(int maxDepth, int maxDepthForSamples, Node[] inds, PairwiseLikelihoodCoreStream2 core, String marginalPath, String lkhdPath, Random rGen, int maxNumNodes, double genTime) throws IOException{
		
		this.numIndiv = inds.length;
		this.maxDepth = maxDepth;
		this.maxDepthForSamples = maxDepthForSamples;
		this.genTime = genTime;
		this.core = core;
		this.rGen = rGen;
		this.curr = 0;
		this.copy = 1;
		

		//initialize everyone to be unrelated
		this.relationships = new Path[2][numIndiv][numIndiv];
		for(int c=0; c<2; c++){
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					this.relationships[c][i][j] = new Path(0,0,0);
					//this.relationships[c][j][i] = new Path(0,0,0);
				}
			}
		}
	 
		//update nodes list
		//add sampled nodes
		for(int c=0; c<2; c++){
			
			//initailize list
			nodes.add(new ArrayList<Node>(maxNumNodes));
			
			//add sample nodes
			for(int i=0; i<inds.length; i++){
				Node ind = inds[i];
				nodes.get(c).add(new Node(true, i, ind.getSex(), ind.getDepth(), ind.getAge()));
			}
			
			//add ghost nodes
			for(int i = inds.length; i<maxNumNodes; i++){
				nodes.get(c).add(new Node(false, i));
			}
			
			//active nodes
			this.nActiveNodes[c] = inds.length;
		}
		

		
		
		//initialize pairwise & marginal likelihoods
		core.setMarginals(marginalPath);
		core.setLikelihoods(lkhdPath);
		
		//compute current likelihood
		for(Node i : inds){
			logLikelihood[curr] += core.getMarginal(i);
		}
		
		
		
	}
	
	
	////// SETTERS ////////
	public void setLogLikelihood(double lkhd){
		//this.logLikelihood = lkhd;
		logLikelihood[curr] = lkhd;
	}
	
	public void setNActiveNodes(int n){
		//this.nActiveNodes = n;
		this.nActiveNodes[curr] = n;
	}
	
	
	////// GETTERS ////////// 
	public double getLogLikelihood(){
		//return logLikelihood;
		return logLikelihood[curr];
	}
	
	
	public Path[][] getRelationships(){
		return relationships[curr];
	}
	
	public int getNActiveNodes(){
		return nActiveNodes[curr];
	}
	

	public Node getNode(int index){
		return nodes.get(curr).get(index);
	}

	
	public Node getRandomNode(){ //works
		return nodes.get(curr).get(rGen.nextInt(nActiveNodes[curr]));
	}
	
	
	public Node[] getNRandomNodes(int sampleSize){//works

		Node[] toReturn = new Node[sampleSize];

		int nSampled = 0;
		int i = 0;
		int N = this.nActiveNodes[curr];

		while(nSampled < sampleSize){
			
			if(rGen.nextDouble()*(N-i) < (sampleSize - nSampled)){
				toReturn[nSampled++] = this.nodes.get(curr).get(i);
			}
			
			i++;
			
		}
		
		return toReturn;
	}
	

	public Node getRandomSampledNode(){
		return nodes.get(curr).get(rGen.nextInt(numIndiv));
	}
	
	
	public Node getRandomPathAncestor(Node node, int targetDepth, int targetSex){
		
		if(node.getDepth() > targetDepth) 
			throw new RuntimeException("starting node has depth greater than target depth");

		Node currNode = node;
		
		//go up to target depth
		int sex;
		while(currNode.getDepth() < targetDepth){

			if(currNode.getDepth()==targetDepth-1){ //sex is determined by targetSex
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
	private void growNodeList(){
		
		int oldCapacity = nodes.get(curr).size();
		int newCapacity = (3 * oldCapacity) / 2 + 1;
		
		for(int i=oldCapacity; i<newCapacity; i++){
			nodes.get(curr).add(new Node(false, i));
		}
				
	}
	
	
	public Node makeNewNode(int depth, int sex){ //works

		//if out of new nodes, make more
		if(nActiveNodes[curr] == nodes.get(curr).size()){
			growNodeList();
		}

		//get new node from the end of array
		Node newNode = nodes.get(curr).get(nActiveNodes[curr]);
		newNode.setDepth(depth);
		newNode.setSex(sex);
		
		//increment pointer
		nActiveNodes[curr]++;
		
		return newNode;
	}
	


	public void deleteNode(Node nodeToDelete){ //works
	
		if(nodeToDelete.sampled)
			throw new RuntimeException("Trying to delete a sampled node");
		
		//clear edges and reset sex, depth
		nodeToDelete.reset();
		
		//update node list
		if(nodeToDelete.getIndex() != nActiveNodes[curr]-1){
			//Node lastNode = nodes.get(curr).get(this.nActiveNodes[curr]-1); //last active node
			Node lastNode = nodes.get(curr).set(this.nActiveNodes[curr]-1, nodeToDelete);
			nodes.get(curr).set(nodeToDelete.getIndex(), lastNode);
			
			//update indicies
			lastNode.setIndex(nodeToDelete.getIndex());
			nodeToDelete.setIndex(nActiveNodes[curr] - 1);
		}

		
		this.nActiveNodes[curr]--;
		
	}
	
	
	
	
	////////// UPDATE PEDIGREE STRUCTURE////////////	
 	public void connect(Node parent, Node child){
 		
 		parent.addChild(child);
 		child.addParent(parent);
 		
 	}
 	
 	
 	public void disconnect(Node parent, Node child){
 		
 		parent.removeChild(child);
 		child.removeParent(parent);
 		
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
		List<Node> nodesBeforeCut = parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeCut);
		
		
		// update graph structure
		child.removeParent(parent);
		parent.removeChild(child);
		
		
		for(Node ind : nodesBeforeCut){
			updateAdjMat(ind);
		}
		
		//update adjmat and likelihood
		if(hasFullSib){ //if child had a full sibling, the original pedigree is not split
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(nodesBeforeCut);
			
		}

		else{ //split pedigrees

			clearVisit();
			List<Node> childPed = child.getConnectedSampledNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> parentPed = parent.getConnectedSampledNodes(new ArrayList<Node>());

			this.logLikelihood[curr] += likelihoodLocalPedigree(childPed);
			this.logLikelihood[curr] += likelihoodLocalPedigree(parentPed);
		}
		
	}
	
	
	//make a ghost copy and randomly assign children to the copy
	public void split(Node parent, Node splitParent, List<Node> splitChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit = parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeSplit);
		
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
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(nodesBeforeSplit);
			
		}
		else{
			
			clearVisit();
			List<Node> parentPed = parent.getConnectedSampledNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> splitPed = splitParent.getConnectedSampledNodes(new ArrayList<Node>());
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(splitPed);
			this.logLikelihood[curr] += likelihoodLocalPedigree(parentPed);

		}
		       
		
		
	}
	

	
	//make a ghost copy and randomly assign children to the copy
	public void split2(Node parent, Node stayParent, List<Node> stayChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit =  parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeSplit);
		
		
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
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(nodesBeforeSplit);
			
		}
		else{
			clearVisit();
			List<Node> parentPed = parent.getConnectedSampledNodes(new ArrayList<Node>());
			clearVisit();
			List<Node> stayPed = stayParent.getConnectedSampledNodes(new ArrayList<Node>());
			
			//setUnrelated(parentPed, stayPed);
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(stayPed);
			this.logLikelihood[curr] += likelihoodLocalPedigree(parentPed);
			
		}
		
		
	}
	
	
	//recipient gets donor's children & parents
	//onePed==true if donor and recipient are already connected
	public void merge(Node donor, Node recipient, boolean onePed) {//works
						
		//subtract current subpedigree likelihoods
		clearVisit();
		List<Node> recipientPed = recipient.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(recipientPed);
		
		
		//boolean onePed = recipientPed.contains(donor);
		if(!onePed){ //if two pedigrees, subtract donor pedigree as well
			clearVisit();
			List<Node> donorPed = donor.getConnectedSampledNodes(new ArrayList<Node>());
			this.logLikelihood[curr] -= likelihoodLocalPedigree(donorPed);
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
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(recipientPed);
		}
		else{ //if two different pedigrees, get new merged pedigree
			
			clearVisit();
			List<Node> mergedPed = recipient.getConnectedSampledNodes(new ArrayList<Node>());
			
			//update
			for(Node i : mergedPed){
				updateAdjMat(i);
			}
			
			this.logLikelihood[curr] += likelihoodLocalPedigree(mergedPed);
		}

		
		
	}

	
	//at least one of them has to be sampled
	
	public void swap(Node child, Node parent){//works
		
		//get cluster
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());


		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);		
		
		//switch nodes
		switchParentChild(parent, child);
	
		//update adj matrix 
		for(Node ind : ped){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
	
		
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
	
	
	
	public void cutOneLinkTwo(Node child){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		
		//cut child from parent
		Node parent = child.getParents().get(0);
		disconnect(parent, child);
		
		//make a new ghost parent for child
		Node newParent = makeNewNode(parent.getDepth(), (parent.getSex()+1)%2); //choose the opposite sex; this way, there will no symmetry 
		connect(newParent, child);
		
		
		//connect parent and newParent to two grand parents
		List<Node> grandParents = parent.getParents();
		for(int i=0; i<2; i++){
			
			Node gp;
			
			if(i < grandParents.size()){ //grand parent exists
				gp = grandParents.get(i);
			}
			else{ // doesn't exist yet
				
				//make gp
				int gpSex = 0;
				
				if(i==1){
					gpSex = (grandParents.get(i-1).getSex()+1)%2;
				}
				
				gp = makeNewNode(parent.getDepth()+1, gpSex);
				
				//connect gp to parent
				connect(gp, parent);
			}
			
			//add to new parent
			connect(gp, newParent);
			
		}
		
		//sanity check
		if(grandParents.get(0).getSex()==grandParents.get(1).getSex()) throw new RuntimeException("Same sex parents");
		
		
		//adjust adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		
		
	}
	
	
	
	public void cutTwoLinkOne(Node parent, Node newParent){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = parent.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		
		
		//cut parent from grandparents & clean up
		for(Node gp : parent.getParents()){
			gp.removeChild(parent);
			clean(gp);
		}
		parent.getParents().clear();
		
		
		//cut child from parent & delete parent
		Node child = parent.getChildren().get(0);
		disconnect(parent, child);
		deleteNode(parent);
		
		//connect child with new parent
		connect(newParent, child);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		
		
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
		
		for(int j=0; j<numIndiv; j++){
		
			//int j = ind.getIndex();
			
			Node ind = nodes.get(curr).get(j);
			
			if(i==j) //skip self
				continue; 
			
			
			if(ind.getNumVisit()==0){ //unrelatedif
				if(i<j) relationships[curr][i][j].updatePath(0,0,0);
				else relationships[curr][j][i].updatePath(0,0,0);
			}
			else{ //update relationship
				if(i<j) relationships[curr][i][j].updatePath(ind.getUp() + offset, ind.getDown(), ind.getNumVisit());
				else relationships[curr][j][i].updatePath(ind.getDown(), ind.getUp() + offset, ind.getNumVisit());	
			}
			
		}
	
		
	}
 

	
	
	
	//////// PATH Update /////////
	//updates the relationship between the given node to everyone else
	public void updatePathFromThisNode(Node node){ //works

		  //record path from source node to every related node
		  clearVisit();
		  node.setNumVisit(1);
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
		int smaller;
		int bigger;
		
		//pairwise
		for(int i=0; i<connectedSamples.size(); i++){
			Node ind1 = connectedSamples.get(i);
			
			for(int j=i+1; j<connectedSamples.size(); j++){
				Node ind2 = connectedSamples.get(j);
				
				if(ind1.getIndex() < ind2.getIndex()){
					smaller = ind1.getIndex();
					bigger = ind2.getIndex();
				}
				else{
					smaller = ind2.getIndex();
					bigger = ind1.getIndex();
				}


				lkhd += core.getLikelihood(ind1, ind2, relationships[curr][smaller][bigger]);

				
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
				
				Path rel = relationships[curr][ind1.getIndex()][ind2.getIndex()];
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
	public void copyCurrPedigree(){
		
		//if not enough nodes in the copy nodeList, make more
		for(int i=nodes.get(copy).size(); i<nActiveNodes[curr]; i++){
			nodes.get(copy).add(new Node(false, i));
		}
		
		//reset extra nodes
		for(int i=nActiveNodes[curr]; i<nActiveNodes[copy]; i++){
			nodes.get(copy).get(i).reset();
		}
		
		
		//copy each active node
		for(int i=0; i<nActiveNodes[curr]; i++){
			
			Node modelNode = nodes.get(curr).get(i);
			Node copyNode = nodes.get(copy).get(i);
			
			copyNode.setDepth(modelNode.getDepth());
			copyNode.setSex(modelNode.getSex());

			
			//node pointers
			copyNode.getParents().clear();
			for(Node parent : modelNode.getParents()){
				copyNode.addParent(nodes.get(copy).get(parent.getIndex()));
			}
			
			copyNode.getChildren().clear();
			for(Node child : modelNode.getChildren()){
				copyNode.addChild(nodes.get(copy).get(child.getIndex()));
			}
			
		}

		//copy relationships
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				
				Path modelRel = relationships[curr][i][j];
				relationships[copy][i][j].updatePath(modelRel.getUp(), modelRel.getDown(), modelRel.getNumVisit());
			}
		}
		
		
		nActiveNodes[copy] = nActiveNodes[curr];
		logLikelihood[copy] = logLikelihood[curr];
		
		
	}
	
	
	public void reverse(){
		this.curr = copy;
		this.copy = (curr+1)%2;
	}
	
	
	//////////// MISC ////////////		
	public void clearVisit(){
		for (int i=0; i<this.nActiveNodes[curr]; i++){
			this.nodes.get(curr).get(i).setNumVisit(0);
		}
	}
	
	
	//marks reachable nodes as visited
	public void performDFS(Node node){
		
		clearVisit();
		dfs(node);
	}
	
	
	
	private void dfs(Node node){
		
		node.setNumVisit(node.getNumVisit()+1);
		
		for(Node i : node.getParents()){
			if(i.getNumVisit() > 0) continue;
			else dfs(i);
		}
		
		for(Node i : node.getChildren()){
			if(i.getNumVisit() > 0) continue;
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
	
	//returns number of full siblings child has with the given sex
	public List<Node> getFullSibs(Node child, int targetSex){
		
		List<Node> toReturn = new ArrayList<Node>();
		
		if(child.getParents().size() < 2){
			return toReturn;
		}
		
		clearVisit();
		
		for(Node parent : child.getParents()){
			
			for(Node sib : parent.getChildren()){
				
				if(sib.getNumVisit() > 0) continue;
				sib.setNumVisit(1);
				
				if(sib.getSex()==targetSex && sib.getChildren().size() > 0 && fullSibs(child, sib))
					toReturn.add(sib);
					
			}
			
		}
		
		return toReturn;
		
		
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
	
		

	public int[] getDepthToCount(Node node, int[] depthToCount){
		ArrayUtility.clear(depthToCount);
		return depthToCount(node, depthToCount);
	}
	
	
	//count the number of descendants at each level
	private int[] depthToCount(Node node, int[] depthToCount){ //works
		
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
		
		for(int i=0; i<numIndiv; i++){
			for(int j=0; j<numIndiv; j++){
				
				if(!(i<j)){
					System.out.print(String.format("(%s, %s, %s) ","-","-","-"));	
				}
				else{
					Path path = relationships[curr][i][j];
					System.out.print(String.format("(%d, %d, %d) ", path.getUp(), path.getDown(), path.getNumVisit()));
				}
			}
			System.out.println();
		}
		

	}

 	
 	
 	public boolean sanityCheck(){

 		for(int i=0; i<nActiveNodes[curr]; i++){
 			
 			Node node = nodes.get(curr).get(i);
 			
 			if(node.getParents().size()==2){
 				
 				Node p1 = node.getParents().get(0);
 				Node p2 = node.getParents().get(1);
 				
 				if(!(p1.getSex()==0 && p2.getSex()==1) && !(p1.getSex()==1 && p2.getSex()==0))
 				
 					return false;
 				
 			}
 			
 			
 		}
 		
 		
 		return true;
 		
 	}



 
 	

}
	 
package dataStructures;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

//import org.apache.commons.math3.distribution.NormalDistribution;




import utility.ArrayUtility;
import utility.DataParser;
import likelihood.PairwiseLikelihoodCoreStreamPed;

//This is the data structure that stores information about all of the individuals in the sample
//Essentially the nodes themselves encode the pedigree graph, and this class contains
//a sort of adjacency matrix, where instead of binary adjacencies, only the individuals are represented
//and the i,j^th entry represents the relationship between the i^th and j^th individuals
public class Pedigree {

	//pedigree variables
	public final int maxDepth;
	public final int maxSampleDepth;
	public final int numIndiv;
	private final PairwiseLikelihoodCoreStreamPed core;
	public final Random rGen;

	
	//for reverse move
	private final Path[][][] relationships; 
	private final List<ArrayList<Node>> nodes = new ArrayList<ArrayList<Node>>(2);
	private int[] nActiveNodes = new int[2];
	private double[] logLikelihood = new double[2];
	public int curr;
	private int copy;


	
	//for prior for unrelatedness
	private final double lambda;
	private final double logLambda;
	private final double[] logFact;
	public final int[] nSingletons;
	
	//for primus
	public boolean looped = false;

	
	
	////// CONSTRUCTOR ///////
	//for primus
	public Pedigree(String inPath, Map<String, Integer> name2Index) throws IOException{
		
		
		this.numIndiv = name2Index.size();
		this.maxDepth = 4;
		this.maxSampleDepth = 4;
		this.core = null;
		this.rGen = null;
		this.curr = 0;
		this.copy = 1;
		nActiveNodes[0] = 200;
		this.lambda = 0;
		this.logLambda = 0;
		logFact = null;
		nSingletons = null;
		
		//relationship
		this.relationships = new Path[2][numIndiv][numIndiv];

		for(int i=0; i<relationships[0][0].length; i++){
			for(int j=i+1; j<relationships[0][0].length; j++){
				this.relationships[0][i][j] = new Path(0,0,0);
			}
		}
		
		
		//set up pedigree

		//initialize list
		nodes.add(new ArrayList<Node>(200));
		
		//fill up nodes
		for(int i=0; i<200; i++){
			nodes.get(0).add(new Node("missing", String.format("%d", i), -1, true, i));
		}
		
		
		BufferedReader reader = DataParser.openReader(inPath);
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			int childIdx = name2Index.get(fields[1]);
			Node child = nodes.get(0).get(childIdx);
			
			if(!fields[2].equals("0")){
				int momIdx = name2Index.get(fields[2]);
				Node mom = nodes.get(0).get(momIdx);
				child.addParent(mom);
				mom.addChild(child);
			}
			if(!fields[3].equals("0")){
				int momIdx = name2Index.get(fields[3]);
				Node mom = nodes.get(0).get(momIdx);
				child.addParent(mom);
				mom.addChild(child);
			}
			
			
		}
		
		reader.close();
		
		
		//record paths
		for(String name : name2Index.keySet()){
			
			updateAdjMat(nodes.get(0).get(name2Index.get(name)));
			
		}

		
	}
	
	
	
	//this is for simulation only
	public Pedigree(String inPath, String outPath, int numIndiv, int[] ids) throws IOException{
		
		this.lambda = 0;
		this.logLambda = 0;
		this.nSingletons = null;
		logFact = null;

		//relationship
		this.relationships = new Path[2][184][184];

		for(int i=0; i<relationships[0][0].length; i++){
			for(int j=i+1; j<relationships[0][0].length; j++){
				this.relationships[0][i][j] = new Path(0,0,0);
			}
		}

	 
		
		this.numIndiv = numIndiv;
		this.maxDepth = 5;
		this.maxSampleDepth = this.maxDepth;
		this.core = null;
		this.rGen = null;
		this.curr = 0;
		this.copy = 1;
		nActiveNodes[0] = 200;
		
		//set up pedigree

		//initialize list
		nodes.add(new ArrayList<Node>(200));
		
		//fill up nodes
		for(int i=0; i<200; i++){

			nodes.get(0).add(new Node("1", String.format("%d", i+1), -1, false, i));
		}
		
		
		BufferedReader reader = DataParser.openReader(inPath);
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			int childIdx = Integer.parseInt(fields[0]);
			int momIdx = Integer.parseInt(fields[1]);
			int dadIdx = Integer.parseInt(fields[2]);
			
			Node child = nodes.get(0).get(childIdx);
			Node mom = nodes.get(0).get(momIdx);
			Node dad = nodes.get(0).get(dadIdx);
			
			mom.setSex(0);
			dad.setSex(1);
			
			
			child.addParent(dad);
			child.addParent(mom);
			mom.addChild(child);
			dad.addChild(child);
			
			
		}
		
		
		//record paths
		for(int i=0; i<ids.length; i++){
			
			updateAdjMat(nodes.get(0).get(ids[i]));
			
		}
		
		
		//write to path
		PrintWriter writer = DataParser.openWriter(outPath);
		
		this.clearVisit();
		for(int i=0; i<this.numIndiv; i++){
			recordFam(this.getNode(ids[i]), writer);
		}
		
		
		
		/*
		for(int i=0; i<ids.length; i++){
			for(int j=i+1; j<ids.length; j++){
				Path rel =  relationships[0][ids[i]][ids[j]];
				writer.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
			}
		}
		*/

		writer.close();
		
	}
	
	
	private void recordFam(Node ind, PrintWriter famWriter){
		
		//visit
		if(ind.getNumVisit()!=0) return;
		ind.setNumVisit(1);
		
		String name = ind.iid;
		String pa = "0";
		String ma = "0";
		
		//get parent ids
		for(Node parent : ind.getParents()){
			
			recordFam(parent, famWriter);
		
			if(parent.getSex()==0)
				ma = parent.iid;
			else if(parent.getSex()==1)
				pa = parent.iid;
			else
				throw new RuntimeException("Parent with unknown sex");
			
		}
		
		//write to file
		famWriter.write(String.format("%s\t%s\t%s\n", name, pa, ma));
		
		
	}
	
	
	
	//TODO handle known relationships
	public Pedigree(String fileName, PairwiseLikelihoodCoreStreamPed core, int maxDepth, int maxSampleDepth, Random rGen, int maxNumNodes, double lambda, int numIndiv, Map<String, Double> name2Age) throws IOException{
		
		this.numIndiv = numIndiv;
		this.maxDepth = maxDepth;
		this.maxSampleDepth = maxSampleDepth;
		this.lambda = lambda;
		this.logLambda = lambda;
		this.core = core;
		this.rGen = rGen;
		this.curr = 0;
		this.copy = 1;
		this.logFact = new double[numIndiv+1];
		this.nSingletons = new int[2];
		nSingletons[curr] = numIndiv;
		
		//log factorials
		logFact[0] = 0;
		for(int i=1; i<logFact.length; i++){
			logFact[i] = logFact[i-1] + Math.log(i);
		}
		

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
		//initialize list
		nodes.add(new ArrayList<Node>(maxNumNodes)); //current
		nodes.add(new ArrayList<Node>(maxNumNodes)); //copy
		
		//add sampled nodes by reading in .tfam
		BufferedReader famFile = DataParser.openReader(fileName+".tfam");
		
		String line;
		int index = 0;
		while((line=famFile.readLine())!=null){
			
			String[] fields = line.split("\\s");
			String fid = fields[0];
			String iid = fields[1];
			int sex = Integer.parseInt(fields[4]) - 1;
			
			//get age
			double age = -1;
			String name = fid+iid;
			if(name2Age!=null && name2Age.containsKey(name)){
				age = name2Age.get(name);
			}
			
			nodes.get(0).add(new Node(fid, iid, sex, true, age, 0, index));
			nodes.get(1).add(new Node(fid, iid, sex, true, age, 0, index));
			index++;
			
			//TODO fix
			if(index == numIndiv) break;
			
		}
		
		
		
		//add ghost nodes
		for(int c=0; c<2; c++){

			//add ghost nodes
			for(int i = numIndiv; i<maxNumNodes; i++){
				nodes.get(c).add(new Node("missing", String.format("%d", i), -1, false, i));
			}
			
			//active nodes
			this.nActiveNodes[c] = numIndiv;
		}
		

		
		
		//initialize pairwise & marginal likelihoods
		core.setMarginals(fileName+".marginal");
		core.setLikelihoods(fileName+".pairwise");
		
		//compute current likelihood
		for(int i=0; i<numIndiv; i++){
			logLikelihood[curr] += core.getMarginal(nodes.get(0).get(i));

		}
		
		logLikelihood[curr] += getSingletonProb();
		
		
		
		
	}
	
	
	public Pedigree(String fileName, PairwiseLikelihoodCoreStreamPed core, int maxDepth, int maxSampleDepth, Random rGen, int maxNumNodes, double lambda, int numIndiv) throws IOException{
	
		this(fileName, core, maxDepth, maxSampleDepth, rGen, maxNumNodes, lambda, numIndiv, null);
		
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
			nodes.get(curr).add(new Node("missing", String.format("%d", i), -1, false, i));
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

		this.logLikelihood[curr] -= getSingletonProb();
		
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
		
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
	}
	
	
	//make a ghost copy and randomly assign children to the copy
	public void split(Node parent, Node splitParent, List<Node> splitChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit = parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeSplit);
		this.logLikelihood[curr] -= getSingletonProb();
		
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
		
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		       
		
		
	}
	

	
	//make a ghost copy and randomly assign children to the copy
	public void split2(Node parent, Node stayParent, List<Node> stayChildren, boolean hasFullSib){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit =  parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeSplit);
		
		//subtract prior
		this.logLikelihood[curr] -= getSingletonProb();
		
		
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

			
			this.logLikelihood[curr] += likelihoodLocalPedigree(stayPed);
			this.logLikelihood[curr] += likelihoodLocalPedigree(parentPed);
			
		}
		
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
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
		
		//subtract prior
		this.logLikelihood[curr] -= getSingletonProb();

		
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
	
		
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();

		
		
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
		

		clean(parent);
		clean(child);
	
		
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
		this.logLikelihood[curr] -= getSingletonProb();
		
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
		//if(grandParents.get(0).getSex()==grandParents.get(1).getSex()) throw new RuntimeException("Same sex parents");
		
		
		//adjust adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
	}
	
	
	
	public void cutTwoLinkOne(Node parent, Node newParent){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = parent.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();
		
		
		//cut parent from grandparents & clean up
		for(Node gp : parent.getParents()){
			gp.removeChild(parent);
			//clean(gp);
			
			//clean gp
			if(!gp.sampled && gp.getNumEdges() < 2)
				deleteNode(gp);
			
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
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	public void shiftCluster(List<Node> cluster, int offset){
		
		//subtract old likelihood
		//this.logLikelihood[curr] -= ageLikelihood(cluster);
		
		//shift cluster
		for(Node i : cluster){
			i.setDepth(i.getDepth() + offset);
		}
		

		//add new likelihood
		//this.logLikelihood[curr] += ageLikelihood(cluster);
		
		
		
	}
	
	
	public void greatUncleToCousin(Node child, Node sibChild){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();
		


		//cut from parents
		Node p1 = child.getParents().get(0);
		Node p2 = child.getParents().get(1);
		disconnect(p1, child);
		disconnect(p2, child);
		
		
		//shift child cluster down
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() - 2);
		}
		
		//make a ghost node for child
		int parentSex = rGen.nextDouble() < .5 ? 0 : 1;
		Node newParent = makeNewNode(child.getDepth() + 1 , parentSex);
		connect(newParent, child);
		
		
		//connect new parent to grand parents
		List<Node> gp = sibChild.getParents();
		for(int i=0; i<2; i++){
			
			if(i < gp.size()){
				connect(gp.get(i), newParent);
			}
			else{
				Node newGP = makeNewNode(gp.get(0).getDepth(), (gp.get(0).getSex()+1)%2);
				connect(newGP, sibChild);
				connect(newGP, newParent);
			}
			
		}
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//clean
		clean(p1);
		clean(p2);
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();

		
		
	}
	
	
	
	public void cousinToGreatUncle(Node child, Node newSib){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		//disconnect parent
		Node parent = child.getParents().get(0);
		Node gp1 = parent.getParents().get(0);
		Node gp2 = parent.getParents().get(1);
		deleteNode(parent);
		
			
		//shift child cluster up
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() + 2);
		}
		
				
		//make sib and child siblings
		List<Node> newParents = newSib.getParents();
		for(int i=0; i<2; i++){
			
			if(i < newParents.size()){
				connect(newParents.get(i), child);
			}
			else{
				
				int parentSex = 0;
				if(i==1){
					parentSex = (newParents.get(0).getSex()+1)%2;
				}
				
				Node newParent = makeNewNode(child.getDepth()+1, parentSex);
				connect(newParent, newSib);
				connect(newParent, child);
			}
			
		}
		
		//clean
		clean(gp1);
		clean(gp2);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
		
	}
	
	
	
	public void POtoFS(Node child, Node parent, boolean goUp){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		//cut child from parent
		disconnect(parent, child);
		
		
		//shift cluster
		int delta = goUp ? 1 : -1;
		Node shiftNode = goUp ? child : parent;

		clearVisit();
		List<Node> shiftCluster = shiftNode.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + delta);
		}	

		
		
		//make sib and parent full siblings
		List<Node> gp = parent.getParents();
		for(int i=0; i<2; i++){
			
			//make new node
			if(i >= gp.size()){
				
				int targetSex = 0;
				if(i==1){
					targetSex = (gp.get(0).getSex()+1) % 2;
				}
				
				Node p1 = makeNewNode(child.getDepth()+1,targetSex);
				connect(p1, child);
				connect(p1, parent);
			}
			
			else{//connect to existing node
				Node p1 = gp.get(i);
				connect(p1, child);
			}
			
		}
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	
	public void FStoPO(Node child, Node parent, boolean goUp){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		//cut child from parent
		Node p1 = child.getParents().get(0);
		Node p2 = child.getParents().get(1);
		disconnect(p1, child);
		disconnect(p2, child);
		
		
		//clean grand parents, if necessary
		List<Node> grandParents = new ArrayList<Node>();
		grandParents.addAll(parent.getParents());
		for(Node gp : grandParents){
			if(!gp.sampled && gp.getNumEdges() < 2)
				deleteNode(gp);
		}
		
		
		
		//shift cluster
		int delta = goUp ? 1 : -1;
		Node shiftNode = goUp ? parent : child;

		clearVisit();
		List<Node> shiftCluster = shiftNode.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + delta);
		}			


		
		//connect child and parent
		connect(parent, child);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	public void halfUncleToCousin(Node child, Node halfSib){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();
		
		
		//cut from parent
		disconnect(child.getParents().get(0), child);
		
		
		//shift child cluster down
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() - 1);
		}
		
		//make a ghost parent node for child
		int parentSex = rGen.nextDouble() < .5 ? 0 : 1;
		Node newParent = makeNewNode(child.getDepth() + 1 , parentSex);
		connect(newParent, child);
		
		
		//connect new parent to grand parents
		List<Node> gp = halfSib.getParents();
		for(int i=0; i<2; i++){
			

			if(i < gp.size()){
				connect(gp.get(i), newParent);
			}
			else{
				Node newGP = makeNewNode(gp.get(0).getDepth(), (gp.get(0).getSex()+1)%2);
				connect(newGP, halfSib);
				connect(newGP, newParent);
			}
			
		}
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	
	public void cousinToHalfUncle(Node child){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		//choose sib
		Node parent = child.getParents().get(0);
		List<Node> gp = new ArrayList<Node>();
		gp.addAll(parent.getParents());


		//disconnect parent
		deleteNode(parent);
			
		//shift child cluster up
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() + 1);
		}
		
		//choose a new parent
		Node newParent = gp.get(rGen.nextInt(2));
				
		//connect new parent to child
		connect(newParent, child);
		

		//clean
		for(Node i : gp){
			clean(i);
		}
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	
	public void halfGreatUncleToHalfCousin(Node child, Node newGP){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();
		
		
		//cut from parent
		Node oldParent = child.getParents().get(0);
		disconnect(oldParent, child);
		
		
		//shift child cluster down
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() - 2);
		}
		
		//make a ghost parent node for child
		int parentSex = rGen.nextDouble() < .5 ? 0 : 1;
		Node newParent = makeNewNode(child.getDepth() + 1 , parentSex);
		connect(newParent, child);
		
		
		//connect new parent to grand parents
		connect(newGP, newParent);
		
		
		//clean old parent
		clean(oldParent);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	public void halfCousinToHalfGreatUncle(Node child, Node newSib){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();


		//disconnect parent
		Node parent = child.getParents().get(0);
		deleteNode(parent);
			
		//shift child cluster up
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() + 2);
		}
		
		
		//choose new parent's sex
		int newParentSex = rGen.nextDouble() < .5 ? 0 : 1;
		
		
		//get new parent with target sex
		Node newParent = newSib.getParentWithSex(newParentSex);
		if(newParent==null){
			newParent = makeNewNode(child.getDepth()+1, newParentSex);
			connect(newParent, newSib);
		}
		
		
		//connect new parent to child
		connect(newParent, child);
		

		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
		
	}
	
	
	public void contract(Node parent, Node child){

		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();


		//disconnect child cluster
		this.disconnect(parent, child);
			
		//shift child cluster up
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() + 1);
		}
		
		//new parents and children
		for(Node p : parent.getParents()){
			this.connect(p, child);
		}
		for(Node c : parent.getChildren()){
			this.connect(child, c);
		}
	
		
		//delete old parent to child
		deleteNode(parent);

		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
		
	}
	
	
	
	public void stretch(Node child){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		Node newParent = this.makeNewNode(child.getDepth(), child.getSex());
		

		//connect new parent, disconnect child cluster
		for(Node p : child.getParents()){
			this.connect(p, newParent);
			p.removeChild(child);
		}
		child.getParents().clear();
		
			
		//shift child cluster down
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() - 1);
		}

		//connect new parent to everyone
		this.connect(newParent, child);
		

		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
		
		
	}
	
	
	
	//both are sampled; not parent-offspring
	public void swapAncDesc(Node anc, Node desc){
		
		//get cluster
		clearVisit();
		List<Node> ped = desc.getConnectedSampledNodes(new ArrayList<Node>());


		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);	

				
		//switch nodes
		//save descendant parents & children
		List<Node> descParents = new ArrayList<Node>();
		List<Node> descChildren = new ArrayList<Node>();
		for(Node p : desc.getParents()){
			descParents.add(p);
			p.removeChild(desc);
			p.addChild(anc);
		}
		for(Node p : desc.getChildren()){
			descChildren.add(p);
			p.removeParent(desc);
			p.addParent(anc);
		}
		int descDepth = desc.getDepth();
		
		desc.getParents().clear();
		desc.getChildren().clear();
		
		//save anc parents & children
		List<Node> ancParents = new ArrayList<Node>();
		List<Node> ancChildren = new ArrayList<Node>();
		for(Node p : anc.getParents()){
			ancParents.add(p);
			p.removeChild(anc);
			p.addChild(desc);
		}
		for(Node p : anc.getChildren()){
			ancChildren.add(p);
			p.removeParent(anc);
			p.addParent(desc);
		}
		
		anc.getParents().clear();
		anc.getChildren().clear();
		
		//update
		anc.setParents(descParents);
		anc.setChildren(descChildren);
		desc.setParents(ancParents);
		desc.setChildren(ancChildren);
		desc.setDepth(anc.getDepth());
		anc.setDepth(descDepth);
		
		
		
	
		//update adj matrix 
		for(Node ind : ped){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
	
	
		
	}
	
	
	public void halfSibstoFullUncle(Node child, Node halfSib){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();
		
		
		//cut from parent
		disconnect(child.getParents().get(0), child);
		
		
		//shift child cluster down
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() - 1);
		}
		
		//make a ghost parent node for child
		int parentSex = rGen.nextDouble() < .5 ? 0 : 1;
		Node newParent = makeNewNode(child.getDepth() + 1 , parentSex);
		connect(newParent, child);
		
		
		//connect new parent to existing gp
		for(Node gp: halfSib.getParents()){
			connect(gp, newParent);	
		}
		
		//make new gp if needed
		if(halfSib.getParents().size()==1){
			
			//connect new parent to new grand parent
			Node newGp = makeNewNode(halfSib.getDepth()+1, (halfSib.getParents().get(0).getSex()+1)%2);
			connect(newGp, newParent);
			
			//connect halfsib to newGP
			connect(newGp, halfSib);
		}
		

		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
	}
	
	
	
	public void fullUncletoHalfSibs(Node child, Node uncle){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		this.logLikelihood[curr] -= getSingletonProb();

		
		//disconnect parent
		Node parent = child.getParents().get(0);
		deleteNode(parent);
			
		//shift child cluster up
		clearVisit();
		List<Node> childCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : childCluster){
			i.setDepth(i.getDepth() + 1);
		}
		
		//choose gp to disconnect from
		int badSex = rGen.nextDouble() < .5 ? 0 : 1;
		int goodSex = (badSex+1) % 2;
		Node badGP = uncle.getParentWithSex(badSex);
		Node goodGP = uncle.getParentWithSex(goodSex);


		
		//disconnect uncle from badGP if it doesn't have any full siblings
		if(getFullSibs(uncle).size()==0){
			disconnect(badGP, uncle);	
		}
		
		//connect to new parent
		connect(goodGP, child);
		

		//clean
		clean(badGP);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updateNumSingletons();
		this.logLikelihood[curr] += getSingletonProb();
		
		
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
		
		for(int j=0; j<relationships[0][0].length; j++){
		
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

		//System.out.print(String.format("%d:\t", node.getIndex()));
		
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
				boolean isLooped = parent.recordPath(up, 0);

				if(isLooped) this.looped = true;
				
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
				boolean isLooped = child.recordPath(up, down);	
				
				if(isLooped) this.looped = true;
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
		
		if(n==1)
			return core.getMarginal(connectedSamples.get(0));
		
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
		


		return lkhd;
		
	}
	
	
	private double getSingletonProb(){
		
		return this.nSingletons[curr]*logLambda - lambda - logFact[this.nSingletons[curr]];

		
		//return 0;
		
		//return this.nSingletons[curr]*Math.log(numIndiv)/numIndiv;
		
	}
	
	

	public double likelihoodAllPedigrees(){
		
		double toReturn = 0d;
		clearVisit();
		
		
		for(int i=0; i<numIndiv; i++){
			
			Node node = nodes.get(curr).get(i);
			
			if(node.getNumVisit() > 0) continue;
	
			List<Node> connected = new ArrayList<Node>();
			//node.getConnectedSampledNodes(connected);
			getConnectedSampledNodes(node, connected);
			
			toReturn += likelihoodLocalPedigree(connected);
			
			/*
			for(Node j  : connected){
				System.out.print(String.format("%d\t", j.getIndex()));
			}
			System.out.println();
			*/
			
			
		}
		
		//System.out.println(String.format("Number of clusters: %d", n));
		
		return toReturn + getSingletonProb();
		
		
	}
	
	//get connectd nodes from relationship matrix
	private List<Node> getConnectedSampledNodes(Node node, List<Node> toReturn){
		
		node.setNumVisit(1);
		if(node.sampled) toReturn.add(node);
		
		//recurse on neighbors
		for(int i=0; i<numIndiv; i++){
			
			if(i==node.getIndex()) continue;
			
			Node neighbor = nodes.get(curr).get(i);
			if(neighbor.getNumVisit() > 0) continue;
			
			int bigger = node.getIndex() > i ? node.getIndex() : i;
			int smaller = node.getIndex() > i ? i : node.getIndex();
			
			
			if(relationships[curr][smaller][bigger].getNumVisit()==0){
				continue;
			}
			
			getConnectedSampledNodes(neighbor, toReturn);
			
		}
		
		return toReturn;
		
	}
	
	
	/*
	private double ageLikelihood(List<Node> connectedSamples){
		
		double toReturn = 0d;
		
		for(Node i : connectedSamples){

				if(i.getAge()==-1) continue;
			
				NormalDistribution normalDist = new NormalDistribution((i.getDepth()+1)*muGenTime, varGenTime); //TODO precompute
				
				toReturn += Math.log(normalDist.density(i.getAge()));
				
			
		}
		
		
		return toReturn;
		
	}
	*/
	
	
	public double totalLikelihood(){
		
		double toReturn = 0d;
		
		clearVisit();
		
		for(int i=0; i<numIndiv; i++){
			
			Node node = nodes.get(curr).get(i);
			
			if(node.getNumVisit() > 0) continue;
			
			List<Node> ped = node.getConnectedSampledNodes(new ArrayList<Node>());
			
			toReturn += likelihoodLocalPedigree(ped);
			
			
		}
		
		return toReturn + getSingletonProb();
		
		
	}
	
	
	public double pairwiseLkhd(){
		
		double toReturn = 0d;
		
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				
				toReturn += core.getLikelihood(nodes.get(curr).get(i), nodes.get(curr).get(j), relationships[curr][i][j]);
				
			}
		}
		
		
		return toReturn/5.0;
		
	}
	
	/////////////////// REVERSE MOVE /////////////////////
	public void copyCurrPedigree(){
		
		//if not enough nodes in the copy nodeList, make more
		for(int i=nodes.get(copy).size(); i<nActiveNodes[curr]; i++){
			nodes.get(copy).add(new Node("missing", String.format("%d", i), -1, false, i));
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
		nSingletons[copy] = nSingletons[curr];
		
		
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
	
	
	public void updateNumSingletons(){
		
		/*
		//get number of singletons
		int k = 0;
		for(int i=0; i<numIndiv; i++){
			
			boolean singleton = true;
			
			//connection to i
			for(int j=0; j<numIndiv; j++){
				
				if(j==i) continue;
				int smaller;
				int bigger;
				if(j<i){
					smaller = j;
					bigger = i;
				}
				else{
					smaller = i;
					bigger = j;
				}
				
				if(this.relationships[curr][smaller][bigger].getNumVisit()!=0){
					singleton = false;
					break;
				}
				
			}
			
			if(singleton==true) k++;
			

			
		}
		
		this.nSingletons[curr] = k;
		*/
		
		//get number of clusters
		
		this.clearVisit();
		
		int numCluster = 0;
		
		for(int i=0; i<numIndiv; i++){
			
			Node ind = nodes.get(curr).get(i);
			
			if(ind.getNumVisit()>0) continue;
			
			numCluster++;
			ind.getConnectedNodes(new ArrayList<Node>());
			
		}
		
		nSingletons[curr] = numCluster;
		
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
	
	//returns number of full siblings child has with the given sex, excluding itself
	public List<Node> getFullSibsWithTargetSex(Node child, int targetSex){
		
		List<Node> toReturn = new ArrayList<Node>();
		
		if(child.getParents().size() < 2){
			return toReturn;
		}
		
		clearVisit();
		
		Node parent = child.getParents().get(0);
			
		for(Node sib : parent.getChildren()){
			
			if(sib.getNumVisit() > 0) continue;
			sib.setNumVisit(1);
			
			if(sib.getSex()==targetSex && sib.getChildren().size() > 0 && fullSibs(child, sib))
				toReturn.add(sib);
				
		}
			
		
		
		return toReturn;
		
		
	}
	
	
	//returns number of full siblings child has with the given sex, excluding itself
	public List<Node> getFullSibs(Node child){
		
		List<Node> toReturn = new ArrayList<Node>();
		
		if(child.getParents().size() < 2){
			return toReturn;
		}

		Node parent = child.getParents().get(0);
			
		for(Node sib : parent.getChildren()){
			
			if(fullSibs(child, sib))
				toReturn.add(sib);
				
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
 			
 			
 			//parent sex
 			if(node.getParents().size()==2){
 				
 				Node p1 = node.getParents().get(0);
 				Node p2 = node.getParents().get(1);
 				
 				if(!(p1.getSex()==0 && p2.getSex()==1) && !(p1.getSex()==1 && p2.getSex()==0)){
 					System.out.println("Parent error!");
 					return false;
 				}
 				
 			}
 			
 			
 			//parent number
 			if(node.getParents().size() > 2){
 				System.out.println("too many parents");
 				return false;
 			}
 				
 			
 			//depth consistency
 			for(Node k : node.getParents()){
 				if(k.getDepth() != node.getDepth()+1){
 					System.out.println("depth error!");
 					return false;
 				}
 			}
 			
 			for(Node k : node.getChildren()){
 				if(k.getDepth() != node.getDepth()-1){
 					return false;
 				}
 			}
 			
 			
 			if(node.getDepth() > maxDepth || node.getDepth() < 0){
 				System.out.println("depth error!");
 				return false;
 			}
 			
 			
 			//ghost nodes
 			if(!node.sampled && node.getNumEdges()<2){
 				System.out.println("ghost error!");
 				return false;
 			}
 			
 			//age consistency
 			if(node.getAge() != -1){
 				
 				for(Node p : node.getParents()){
 					
 					if(p.getAge()!=-1 && p.getAge() < node.getAge()){
 						System.out.println("Age error!");
 						return false;
 					}
 
 				}
 				
 			}
 			
 			
 		}
 		
 		
 		//likelihood consistency
 		if(Math.abs(totalLikelihood() - logLikelihood[curr]) > 1e-1){
 			System.out.println("lkhd error!");
 			System.out.println(Math.abs(totalLikelihood() - logLikelihood[curr]));
 			return false;
 		}

 		
 		
 		return true;
 		
 	}
 	
 	
	public int getMaxDepth(Node node){

		int maxDepth = node.getDepth();
		node.setNumVisit(1);
		
		for(Node c : node.getChildren()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMaxDepth(c);
			if(currDepth > maxDepth)
				maxDepth = currDepth;
		}
		for(Node c : node.getParents()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMaxDepth(c);
			if(currDepth > maxDepth)
				maxDepth = currDepth;
		}
		
		
		return maxDepth;
		
		
	}
	
	
	public int getMinDepth(Node node){

		int minDepth = node.getDepth();
		node.setNumVisit(1);
		
		for(Node c : node.getChildren()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMinDepth(c);
			if(currDepth < minDepth)
				minDepth = currDepth;
		}
		for(Node c : node.getParents()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMinDepth(c);
			if(currDepth < minDepth)
				minDepth = currDepth;
		}

		
		return minDepth;
		
		
	}
	
	public boolean isInbred(Node node){
		
		clearVisit();
		visitAncestors(node);
		List<Node> ancs = node.getAncestors(new ArrayList<Node>());
		
		for(Node i : ancs){
			if(i.getNumVisit() > 1)
				return true;
		}
		
		return false;
		
		
	}
	
	
	public void visitAncestors(Node node){
		
		node.setNumVisit(node.getNumVisit()+1);
		
		for(Node p : node.getParents()){
			visitAncestors(p);
		}
		
		
	}
	
	
	//returns true if the parent's sex cannot be changed
	public boolean sexLocked(Node parent){

		parent.setNumVisit(parent.getNumVisit()+1);
		
		if(parent.sampled) 
			return true;
		
		//recurse on neighbor parents
		else{
		
			for(Node c : parent.getChildren()){
				
				if(c.getNumVisit() > 0) continue;
				c.setNumVisit(c.getNumVisit()+1);
				
				
				for(Node p : c.getParents()){
					
					if(p.getNumVisit() > 0) continue;
					
					if(sexLocked(p))
						return true;
					
				}
				
								
			}
			
			return false;
		}

		
	}
 
 	

}
	 
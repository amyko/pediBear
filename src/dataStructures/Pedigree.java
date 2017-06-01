package dataStructures;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;

//import org.apache.commons.math3.distribution.NormalDistribution;








import java.util.Set;

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
	private double[] prior = new double[2];
	public int curr;
	private int copy;


	
	//for prior for unrelatedness
	private final double lambda;
	private final double logLambda;
	private final double[] logFact;
	public final int[] nSingletons;
	private final double beta; //multiplier for poisson 
	
	
	//NEW PRIOR
	private int effectivePop;
	private int minN;
	private int maxN;
	private int[] nNodes;
	private int[] nUnsampledDads;
	private int[] nUnsampledMoms;
	private int[] nSampledMoms;
	private int[] nSampledDads;
	private final double[] log;
	private int totalUnits;
	
	
	//for primus
	public boolean looped = false;

	
	
	////// CONSTRUCTOR ///////

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
	
	
	
	//TODO handle known relationships; effective pop size
	public Pedigree(String fileName, PairwiseLikelihoodCoreStreamPed core, int maxDepth, int maxSampleDepth, Random rGen, int maxNumNodes, double lambda, int numIndiv, Map<String, Double> name2Age, double beta, int minN, int maxN) throws IOException{
		
		this.numIndiv = numIndiv;
		this.maxDepth = maxDepth;
		this.maxSampleDepth = maxSampleDepth;
		this.lambda = lambda;
		this.logLambda = Math.log(lambda);
		this.core = core;
		this.rGen = rGen;
		this.curr = 0;
		this.copy = 1;
		this.nSingletons = new int[2];
		this.logFact = new double[numIndiv+1];
		nSingletons[curr] = numIndiv;
		nNodes = new int[maxDepth+1];
		nUnsampledDads = new int[maxDepth+1];
		nUnsampledMoms = new int[maxDepth+1];
		nSampledMoms = new int[maxDepth+1];
		nSampledDads = new int[maxDepth+1];
		log = new double[5000]; //TODO set this better
		this.beta = beta;
		this.minN = minN;
		this.maxN = maxN;
		this.effectivePop = (minN+maxN) / 2;
		
		for(int i=0; i<=maxDepth; i++) totalUnits+=(i+1);
		
		

		//log factorials
		logFact[0] = 0;
		for(int i=1; i<logFact.length; i++){
			logFact[i] = logFact[i-1] + Math.log(i);
		}
		
		//log
		log[0] = 0;
		for(int i=1; i<log.length; i++){
			log[i] = Math.log(i);
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
			String name = fid+"_"+iid;
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
		
		
		//prior
		//logLikelihood[curr] += getSingletonProb();
		updatePrior();
		logLikelihood[curr] += prior[curr];
		
		
	}
	
	
	public Pedigree(String fileName, PairwiseLikelihoodCoreStreamPed core, int maxDepth, int maxSampleDepth, Random rGen, int maxNumNodes, double lambda, int numIndiv, double beta, int minN, int maxN) throws IOException{
	
		this(fileName, core, maxDepth, maxSampleDepth, rGen, maxNumNodes, lambda, numIndiv, null, beta, minN, maxN);
		
	}

	
	
	//this is for simulation only
	public Pedigree(String inPath, String truePath) throws IOException{
		
		this.lambda = 0;
		this.logLambda = 0;
		this.nSingletons = null;
		logFact = null;
		log = new double[500]; //TODO set this better
		beta = 0;
		effectivePop = 1000;

	 
		
		this.numIndiv = 20;
		this.maxDepth = 5;
		this.maxSampleDepth = this.maxDepth;
		this.core = null;
		this.rGen = null;
		this.curr = 0;
		this.copy = 1;
		nodes.add(new ArrayList<Node>(300));
		

		//fill up nodess
		BufferedReader reader = DataParser.openReader(inPath);
		reader.readLine();//header
		String line;
		int idx = 0;
		Map<String, Node> name2node = new HashMap<String, Node>();
		Map<Integer, Node> idx2node = new HashMap<Integer, Node>();
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			String childIdx = fields[0];
			String fid = fields[0].split("_")[0];
			String iid = fields[0].split("_")[1];
			int sex = fields[4].equals("1") ? 0 : 1;
			boolean sampled = fields[4].equals("000000") ? true : false;
			
			
			Node child = new Node(fid, iid, sex, sampled, idx);
			nodes.get(0).add(child);
			name2node.put(childIdx, child);
			idx2node.put(idx, child);
			
			//increment
			idx++;

		}
		reader.close();
		
		//relationship
		this.relationships = new Path[2][idx][idx];
		nActiveNodes[0] = idx;

		for(int i=0; i<relationships[0][0].length; i++){
			for(int j=i+1; j<relationships[0][0].length; j++){
				this.relationships[0][i][j] = new Path(0,0,0);
			}
		}


		
		//make pedigree
		reader = DataParser.openReader(inPath);
		reader.readLine();//header
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			Node child = name2node.get(fields[0]);
			Node dad = name2node.get(fields[1]);
			Node mom = name2node.get(fields[2]);
			
			
			if(dad!=null){
				child.addParent(dad);
				dad.addChild(child);
			}
			if(mom!=null){
				child.addParent(mom);
				mom.addChild(child);
			}
			

		}
		reader.close();
			
		
		//un-inbred: disconnect inbred connections
		for(Node x : nodes.get(0)){
			
			if(x.getParents().size()!=2) continue;
			
			
			Node p1 = x.getParents().get(0);
			Node p2 = x.getParents().get(1);
			
			List<Node> gp1 = new ArrayList<Node>();
			gp1.addAll(p1.getParents());
			List<Node> gp2 = new ArrayList<Node>();
			gp2.addAll(p2.getParents());
			
			for(Node g1 : gp1){
				
				for(Node g2 : gp2){
					
					//disconnect inbred
					if(g1.getIndex()==g2.getIndex()){
						
						disconnect(g1, p1);
						disconnect(g1, p2);
						
					}
					
				}
				
			}
			
			
		}
		
	
		
		//record paths
		for(Node x : nodes.get(0)){
			
			if(!x.sampled) continue;
			
			updateAdjMat(x);
			
		}
		
		
		//write to path
		PrintWriter writer = DataParser.openWriter(truePath);

		
		//print relationships
		for(int i=0; i<relationships[0].length; i++){

			Node node1 = idx2node.get(i);
			if(node1.sampled==false) continue;
			
			
			for(int j=i+1; j<relationships[0].length; j++){
				
				Node node2 = idx2node.get(j);
				if(node2.sampled==false) continue;
				
				Path rel = relationships[0][i][j];
				if(rel.getNumVisit()==-1){
					//System.out.println("Bad");
					looped = true;
				}
				
				if(rel.getNumVisit()!=0){
					//System.out.println(String.format("%d %d %d", rel.getUp(),rel.getDown(),rel.getNumVisit()));
				}
				
				//write name1 name2 relationship
				String name1 = String.format("%s_%s", node1.fid, node1.iid);
				String name2 = String.format("%s_%s", node2.fid, node2.iid);
				writer.write(String.format("%s\t%s\t%d\t%d\t%d\n", name1, name2, rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
			}
			
		}

		writer.close();
		

		
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
	
	public double getPrior(){
		return prior[curr];
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
	
	
	public void makeFullSibs(Node nodeWithParents, Node loneNode){
		
		if(loneNode.getParents().size()!=0) throw new RuntimeException("Error in make fullSibs");
		
		//make ghost parents and nephew fullSibs
		int nGP = nodeWithParents.getParents().size();
		if(nGP==2){
			for(Node x : nodeWithParents.getParents()){
				connect(x, loneNode);
				
			}
		}
		else if(nGP==1){
			
			connect(nodeWithParents.getParents().get(0), loneNode);
			
			Node gp = makeNewNode(nodeWithParents.getDepth()+1, (nodeWithParents.getParents().get(0).getSex()+1)%2);
			connect(gp, nodeWithParents);
			connect(gp, loneNode);
		}
		else{
			for(int sex=0; sex<2; sex++){
				Node gp = makeNewNode(nodeWithParents.getDepth()+1, sex);
				connect(gp, nodeWithParents);
				connect(gp, loneNode);
			}
		}
		
		
	}
	
	
	
	//returns highest node touched
	public void cut(Node child, Node parent, boolean hasFullSib, Node[] ijPrime){ //works
		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeCut = parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeCut);

		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
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
		
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		
		ijPrime[0] = clean(child);
		ijPrime[1] = clean(parent);
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	
		
	}
	
	
	//make a ghost copy and randomly assign children to the copy
	public void split(Node parent, Node splitParent, List<Node> splitChildren, boolean hasFullSib, Node[] ijPrime){ //works

		
		//subtract old likelihood
		clearVisit();
		List<Node> nodesBeforeSplit = parent.getConnectedSampledNodes(new ArrayList<Node>());
		this.logLikelihood[curr] -= likelihoodLocalPedigree(nodesBeforeSplit);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
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
		
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		ijPrime[0] = clean(parent);
	    ijPrime[1] = clean(splitParent);
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		       
		
		
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
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];

		
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
	
		
		deleteNode(donor);
		
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
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
	
	
	
	public void op2po(Node lowerNode, Node middleNode, Node upperNode){
		
		//get cluster
		clearVisit();
		List<Node> ped = middleNode.getConnectedSampledNodes(new ArrayList<Node>());


		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);	
		this.logLikelihood[curr] -= prior[curr];
		
		
		//switch nodes
		if(upperNode==null){
			upperNode = makeNewNode(lowerNode.getDepth()+2, lowerNode.getSex());
			connect(upperNode, middleNode);
		}
		if(lowerNode==null){
			lowerNode = makeNewNode(upperNode.getDepth()-2, upperNode.getSex());
			connect(middleNode, lowerNode);
		}
		
		
		//save gp's edges
		List<Node> gpParents = new ArrayList<Node>(upperNode.getParents());
		List<Node> gpChildren = new ArrayList<Node>(upperNode.getChildren());
		
		//update incoming edges
		for(Node p : upperNode.getParents()){
			p.removeChild(upperNode);
			p.addChild(lowerNode);
		}
		
		for(Node c : upperNode.getChildren()){
			c.removeParent(upperNode);
			c.addParent(lowerNode);
		}
		for(Node p : lowerNode.getParents()){
			p.removeChild(lowerNode);
			p.addChild(upperNode);
		}
		
		for(Node c : lowerNode.getChildren()){
			c.removeParent(lowerNode);
			c.addParent(upperNode);
		}
		
		//update outgoing edges
		upperNode.setDepth(lowerNode.getDepth());
		upperNode.setParents(lowerNode.getParents());
		upperNode.setChildren(lowerNode.getChildren());
		lowerNode.setDepth(lowerNode.getDepth() + 2);
		lowerNode.setParents(gpParents);
		lowerNode.setChildren(gpChildren);
			
			
			
		clean(upperNode);
		clean(lowerNode);
		
	
		//update adj matrix 
		for(Node ind : ped){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
	
		
		
	}
	
	

	
	
	public void shiftCluster(List<Node> cluster, int offset){
		

		this.logLikelihood[curr] -= prior[curr];
		

		//shift cluster
		for(Node i : cluster){
			i.setDepth(i.getDepth() + offset);
		}
		
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
	
	}
	

	public void switchSex(Node parent){
	
		clearVisit();
		
		this.logLikelihood[curr] -= prior[curr];
	
		//shift cluster
		switchSexHelper(parent);
		
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}
	
	
	
	private void switchSexHelper(Node parent){ //works
		
		if(parent.getNumVisit() > 0) return;
		
		//switch sex and mark visit
		parent.setSex((parent.getSex()+1) % 2);
		parent.setNumVisit(parent.getNumVisit()+1);
		
		//recurse on neighbor parents
		for(Node i : parent.getChildren()){
			
			if(i.getNumVisit() > 0) continue;
			i.setNumVisit(i.getNumVisit()+1);

			for(Node j : i.getParents()){
				switchSexHelper(j);
			}
							
		}
	}
	
	
	public void POtoFS(Node child, Node parent, int shift){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		//cut child from parent
		disconnect(parent, child);
		
		
		//shift cluster
		if(shift==1){ //case1: shift child cluster up
			clearVisit();
			List<Node> cluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() + 1);
			}			
		}
		else{ //case2: shift parent cluster down
			clearVisit();
			List<Node> cluster = parent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}	
		}

		
		//make child and parent full siblings
		makeFullSibs(parent, child);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
	}
	
	
	
	public void FStoPO(Node child, Node parent, int shift){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		//cut child from parent
		Node p1 = child.getParents().get(0);
		Node p2 = child.getParents().get(1);
		disconnect(p1, child);
		disconnect(p2, child);
		
		/*
		//clean grand parents, if necessary
		List<Node> grandParents = new ArrayList<Node>();
		grandParents.addAll(parent.getParents());
		for(Node gp : grandParents){
			if(!gp.sampled && gp.getNumEdges() < 2)
				deleteNode(gp);
		}
		*/
		
		
		//shift cluster
		if(shift==-1){
			clearVisit();
			List<Node> cluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}			
		}
		else{
			clearVisit();
			List<Node> cluster = parent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() + 1);
			}	
		}

		
		//connect child and parent
		connect(parent, child);
		
		
		//clean grand parents, if necessary
		List<Node> grandParents = new ArrayList<Node>();
		grandParents.addAll(parent.getParents());
		for(Node gp : grandParents){
			if(!gp.sampled && gp.getNumEdges() < 2)
				deleteNode(gp);
		}
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}
	
	
	public void FStoSelf(Node donor, Node recipient){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = donor.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
		//donate children
		for(Node x : donor.getChildren()){
			x.removeParent(donor);
			x.addParent(recipient);
			recipient.addChild(x);
		}
		donor.getChildren().clear();
		
		clean(donor);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}
	
	
	
	public void selfToFS(Node parent, List<Node> splitChildren){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = parent.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
		//make split parent
		Node splitParent = makeNewNode(parent.getDepth(), parent.getSex());
		
		
		//donate children
		for(Node x : splitChildren){
			disconnect(parent, x);
			connect(splitParent, x);
		}

		
		//make parent & splitParents fullSibs
		makeFullSibs(parent, splitParent);
		
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
		
	}
	
	
	public void HStoPO(Node child, Node parent, Node hs, List<Node> fullSibs, int shift){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
		//cut child cluster from parent
		disconnect(parent, child);
		for(Node fs : fullSibs) disconnect(parent, fs);
		
		//shift
		if(shift==-1){//case1: shift down child & its full sib cluster
			clearVisit();
			List<Node> cluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}	
		}
		else{ //case2 : shift up parent cluster
			clearVisit();
			List<Node> cluster = parent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() + 1);
			}	
		}
		
		

		//connect child and parent
		connect(hs, child);
		for(Node fs : fullSibs) connect(hs, fs);
		
		clean(parent);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}

	
	public void POtoHS(Node child, Node parent, int targetSex, List<Node> fullSibs, int shift){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
			
		//cut child from parent
		disconnect(parent, child);
		for(Node fs : fullSibs) disconnect(parent, fs);
		
		
		//shift cluster
		if(shift==1){ //case1 : shift child cluster up
			clearVisit();
			//parent.setNumVisit(1);
			List<Node> shiftCluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : shiftCluster){
				i.setDepth(i.getDepth() + 1);
			}	
		}
		else{//case2 : shift parent cluster down
			clearVisit();
			List<Node> cluster = parent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}		
		}
		
		
		//get future parent
		Node gp = parent.getParentWithSex(targetSex);
		
		//make new node
		if(gp==null){
			gp = makeNewNode(parent.getDepth()+1, targetSex);
			connect(gp, parent);
		}
		
		//connect 
		connect(gp, child);
		for(Node fs : fullSibs) connect(gp, fs);
		
		clean(parent);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}
	
	
	public void HStoGP(Node child, Node parent, Node hs, List<Node> fullSibs, int targetSex, int childShift, int hsShift){

		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		clearVisit();
		parent.setNumVisit(1);
		List<Node> shiftCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + childShift);
		}	
		
		clearVisit();
		child.setNumVisit(1);
		for(Node fs : fullSibs) fs.setNumVisit(1);
		shiftCluster = hs.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + hsShift);
		}
		
		
		//cut child from parent
		disconnect(parent, child);
		for(Node fs : fullSibs) disconnect(parent, fs);
				
		
		//make ghost parent for child
		Node ghostParent = makeNewNode(child.getDepth()+1, targetSex);
		
		//connect ghost parents
		connect(ghostParent, child);
		for(Node fs : fullSibs) connect(ghostParent, fs);
		connect(hs, ghostParent);
		

		clean(parent);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}

	
	public void GPtoHS(Node child, Node parent, Node gp, List<Node> fullSibs, int childShift, int gpShift){
		
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
	
		

		//shift cluster;
		clearVisit();
		parent.setNumVisit(1);
		List<Node> shiftCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + childShift);
		}	
		
		clearVisit();
		parent.setNumVisit(1);
		shiftCluster = gp.getConnectedNodes(new ArrayList<Node>());
		for(Node i : shiftCluster){
			i.setDepth(i.getDepth() + gpShift);
		}	
		
		//disconnect
		disconnect(parent, child);
		for(Node fs : fullSibs) disconnect(parent, fs);
		
		int targetSex = parent.getSex();
		
		//get new parent
		Node newParent = gp.getParentWithSex(targetSex);
		
		//make new node
		if(newParent==null){
			newParent = makeNewNode(gp.getDepth()+1, targetSex);
			connect(newParent, gp);
		}
		
		//connect 
		connect(newParent, child);
		for(Node fs : fullSibs) connect(newParent, fs);
		
		
		clean(parent);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
	}
	
	
	
	
	public void contract(Node parent, Node child, int shift){

		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];

		//disconnect child cluster
		this.disconnect(parent, child);
			
		//shift child cluster up
		if(shift==1){
			clearVisit();
			List<Node> cluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() + 1);
			}
		}
		else{
			clearVisit();
			List<Node> cluster = parent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}
		}
		
		//new parents and children
		for(Node p : parent.getParents()){
			this.connect(p, child);
		}
		for(Node c : parent.getChildren()){
			this.connect(child, c);
		}
	
		
		//clean
		clean(parent);

		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
	}
	
	
	
	public void stretch(Node child, int shift){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		Node newParent = this.makeNewNode(child.getDepth(), child.getSex());
		

		//connect new parent, disconnect child cluster
		for(Node p : child.getParents()){
			connect(p, newParent);
			p.removeChild(child);
		}
		child.getParents().clear();
		
			
		if(shift==-1){ //case 1 : shift child cluster down
			clearVisit();
			List<Node> cluster = child.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() - 1);
			}
		}
		else{//case 2 : shift parent cluster up
			clearVisit();
			List<Node> cluster = newParent.getConnectedNodes(new ArrayList<Node>());
			for(Node i : cluster){
				i.setDepth(i.getDepth() + 1);
			}
		}

		//connect new parent to child
		this.connect(newParent, child);
		

		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
		
	}
	
	
	
	//both are sampled NEW: at least one of them is sampled; no node is deleted
	public void swapAncDesc(Node anc, Node desc){
		
		//get cluster
		clearVisit();
		List<Node> ped = desc.getConnectedSampledNodes(new ArrayList<Node>());


		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);	
		this.logLikelihood[curr] -= prior[curr];
		
		if(anc.getDepth() == desc.getDepth()+1){
			switchParentChild(anc,desc);
		}
		
		else{
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
		}
		
		
	
		//update adj matrix 
		for(Node ind : ped){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
	
	
		
	}
	
	

	
	public void uncle2nephew(Node uncle, Node nephew, Node fs, int targetSex, int uncleShift, int nephewShift){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = uncle.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
		//cut child from parents
		for(Node x : fs.getParents()){
			disconnect(x, uncle);
		}
		
		
		
		//shift clusters
		clearVisit();
		List<Node> shiftCluster = uncle.getConnectedNodes(new ArrayList<Node>());
		for(Node x : shiftCluster) x.setDepth(x.getDepth() + uncleShift);
		
		clearVisit();
		shiftCluster = nephew.getConnectedNodes(new ArrayList<Node>());
		for(Node x : shiftCluster) x.setDepth(x.getDepth() + nephewShift);
		
		
		//make ghost parent for child
		Node ghostParent = makeNewNode(uncle.getDepth()+1, targetSex);
		connect(ghostParent, uncle);
		
		//make ghost parents and nephew fullSibs
		makeFullSibs(nephew, ghostParent);
		
		//clean gp
		Node gp1 = fs.getParents().get(0);
		Node gp2 = fs.getParents().get(1);
		clean(gp1);
		clean(gp2);
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
	}
	
	
	
	public void nephew2uncle(Node child, Node gp, int childShift, int gpShift){
		
		//cluster containing child
		clearVisit();
		List<Node> ped = child.getConnectedSampledNodes(new ArrayList<Node>());
		
		//subtract likelihood for old cluster
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);
		//this.logLikelihood[curr] -= getSingletonProb();
		this.logLikelihood[curr] -= prior[curr];
		
		
		//cut child from parents
		Node parent = child.getParents().get(0);
		disconnect(parent, child);



		//shift clusters
		clearVisit();
		List<Node> shiftCluster = child.getConnectedNodes(new ArrayList<Node>());
		for(Node x : shiftCluster) x.setDepth(x.getDepth() + childShift);
		
		clearVisit();
		shiftCluster = gp.getConnectedNodes(new ArrayList<Node>());
		for(Node x : shiftCluster) x.setDepth(x.getDepth() + gpShift);

		
		//make child and gp fullSibs
		makeFullSibs(gp, child);
		
		//clean parent
		clean(parent);
		
		
		
		//update adj matrix
		for(Node ind : ped){
			updateAdjMat(ind);
		}
		
		
		//add new likelihood
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		//updateNumSingletons();
		//this.logLikelihood[curr] += getSingletonProb();
		updatePrior();
		this.logLikelihood[curr] += prior[curr];
		
		
	}
	
	
	
	
	
	////////////////////////   PRIOR   ////////////////////////
	//prior: random mating model
	public void updatePrior(){	
		
		//uniform sample for effective pop
		effectivePop = rGen.nextInt(maxN - minN) + minN;
		
		//clear everything
		clearVisit();
		for(int i=0; i<nNodes.length; i++) nNodes[i] = 0;
		for(int i=0; i<nNodes.length; i++) nUnsampledDads[i] = 0;
		for(int i=0; i<nNodes.length; i++) nUnsampledMoms[i] = 0;
		for(int i=0; i<nNodes.length; i++) nSampledDads[i] = 0;
		for(int i=0; i<nNodes.length; i++) nSampledMoms[i] = 0;
		
		
		// Count relevant quantities
		// for every cluster
		for(int i=0; i<numIndiv; i++){
			
			Node node = nodes.get(curr).get(i);
			
			if(node.getNumVisit() > 0) continue;
			
			List<Node> ped = node.getConnectedNodes(new ArrayList<Node>());

			//for every node in cluster
			for(Node x : ped){
				
				int xDepth = x.getDepth();
				
				nNodes[xDepth]++;
				
				if(x.sampled){
					
					if(x.getSex()==0) nSampledMoms[xDepth]++;
					else nSampledDads[xDepth]++;
					
				}
				
				else{
					
					if(x.getSex()==0) nUnsampledMoms[xDepth]++;
					else nUnsampledDads[xDepth]++;
				}
				
				
				//increment unlabeled parents even if they're not actually represented
				if(x.getParents().size()==0){
					countGhostNode(xDepth+1, 0);
					countGhostNode(xDepth+1, 1);
				}
				else if(x.getParents().size()==1){
					int ghostParentSex = (x.getParents().get(0).getSex()+1) % 2;
					countGhostNode(xDepth+1, ghostParentSex);
				}
				
				x.setNumVisit(1);
				
				
			}
			
			
		}
		
		//TODO testing
		updateNumSingletons();
		prior[curr] = computePrior(effectivePop) + sampleDepthPrior(nNodes) + getSingletonProb();
	
		
	}
	
	
	
	public double computePrior(int popSize){	
	
		//compute prior
		double fa = 0d;
		double ma = 0d;
		
		//for every generation
		for(int i=0; i<maxDepth; i++){
		
			//shared terms: 1/N^n
			double oneOverNn = -nNodes[i]*log[popSize/2];
			
			//father probability
			int start = popSize/2 - nSampledDads[i+1];
			int end = start - nUnsampledDads[i+1] + 1;
			
			//TODO testing
			if(end<0) return Double.NEGATIVE_INFINITY;
			
			for(int j=start; j>=end; j--)
				fa += log[j];	

			fa += oneOverNn;
			
			
			//mother probability
			start = popSize/2 - nSampledMoms[i+1];
			end = start - nUnsampledMoms[i+1] + 1;
			
			//TODO testing
			if(end<0) return Double.NEGATIVE_INFINITY;
			
			for(int j=start; j>=end; j--)
				ma += log[j];
			
			ma += oneOverNn;
	
		}
		
		
		return fa + ma + sampleDepthPrior(nNodes);
	
		
	}
	
	
	private double sampleDepthPrior(int[] nNodes){
		
		double toReturn = 0d;
		
		// P(max) = 1/totalUnits, P(max-1) = 2/totalUnits, ..., P(0)=(maxDepth+1)/totalUnits
		for(int i=0; i<nNodes.length; i++){

			toReturn += log[nNodes[i]] * (log[maxDepth+1-i] - log[totalUnits]);
			
			//toReturn += (log[nNodes[i]] - (i+1)*log[2]);
			
		}

		return toReturn;
		
	}
	
	
	private void countGhostNode(int currDepth, int sex){
		
		if(currDepth > maxDepth) return;
		
		nNodes[currDepth]++;
		
		if(sex==0) nUnsampledMoms[currDepth]++;
		else nUnsampledDads[currDepth]++;
		
		//recurse
		countGhostNode(currDepth+1, 0);
		countGhostNode(currDepth+1, 1);
		
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
		
		return beta * (nSingletons[curr]*logLambda - lambda - logFact[this.nSingletons[curr]]);
		
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
		
		//return toReturn + getSingletonProb();
		return toReturn;
		
		
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
	
	
	
	public double totalLikelihood(){
		
		double toReturn = 0d;
		
		clearVisit();
		
		for(int i=0; i<numIndiv; i++){
			
			Node node = nodes.get(curr).get(i);
			
			if(node.getNumVisit() > 0) continue;
			
			List<Node> ped = node.getConnectedSampledNodes(new ArrayList<Node>());
			
			toReturn += likelihoodLocalPedigree(ped);
			
			
		}
		
		//return toReturn + getSingletonProb();
		updatePrior();
		return toReturn + prior[curr];
		
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
		prior[copy] = prior[curr];
		
		
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
	
	
	
	public void dfs(Node node){
		
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
	
	//returns full siblings child has with the given sex, excluding itself
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
			
			if(sib.getSex()==targetSex && fullSibs(child, sib))
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
 		
 		/*
 		//pairwise error
 		//copy current relationships
 		Path[][] rel = new Path[numIndiv][numIndiv];
 		for(int j=0; j<numIndiv; j++){
 			for(int k=j+1; k<numIndiv; k++){
 				Path old = relationships[curr][j][k];
 				rel[j][k] = new Path(old.getUp(), old.getDown(), old.getNumVisit());
 			}
 		}
		for(int k=0; k<numIndiv; k++){
			Node ind = nodes.get(curr).get(k);
			updateAdjMat(ind);
		}
 		for(int j=0; j<numIndiv; j++){
 			for(int k=j+1; k<numIndiv; k++){
 				Path old = rel[j][k];
 				Path newRel = relationships[curr][j][k];
 				if(old.getUp()!=newRel.getUp() || old.getDown()!=newRel.getDown() || old.getNumVisit()!=newRel.getNumVisit()){
 					System.out.println("Pairwise error!");
 					return false;
 				}
 			}
 		}
 		*/
		
		
 		
 		
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
	public boolean sexLocked(Node node){

		if(node.sampled) return true;
		
		for(Node x : node.getChildren()){
			if(x.getParents().size()==2) return true;
		}
		
		return false;
		
	}
	
	
	//get spouses, including itself
	public List<Node> getSpouses(Node node, List<Node> spouses){

		if(node.getNumVisit()>0) return spouses;	
		
		//add this spouse
		spouses.add(node);
		node.setNumVisit(1);

		//recurse on neighboring parents
		for(Node c : node.getChildren()){
			
			if(c.getNumVisit()>0) continue;
			c.setNumVisit(1);
			
			for(Node p : c.getParents()){

				//spouse already visited
				if(p.getNumVisit()>0) continue;
				
				
				//new spouse
				spouses.addAll(getSpouses(p, spouses));

				
				
				
			}
			
		}
		
		return spouses;
		
	}
 
 	

}
	 
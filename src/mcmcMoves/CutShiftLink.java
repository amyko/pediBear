package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class CutShiftLink extends Move {

	
	public CutShiftLink(String name, double moveProb) {
		super(name, moveProb);
	}

	//for cut
	private int[] iDepthToCount = new int[maxDepth];
	private int[] jDepthToCount = new int[maxDepth];

	
	//for link
	private Node donor;
	private Node recipient;
	private boolean specialMerge;
	private boolean mergingFormsFullSibs;
	int nCuttableNode;
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//////////// CUT ///////////////
		
		//get a random child
		Node child = currPedigree.getRandomNode();

		//reject if not exactly 1 parent
		if(child.getParents().size()!=1) return REJECT;
		
		//parent to cut from
		Node parent = child.getParents().get(0);
		if(isSplitNode(parent)) //reject if parent is a splitNode
			return REJECT;


		//determine if the child has full siblings
		boolean hasFullSib = false; //because only has one parent

		
		//save current pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new  via split
		double oldToNewCut = 0d;
		Node highestNode = currPedigree.getHighestNode(parent);
		int targetDepth = highestNode.getDepth();

		if(highestNode.getChildren().size() > 1){ //split prob of the highest node
			int symm = (!highestNode.sampled && highestNode.getParents().size()==0) ? 1 : 0;
			oldToNewCut += (1+symm) * getPowersOfHalf2(highestNode.getChildren().size()) * moveProbs.get("splitLink");
		}
		
		//cut 
		int nBefore = currPedigree.getNActiveNodes();
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.cut(child, parent, hasFullSib);
		Node iPrime = currPedigree.clean(child);
		Node jPrime = currPedigree.clean(parent);
		int nAfter = currPedigree.getNActiveNodes();
		
		
		//old to new via cut
		oldToNewCut += (1 + nBefore - nAfter) * .5 * moveProbs.get("cutLink");
		oldToNewCut = getLogChooseOne(nBefore) + Math.log(oldToNewCut);
		

		//new to old
		iDepthToCount = currPedigree.getDepthToCount(iPrime, iDepthToCount);
		jDepthToCount = currPedigree.getDepthToCount(jPrime, jDepthToCount);
		
	
		double outerSum = 0d;
		double innerSum = 0d;
		for(int l1=0; l1<=iPrime.getDepth(); l1++){
			
			if(iDepthToCount[l1]==0) continue;
			
			innerSum = 0d;
			
			for(int l2=0; l2<=jPrime.getDepth(); l2++){
				
				if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth  - Math.max(l1,l2) - l1 - l2);
			}
			
			outerSum += iDepthToCount[l1] * innerSum;
		}
		
		double newToOldCut =  getLogChooseTwo(nAfter) + Math.log(outerSum);


		

		
		////////////// SHIFT //////////////////////////
		//get the cluster the node belongs to
		List<Node> cluster = new ArrayList<Node>();
		currPedigree.clearVisit();
		iPrime.getConnectedNodes(cluster);
		
		
		
		//get the lowest level of the cluster
		int oldLowestLevel = getLowestLevel(currPedigree, cluster);
		
		
		//pick the new lowest level
		int k = geometricDist(currPedigree.rGen);
		int newLowestLevel = k - 1;
		int offset = newLowestLevel - oldLowestLevel;
		
		//reject if highest node goes over max depth or no change
		int highestLevel = getHighestLevel(currPedigree, cluster);
		if(highestLevel + offset > currPedigree.maxDepth)
			return REJECT;
		

		
		//shift cluster
		currPedigree.shiftCluster(cluster, offset);
		
		
		//hastings ratio
		double oldToNewShift = Math.log(getPowersOfHalf(k));
		double newToOldShift = Math.log(getPowersOfHalf(oldLowestLevel+1));

		
		
		
		////////////////////// LINK /////////////////////
		//choose nodes i and j
		Node i = iPrime;
		Node j = currPedigree.getRandomNode();
		
		
		if(i.getIndex()==j.getIndex()) return REJECT;

		//choose target depth
		k = geometricDist(currPedigree.rGen);
		targetDepth = k - 1 + Math.max(i.getDepth(), j.getDepth());
		
		
		if(targetDepth > currPedigree.maxDepth || (i.getDepth()==j.getDepth() && i.getDepth()==targetDepth)){ 
			reverseMove(currPedigree);
			return REJECT;
		}

		
		
		//determine sex
		int targetSex;
		if(j.getDepth()==targetDepth){ //if j is at targetDepth, it determines the sex
			targetSex = j.getSex();
		}
		else if(i.getDepth()==targetDepth){
			targetSex = i.getSex();
		}
		else{//choose randomly
			targetSex = currPedigree.rGen.nextDouble() <.5 ? 0 : 1;
		}
		
		
		//take a random path to targetDepth-1
		nCuttableNode = 0;
		Node[] iCluster = getRandomPathAncestor(currPedigree, i, targetDepth, targetSex);
		Node[] jCluster = getRandomPathAncestor(currPedigree, j, targetDepth, targetSex);
		Node iAnc = iCluster[0];
		Node jAnc = jCluster[0];
		iPrime = iCluster[1];
		jPrime = jCluster[1];
		
		if(nCuttableNode==0){
			System.out.println("NO cuttable node");
		}
		
		//reject if both merging nodes are sampled, or both have ancestors, or they're the same node
		if((iAnc==jAnc) || (iAnc.sampled && jAnc.sampled) || (iAnc.getParents().size()>0 && jAnc.getParents().size()>0)){
			reverseMove(currPedigree);
			return REJECT;
		}
		if(createsIllegalCycle(currPedigree, iAnc, jAnc)){
			reverseMove(currPedigree);
			return REJECT;
		}
		
		
		
		
		//assign donor & recipient; recipient is sampled or has parents
		if(iAnc.sampled){
			donor = jAnc;
			recipient = iAnc;
		}
		else{
			donor = iAnc;
			recipient = jAnc;
		}


		//for later //TODO
		specialMerge = recipient.sampled && donor.getParents().size() > 0;
		if(specialMerge){
			reverseMove(currPedigree);
			return REJECT;
		}

		
		//reject bad cases
		/*
		if(violatesAgeConstraints(currPedigree, donor, recipient)){
			reverseMove(currPedigree);
			return REJECT;
		}
		*/

		

		//old to new via link
		iDepthToCount = currPedigree.getDepthToCount(iPrime, iDepthToCount);
		jDepthToCount = currPedigree.getDepthToCount(jPrime, jDepthToCount);
		
		outerSum = 0d;
		innerSum = 0d;
		for(int l1=0; l1<=iPrime.getDepth(); l1++){
			
			if(iDepthToCount[l1]==0) continue;
			
			innerSum = 0d;
			for(int l2=0; l2<=jPrime.getDepth(); l2++){
				
				if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth  - Math.max(l1,l2) - l1 - l2);
			}
			outerSum += iDepthToCount[l1] * innerSum;
		}

		double oldToNewLink = getLogChooseTwo(nBefore) + Math.log(outerSum);
		
		
		//merge
		currPedigree.merge(donor, recipient, mergingFormsFullSibs);
		currPedigree.clean(donor);
		nAfter = currPedigree.getNActiveNodes(); //add number of nodes created; subtract donor node (which will be deleted)
	

		//new to old via cut/split
		double newToOldLink = 0d;
		double cutProb = 0d;
	
		cutProb += nCuttableNode * moveProbs.get("cutShiftLink");


		newToOldLink = getLogChooseOne(nAfter) + Math.log(cutProb);
		
			
		double oldToNew = oldToNewCut + oldToNewLink + oldToNewShift;
		double newToOld = newToOldCut + newToOldLink + newToOldShift;
		

		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, oldToNew, newToOld, heat);

	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}
	
	

	
	//returns true if removing this node splits the pedigree into two
	private boolean isSplitNode(Node node){
		
		if(node.sampled || node.getParents().size()==0 || node.getChildren().size() > 1) return false;
		if(node.getParents().size()==2) return true;
		
		//ghost and has 1 parent; recurse on parent
		return isSplitNode(node.getParents().get(0));
		
	}
	

	

	
	
	public Node[] getRandomPathAncestor(Pedigree currPedigree, Node node, int targetDepth, int targetSex){
		
		if(node.getDepth() > targetDepth) throw new RuntimeException("starting node has depth greater than target depth");

		Node currNode = node;
		Node lastExistingNode = node;
		if(node.getDepth() == targetDepth){
			if(node.getSex()!=targetSex) 
				throw new RuntimeException("node's sex at targetDepth does not match targetSex");
			else{
				currNode = node;
				lastExistingNode = node;
			}
		}
		
		else{
			//go up to target depth
			int sex;
			while(currNode.getDepth() < targetDepth){
	
				if(currNode.getDepth()==targetDepth-1){ //target sex 
					sex = targetSex;
				}
				else{//randomly choose mom or dad
					sex = currPedigree.rGen.nextInt(2);
				}
				
				Node parent = currNode.getParentWithSex(sex);
				
				//if there was no parent with the given sex, make one
				if(parent==null){
					parent = currPedigree.makeNewNode(currNode.getDepth() + 1, sex);
					currNode.addParent(parent);
					parent.addChild(currNode);
					nCuttableNode++;
				}
				else{
					lastExistingNode = parent;
				}
				
				currNode = parent;
						
			}
			
		}
		
		
		return new Node[]{currNode, lastExistingNode};

		
	}
	
	private boolean createsIllegalCycle(Pedigree currPedigree, Node node1, Node node2){		

		//creates cycles
		//there's a cycle if donor can reach recipient before merging
		//technically, two merging nodes can have a common parent, but ignore this case for now
		currPedigree.performDFS(node1);
		
		if(node2.getNumVisit() > 1){
			return true;
		}
		
		if(node2.getNumVisit() == 1){ //okay only if merging creates FS
			
			mergingFormsFullSibs = formsFullSibs(node1.getChildren(), node2.getChildren());
			
			if(mergingFormsFullSibs)
				return false;
			else
				return true;
		}		
		
		else{
			mergingFormsFullSibs = false;
			return false;
		}
		
	}
	

	private boolean violatesAgeConstraints(Pedigree currPedigree, Node donor, Node recipient){

		//check maxD < minR
		Node maxDonorDesc = currPedigree.getDescendantWithMaxAge(donor);
		Node minRecipientAnc = null;
		
		if(recipient.sampled){
			minRecipientAnc = recipient;
		}
		else{
			minRecipientAnc = currPedigree.getAncestorWithMinAge(recipient);
		}
		
		if(maxDonorDesc!=null && minRecipientAnc!=null && (minRecipientAnc.getAge() <= maxDonorDesc.getAge())){
			return true;
		}

		

		//check minD > maxR
		Node minDonorAnc = currPedigree.getAncestorWithMinAge(donor);
		Node maxRecipientDesc = null;
		
		if(recipient.sampled){
			maxRecipientDesc = recipient;
		}
		else{
			maxRecipientDesc = currPedigree.getDescendantWithMaxAge(recipient);
		}
		
		if(minDonorAnc!=null && maxRecipientDesc!=null && (minDonorAnc.getAge() <= maxRecipientDesc.getAge())){
			return true;
		}
		
		
		return false;
		
		
	}

	
	
	//returns true if merging form a FS
	private boolean formsFullSibs(List<Node> donorChildren, List<Node> recipientChildren){
		
		for(Node i : donorChildren){
			if(i.getParents().size()!=2) continue;
			for(Node j : recipientChildren){
				if(j.getParents().size()!=2) continue;
				
				if(i==j) continue;
				
				//already share a parent
				if(i.getParents().get(0) == j.getParents().get(0) || i.getParents().get(1) == j.getParents().get(1) || i.getParents().get(0) == j.getParents().get(1) || i.getParents().get(1) == j.getParents().get(0)) 
					return true;
			}
		}
		
		return false;
	}
	


	
	private int getLowestLevel(Pedigree currPedigree, List<Node> cluster){
		
		int lowest = currPedigree.maxDepth;
		
		for(Node i : cluster){
			if(i.getDepth() < lowest){
				lowest = i.getDepth();
			}
			if(lowest==0) break;
		}
		
		
		return lowest;
		
		
	}

	
	private int getHighestLevel(Pedigree currPedigree, List<Node> cluster){
		
		int highest = 0;
		
		for(Node i : cluster){
			if(i.getDepth() > highest){
				highest = i.getDepth();
			}
			if(highest==currPedigree.maxDepth) break;
		}
		
		
		return highest;
		
		
	}
	

}
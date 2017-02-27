package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class SplitLink extends Move {//WORKS

	
	public SplitLink(String name, double moveProb) {
		super(name, moveProb);
	}

	//for split
	private boolean hasFullSib;
	private int[] iDepthToCount = new int[maxDepth];
	private int[] jDepthToCount = new int[maxDepth];
	private List<Node> splitChildren = new ArrayList<Node>();
	private List<Node> stayChildren = new ArrayList<Node>();
	
	//for link
	private Node donor;
	private Node recipient;
	int nCuttableNode;
	private boolean specialMerge;
	private boolean mergingFormsFullSibs;
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//////////// SPLIT //////////////
		
		//get a random node to split
		Node parent = currPedigree.getRandomNode();
		
		//reject if node has < 2 children
		int nChildren = parent.getChildren().size();
		if(nChildren < 2)
			return REJECT;
		
	
		//randomly assign children to a clone (1 to nChildren-1); choose between [0,2^n - 2] and add 1
		int powerSetInd = currPedigree.rGen.nextInt((int) getPowersOfTwo(nChildren)-1) + 1;

		
		//no change if all children are assigned to clone and the original parent is a ghost
		if(!parent.sampled && powerSetInd==(int) getPowersOfTwo(nChildren)-1) 
			return REJECT;
		
		
		splitChildren.clear();
		stayChildren.clear();
		List<Node> children = parent.getChildren();
	    for (int i = 0; i < nChildren; i++) {
	        if((1 << i & powerSetInd) != 0){ //if corresponding bit=1, add child
	        	splitChildren.add(children.get(i));
	        }  
	        else{
	        	stayChildren.add(children.get(i));
	        }
	    }
	    
	    
	    
	    //TODO
	    //reject if split creates illegal cycles
	    if(createsIllegalCycle(currPedigree, splitChildren, stayChildren, parent)) 
	    	return REJECT;
	    	
	    



	    //check if any of the children between split and stay form full sibs
	    hasFullSib = false;
	    for(Node c1 : splitChildren){
	    	if(c1.getParents().size()!=2) continue;
	    	for(Node c2 : stayChildren){
	    		hasFullSib = currPedigree.fullSibs(c1, c2);
	    		if(hasFullSib) break;
	    	}
	    	if(hasFullSib) break;
	    }
	    
	    
	    
		//old to new
		//via split
	    int symmetric = (!parent.sampled && parent.getParents().size()==0) ? 1 : 0;
		double oldToNewSplit = (1+symmetric) * getPowersOfHalf2(nChildren) * moveProbs.get("splitLink");
	    
		//save current pedigree
		currPedigree.copyCurrPedigree();
	    double prevLogLikelihood = currPedigree.getLogLikelihood();
	    int nBefore = currPedigree.getNActiveNodes();
	    int targetDepth = parent.getDepth();
	    
	    //split
	    Node splitParent = currPedigree.makeNewNode(parent.getDepth(), parent.getSex());
	    currPedigree.split(parent, splitParent, splitChildren, hasFullSib);
	    Node iPrime = currPedigree.clean(parent);
	    Node jPrime = currPedigree.clean(splitParent);	    
	    int nAfter = currPedigree.getNActiveNodes();
	    

		//via cut
		oldToNewSplit += .5 * (nBefore + 1 - nAfter) * moveProbs.get("cutLink");
		oldToNewSplit = getLogChooseOne(nBefore) + Math.log(oldToNewSplit);
		
		
		//new to old	
		iDepthToCount = currPedigree.getDepthToCount(iPrime, iDepthToCount);
		jDepthToCount = currPedigree.getDepthToCount(jPrime, jDepthToCount);
		
		double outerSum = 0d;
		double innerSum = 0d;
		for(int l1=0; l1<=iPrime.getDepth(); l1++){
			innerSum = 0d;
			for(int l2=0; l2<=jPrime.getDepth(); l2++){
				
				//if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth  - Math.max(l1,l2) - l1 - l2);
			}
			outerSum += iDepthToCount[l1] * innerSum;
		}
		double newToOldSplit =  getLogChooseTwo(nAfter) + Math.log(outerSum);

		
		
		
		////////////////////// LINK /////////////////////
		//choose nodes i and j
		nBefore = currPedigree.getNActiveNodes();
		Node[] nodes = currPedigree.getNRandomNodes(2);
		Node i = nodes[0];
		Node j = nodes[1];
		//int offset = currPedigree.rGen.nextInt(currPedigree.getNActiveNodes()-1) + 1;

		
		//choose target depth
		int k = geometricDist(currPedigree.rGen);
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
		if((iAnc.sampled || iAnc.getParents().size() > 0) && !jAnc.sampled){
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
			innerSum = 0d;
			for(int l2=0; l2<=jPrime.getDepth(); l2++){
				
				//if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth  - Math.max(l1,l2) - l1 - l2);
			}
			outerSum += iDepthToCount[l1] * innerSum;
		}

		double oldToNewLink = getLogChooseTwo(nBefore) + Math.log(outerSum);

		
		//merge
		currPedigree.merge(donor, recipient, mergingFormsFullSibs);
		currPedigree.clean(donor);
		nAfter = currPedigree.getNActiveNodes(); //add number of nodes created; subtract donor node (which will be deleted)

		//new to old
		double newToOldLink = 0d;
		double cutProb = 0d;
		double splitProb = 0d;

		//via cut
		cutProb = nCuttableNode * .5 * moveProbs.get("cutLink");

		//via split
		if(recipient.getChildren().size()>=2){ 
			int symm = !recipient.sampled && recipient.getParents().size()==0 ? 1 : 0;
			splitProb = (1+symm) * getPowersOfHalf2(recipient.getChildren().size()) * moveProbs.get("splitLink");		
		}
		


		newToOldLink = getLogChooseOne(nAfter) + Math.log(cutProb + splitProb);
		
		
		
		double oldToNew = oldToNewSplit + oldToNewLink;
		double newToOld = newToOldSplit + newToOldLink;
		

		
		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLogLikelihood, heat);
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
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
	
	/*
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
	*/

	
	
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

	
	private boolean createsIllegalCycle(Pedigree currPedigree, List<Node> splitChildren, List<Node> stayChildren, Node parent){		

	    for(Node x : splitChildren){
	    	
	    	if(currPedigree.getFullSibs(x).size()==0) continue;
	    	
	    	for(Node y : splitChildren){
	    		
	    		if(x.getIndex()==y.getIndex()) continue;
	    		
	    		if(currPedigree.getFullSibs(y).size()==0) continue;
	    		
	    		if(!currPedigree.fullSibs(x,y)) return true;
	    		
	    	}
	    		
	    		
	    	
	    }
	    
	    return false;

	}
	

}
package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class Cut extends Move {//WORKS

	
	public Cut(String name, double moveProb) {
		super(name, moveProb);
	}


	private int[] iDepthToCount = new int[maxDepth];
	private int[] jDepthToCount = new int[maxDepth];
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		

		//get a random child
		Node child = currPedigree.getRandomNode();
		
		
		//choose mom or dad to cut from
		int parentSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		Node parent = child.getParentWithSex(parentSex);
		if(parent==null || isSplitNode(parent)) //reject if parent not available or parent is a splitNode
			return REJECT;
		
		


		//determine if the child has full siblings; if so, cutting doens't split the pedigree
		boolean hasFullSib = hasFullSib(currPedigree, child);


		//copy current pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new  via split
		double oldToNew = 0d;
		Node highestNode = currPedigree.getHighestNode(parent);
		int targetDepth = highestNode.getDepth();


		if(highestNode.getChildren().size() > 1){ //split prob of the highest node
			int symm = (!highestNode.sampled && highestNode.getParents().size()==0) ? 1 : 0;
			oldToNew += (1+symm) * getPowersOfHalf2(highestNode.getChildren().size()) * moveProbs.get("split");
		}

		
		//cut 
		int nBefore = currPedigree.getNActiveNodes();
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.cut(child, parent, hasFullSib);
		Node iPrime = currPedigree.clean(child);
		Node jPrime = currPedigree.clean(parent);
		int nAfter = currPedigree.getNActiveNodes();
		
		
		//old to new via cut
		oldToNew += (1 + nBefore - nAfter) * .5 * moveProbs.get("cut");
		oldToNew = getLogChooseOne(nBefore) + Math.log(oldToNew);
		
		
		//new to old via link
		currPedigree.getDepthToCount(iPrime, iDepthToCount);
		currPedigree.getDepthToCount(jPrime, jDepthToCount);
		
		double innerSum = 0d;
		double outerSum = 0d;
		for(int l1=0; l1 <= iPrime.getDepth(); l1++){
			
			if(iDepthToCount[l1]==0) continue;
			
			innerSum = 0d;
			for(int l2=0; l2 <= jPrime.getDepth(); l2++){
				
				if(l1==targetDepth && l2==targetDepth){
					continue;
				}
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth-Math.max(l1,l2)-l1-l2);
			}
			outerSum += iDepthToCount[l1] * innerSum;
		}
		double newToOld =  getLogChooseTwo(nAfter) + Math.log(moveProbs.get("link") * outerSum);
		
		
		
		//accept ratio
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
	
	
	//returns true if child has at leaset one full sibling
	private boolean hasFullSib(Pedigree currPedigree, Node child){
		
		if(child.getParents().size() != 2){
			return false;
		}
		else{
			for(Node candidate : child.getParents().get(0).getChildren()){ //mom's children
				if(currPedigree.fullSibs(child, candidate))
					return true;
			}

				
			for(Node candidate : child.getParents().get(1).getChildren()){ //dad's children
				if(currPedigree.fullSibs(child, candidate))
					return true;
			}
			
			
		}
		
		return false;
	}

	

	

}


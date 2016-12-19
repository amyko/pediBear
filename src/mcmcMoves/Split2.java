package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//Handles the special case where the chosen node is sampled and has parents & children
//This move replaces the sampled node with a ghost, and splits off the sampled node from the pedigree. 
//Children are randomly assigned to the sampled node, but must not assign all of them. 

public class Split2 extends Move {//TODO make this move more efficient (can we combine this with split?)

	
	public Split2(String name, double moveProb) {
		super(name, moveProb);
	}

	//this is stuff to store each move so that if it gets accepted it can then be performed
	private Node parent;
	private Node stayParent; //ghost node
	private boolean hasFullSib;
	private int[] iDepthToCount = new int[maxDepth];
	private int[] jDepthToCount = new int[maxDepth];
	private List<Node> stayChildren = new ArrayList<Node>();
	private List<Node> splitChildren = new ArrayList<Node>();
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
				
		//get a random sampled node
		parent = currPedigree.getRandomSampledNode();

		
		//reject if node has no parents or has < 2 children
		int nChildren = parent.getChildren().size();
		if(nChildren == 0 || parent.getParents().size()==0)
			return REJECT;
		
	
		//randomly assign children to stay (1 to 2^n-1); choose between [0,2^n - 2] and add 1
		int powerSetInd = currPedigree.rGen.nextInt((int) getPowersOfTwo(nChildren)-1) + 1;
		

		
		stayChildren.clear();
		splitChildren.clear();
		List<Node> children = parent.getChildren();
	    for (int i = 0; i < nChildren; i++) {
	        if((1 << i & powerSetInd) != 0){ //if corresponding bit=1, add child
	        	stayChildren.add(children.get(i));
	        }
	        else{
	        	splitChildren.add(children.get(i));
	        }
	    }
	    
	    
	    //reject if no children are left behind
	    if(stayChildren.size()==0) return REJECT;
	    
		
	    
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
	    
	
	    

		//split
	    double prevLkhd = currPedigree.getLogLikelihood();
	    int nBefore = currPedigree.getNActiveNodes();
	    stayParent = currPedigree.makeNewNode(parent.getDepth(), parent.getSex());
	    currPedigree.split2(parent, stayParent, stayChildren, hasFullSib);
	    
		
	    
		//old to new via split2
	    double oldToNew = 0d;
		oldToNew = getPowersOfHalf2(nChildren) * moveProbs.get("split2");
		oldToNew = getLogChooseOne(nBefore) + Math.log(oldToNew);
		
		
		//new to old via link
		Node iPrime = parent;
		Node jPrime = stayParent;
		currPedigree.getDepthToCount(iPrime, iDepthToCount);
		currPedigree.getDepthToCount(jPrime, jDepthToCount);
		int targetDepth = parent.getDepth();
		
		
		double newToOld = 0d;
		double innerSum = 0d;
		double outerSum = 0d;
		
		for(int l1=0; l1<=iPrime.getDepth(); l1++){
			
			if(iDepthToCount[l1]==0) continue;
			
			innerSum = 0d;
			for(int l2=0; l2<=jPrime.getDepth(); l2++){
				
				if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth  - Math.max(l1,l2) - l1 - l2);
			}
			
			
			outerSum += iDepthToCount[l1] * innerSum;
		
		
		}
		int nAfter = currPedigree.getNActiveNodes();
		newToOld = getLogChooseTwo(nAfter) + Math.log(outerSum * moveProbs.get("link"));

		//return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLogLikelihood, heat);
		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, oldToNew, newToOld, heat);
		
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {

		
		currPedigree.merge(stayParent, parent, hasFullSib);
		currPedigree.clean(stayParent);


		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}

	
}
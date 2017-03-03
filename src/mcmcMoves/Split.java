package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;


public class Split extends Move {

	
	public Split(String name, double moveProb) {
		super(name, moveProb);
	}

	//this is stuff to store each move so that if it gets accepted it can then be performed
	private boolean hasFullSib;
	private List<Node> splitChildren = new ArrayList<Node>();
	private List<Node> stayChildren = new ArrayList<Node>();
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
				
		//get a random node to split
		Node parent = currPedigree.getRandomNode();
		

		//reject if node has < 2 children
		int nChildren = parent.getChildren().size();
		if(nChildren < 2)
			return REJECT;
		
	
		//randomly assign children to a clone (1 to 2^n-1); choose between [0,2^n - 2] and add 1
		int powerSetInd = currPedigree.rGen.nextInt((int) getPowersOfTwo(nChildren)-1) + 1;
		
		//no change if all children are assigned to clone
		if(!parent.sampled && powerSetInd==(int) getPowersOfTwo(nChildren)-1) 
			return REJECT;
		
		
		//assign children to split or stay
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
	    
	    
	    //reject if creates illegal cycle
	    if(createsIllegalCycle(currPedigree, splitChildren, stayChildren, parent))
	    	return REJECT;
	    
	    
    
		//save current pedigree
		currPedigree.copyCurrPedigree();
	    double prevLogLikelihood = currPedigree.getLogLikelihood();

	    //split
	    Node splitParent = currPedigree.makeNewNode(parent.getDepth(), parent.getSex());
	    currPedigree.split(parent, splitParent, splitChildren, hasFullSib);
	    currPedigree.clean(parent);
	    currPedigree.clean(splitParent);	    


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

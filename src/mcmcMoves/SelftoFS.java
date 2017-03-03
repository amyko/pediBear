package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class SelftoFS extends Move{

	
	private List<Node> splitChildren = new ArrayList<Node>();
	private List<Node> stayChildren = new ArrayList<Node>();
	
	public SelftoFS(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a node to split
		Node parent = currPedigree.getRandomNode();

		
		//reject if no children
		int nChildren = parent.getChildren().size();
		if(nChildren==0) return REJECT;
		
		
		//depth constraint
		if(parent.getDepth()==currPedigree.maxDepth) return REJECT;
		
		
		//randomly assign children to a clone (1 to 2^n-1); choose between [0,2^n - 2] and add 1
		int powerSetInd = currPedigree.rGen.nextInt((int) getPowersOfTwo(nChildren)-1) + 1;
		
		
		//reject if all children are assigned to clone and the parent would be deleted
		if(!parent.sampled && powerSetInd==(int) getPowersOfTwo(nChildren)-1){
			return REJECT;
		}

		
		
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
		
		
	
	    //reject if fullSibs exist between stay and split clusters
	    for(Node x : splitChildren){
	    	
	    	for(Node y : stayChildren)
	    		if(currPedigree.fullSibs(x, y)) return REJECT;
	    	
	    }
	    
	    
		
		//copy pedigree
		currPedigree.copyCurrPedigree();

		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.selfToFS(parent, splitChildren);
		

		//accept ratio
		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, heat);

		
		
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}

}
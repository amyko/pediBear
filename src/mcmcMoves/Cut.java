package mcmcMoves;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class Cut extends Move {//WORKS

	
	public Cut(String name, double moveProb) {
		super(name, moveProb);
	}

	
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
		//cut 
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.cut(child, parent, hasFullSib);
		currPedigree.clean(child);
		currPedigree.clean(parent);
		

		
		//accept ratio
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


package Unused;


import mcmcMoves.Move;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class Compress extends Move {

	
	public Compress(double moveProb) {
		super(moveProb);
	}

	//this is stuff to store each move so that if it gets accepted it can then be performed
	private Node grandParent;
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//get a node
		grandParent = currPedigree.getRandomNode();
		
		
		//reject if depth is too low or has parents
		if(grandParent.getDepth() < 2 || grandParent.getParents().size() > 0)
			return REJECT;
		
		//reject if any of the children is sampled or has two parents or has different sex from grandparent or has more than 1 grandchildren per child
		for(Node c : grandParent.getChildren()){
			if(c.sampled || c.getParents().size()==2 || c.getSex()!=grandParent.getSex() || c.getChildren().size()!=1) 
				return REJECT;
		}

		


		
		//compress
		int nBefore = currPedigree.getNActiveNodes();
		int nAfter = nBefore - grandParent.getChildren().size();
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		
		currPedigree.compress(grandParent);
		
		
		//old to new
		double oldToNew = logChooseOne[nBefore] + Math.log(moveProbs[7]);

		double newToOld =  logChooseOne[nAfter] + Math.log(moveProbs[6]);
		

		double acceptRatio = heat * (currPedigree.getLogLikelihood() - prevLogLikelihood) + newToOld - oldToNew;
		
		
		if(acceptRatio > 0){
			return 1;
		}
		else{
			return Math.exp(acceptRatio);
		}
	}
	
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
	}
	
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {
		
		currPedigree.stretch(grandParent);
		
	}
	
	

	

}

	
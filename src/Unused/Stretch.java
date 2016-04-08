package Unused;


import mcmcMoves.Move;
import dataStructures.Pedigree;
import dataStructures.Node;

//To check if the edge is cuttable, check if the parent is 1) ghost, 2) has two parents, 3) has no other children. Recurse up. 

public class Stretch extends Move {//works

	
	public Stretch(double moveProb) {
		super(moveProb);
	}

	//this is stuff to store each move so that if it gets accepted it can then be performed
	private Node parent;
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//get a node
		parent = currPedigree.getRandomNode();
		
		
		//reject if highest node, has fewer than 2 children, or has parents
		if(parent.getDepth()==currPedigree.depth || parent.getChildren().size()<2 || parent.getParents().size()>0)
			return REJECT;
		
		//reject if any of the children from full sibs
		for(Node c1 : parent.getChildren()){
			for(Node c2 : parent.getChildren()){
				if(currPedigree.fullSibs(c1, c2))
					return REJECT;
			}
		}
		
		//reject if stretching would violate age constraints
		if(violatesAgeConstraints(currPedigree, parent))
			return REJECT;
		
		

		
		//stretch
		int nBefore = currPedigree.getNActiveNodes();
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		
		currPedigree.stretch(parent);
		
		
		//old to new
		double oldToNew = logChooseOne[nBefore] + Math.log(moveProbs[6]);

		double newToOld =  logChooseOne[currPedigree.getNActiveNodes()] + Math.log(moveProbs[7]);
		

		double acceptRatio = heat * (currPedigree.getLogLikelihood() - prevLogLikelihood) + newToOld - oldToNew;

		if(acceptRatio > 0){
			return 1;
		}
		else{
			return Math.exp(acceptRatio);
		}
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {
		
		currPedigree.compress(parent);
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}
	

	
	private boolean violatesAgeConstraints(Pedigree currPedigree, Node parent){

		if(!parent.sampled) return false;
		
		//check maxD < minR
		Node maxDonorDesc = currPedigree.getDescendantWithMaxAge(parent);
		

		if(maxDonorDesc!=null && ((parent.getAge() - maxDonorDesc.getAge()) < (parent.getDepth() - maxDonorDesc.getDepth() + 1) * currPedigree.genTime)){
			return true;
		}


		return false;
		
		
	}

	

	

	

}

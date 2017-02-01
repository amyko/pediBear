package mcmcMoves;


import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;


// 1) swap two sampled nodes

public class OPtoPO extends Move {

	
	public OPtoPO(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose a sampled node
		Node lowerNode = currPedigree.getRandomSampledNode();

		
		//depth
		if(lowerNode.getDepth() + 2 > currPedigree.maxDepth)
			return REJECT;
		
		
		int numParents = lowerNode.getParents().size();
		if(numParents!=1)
			return REJECT;
		
		if(lowerNode.getChildren().size()!=0)
			return REJECT;
		
		
		//choose a parent
		Node middleNode = lowerNode.getParents().get(0);
		
		
		//get grand parents to switch with
		Node upperNode = middleNode.getParentWithSex(lowerNode.getSex());
		
		//to prevent overlap with swapAncDesc
		if(upperNode!=null && upperNode.sampled)
			return REJECT;
		


		if(isSplitNode(middleNode, upperNode)) 
			return REJECT;	
		
		
		/*
		//TODO testing
		boolean testing = false;
		if(lowerNode.getIndex()==0 && upperNode!=null && upperNode.getChildren().size()==3){
			
			Path rel = currPedigree.getRelationships()[0][8];
			if(rel.getUp()==1 && rel.getDown()==3 && rel.getNumVisit()==1){
				
				
				rel = currPedigree.getRelationships()[0][5];
				
				if(rel.getUp()==2 && rel.getDown()==3 && rel.getNumVisit()==1){
					
					System.out.println();
					
					testing = true;
				}
				
				
			}
			
			
		}
		*/
		
		
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		

		
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.op2po(lowerNode, middleNode, upperNode);
		
		
		boolean deleteMiddleNode = false;
		if(middleNode.getSex()==-1)//middle node deleted
			deleteMiddleNode = true;
		
		double oldToNew = Math.log(moveProbs.get("OPtoPO"));
		

		double newToOld = 0d;
		int numChildren = lowerNode.getChildren().size();
		
		if(deleteMiddleNode){ //make a new ghost lower node
			newToOld = 1.0 / (1+numChildren) * .5; //choose sex
		}
		else{
			newToOld = numChildren/(1.0+numChildren) * (1d/numChildren);
		}
		
		newToOld = Math.log(newToOld * moveProbs.get("POtoOP"));

		
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
	private boolean isSplitNode(Node parent, Node gp){
		
		//parent would be ghost node after cutting child
		if(!parent.sampled && parent.getChildren().size() < 2){
			
			//2 parents
			if(parent.getParents().size()==2) 
				return true;
			
			//effectively 2 parents
			if(parent.getParents().size()==1 && gp==null)
				return true;
			
		}
		
			
		
		return false;
		
	}
	
	
}
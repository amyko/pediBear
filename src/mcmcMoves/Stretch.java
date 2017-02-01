package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;

//TODO this is not exactly the reverse move of contract; need to selectively shift down descendants of child

public class Stretch extends Move{ //WORKS; special merge not tested

	
	public Stretch(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomNode();
				
		
		//reject if child has no parents
		if(child.getParents().size()==0)
			return REJECT;
		
		
		//prevents overlap with swap down
		if(child.getChildren().size()==0)
			return REJECT;
		

		
		//reject if shifting down the child cluster violates depth constraint
		currPedigree.clearVisit();
		for(Node p : child.getParents()) p.setNumVisit(1);
		int minCluDepth = currPedigree.getMinDepth(child);
		if(minCluDepth < 1) return REJECT;
		
		/*
		boolean testing = false;
		if(child.getIndex()==3){
			
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
		
		
		//old to new
		double oldToNew = Math.log(moveProbs.get("stretch")) + getLogChooseOne(currPedigree.getNActiveNodes());
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		

		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.stretch(child);
		
		
		//newToOld
		double newToOld = Math.log(moveProbs.get("contract")) + getLogChooseOne(currPedigree.getNActiveNodes());

		
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

}

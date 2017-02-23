package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;

//TODO this is not exactly the reverse move of contract; need to selectively shift down descendants of child

public class Stretch extends Move{ //WORKS; special merge not tested

	
	public Stretch(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomSampledNode();
				
		
		//reject if child has no parents
		if(child.getParents().size()==0)
			return REJECT;
		
		/*
		//prevents overlap with swap down
		if(child.getChildren().size()==0)
			return REJECT;
		*/
		
		//choose 1) shift child cluster down or 2) shift parent cluster up
		int shift = currPedigree.rGen.nextDouble() < .5 ? -1 : 1;

		
		//depth constraint
		if(shift==-1){ //case 1: shift child cluster down
			currPedigree.clearVisit();
			for(Node p : child.getParents()) p.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(child);
			if(minDepth == 0) return REJECT;
		}
		else{ //case 2: shift parent cluster up
			
			currPedigree.clearVisit();
			for(Node x : child.getChildren()) x.setNumVisit(1); 
			
			int maxDepth = currPedigree.getMaxDepth(child);
			if(maxDepth==currPedigree.maxDepth) return REJECT;
			
			
		}

		
		//old to new
		double oldToNew = Math.log(moveProbs.get("stretch"));
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		

		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.stretch(child, shift);
		
		
		//newToOld
		double newToOld = Math.log(moveProbs.get("contract"));

		
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

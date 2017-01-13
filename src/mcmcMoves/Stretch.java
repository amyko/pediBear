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
		Node child = currPedigree.getRandomNode();
				
		/*
		//reject if child is a ghost AND doesn't have any descendants
		if(!child.sampled && child.getChildren().size()==0)//this should never happen
			return REJECT;
		*/
		
		
		//reject if child has no parents
		if(child.getParents().size()==0)
			return REJECT;
		
		
		//prevents overlap with swap down
		if(child.sampled && child.getChildren().size() == 0)
			return REJECT;
		
		//reject if shifting down the child cluster violates depth constraint
		currPedigree.clearVisit();
		int minCluDepth = currPedigree.getMinDepth(child);
		if(minCluDepth < 1) return REJECT;
		
		
		//old to new
		double oldToNew = Math.log(moveProbs.get("stretch")) + getLogChooseOne(currPedigree.getNActiveNodes());
		
		/*
		//via swap down
		if(child.sampled && child.getChildren().size()==0){
			oldToNew += moveProbs.get("swapDown");
		}
		*/
		
		oldToNew = Math.log(oldToNew);
		
		
		
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

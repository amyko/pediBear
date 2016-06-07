package mcmcMoves;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class POtoFS extends Move{

	
	public POtoFS(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		//get parent
		Node parent = child.getParents().get(0);
		
		//reject if parent is ghost and has no children
		if(!parent.sampled && parent.getChildren().size()==1)
			return REJECT;
		
		//with prob=.5, shift parent cluster down or shift child cluster up
		boolean goUp = currPedigree.rGen.nextDouble() < .5 ? true : false;
		
		//reject if going up violates depth constraint
		if(goUp){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int maxDepth = currPedigree.getMaxDepth(child);
			if(maxDepth+2 > currPedigree.maxDepth)
				return REJECT;
		}
		else{
			
			currPedigree.clearVisit();
			child.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(parent);
			if(minDepth==0)
				return REJECT;
			
		}
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) +  Math.log(moveProbs.get("POtoFS"));
		
		//modify pedigree
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.POtoFS(child, parent, goUp);
		
		
		//new to old
		int nSibs = currPedigree.getFullSibs(child).size();
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nSibs) +Math.log(moveProbs.get("FStoPO"));
		

		
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

}
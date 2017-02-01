package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class POtoHS extends Move{

	
	public POtoHS(String name, double moveProb) {
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
		
		//reject if parent is ghost and has no other children
		if(!parent.sampled && parent.getChildren().size()<2)
			return REJECT;
		
			
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		int maxDepth = currPedigree.getMaxDepth(child);
		if(maxDepth+1 > currPedigree.maxDepth || child.getDepth()+2 > currPedigree.maxDepth)
			return REJECT;


		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0 : 1;
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) +  Math.log(.5*moveProbs.get("POtoHS"));
		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.POtoHS(child, parent, targetSex);
		

		
		//new to old
		int nSibs = child.getParents().get(0).getChildren().size() - 1;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nSibs) + Math.log(moveProbs.get("HStoPO"));
		

		
		//accept ratio
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
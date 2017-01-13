package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; make a ghost parent; link ghost parent to two parents

public class CutOneLinkTwo extends Move{

	
	public CutOneLinkTwo(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child's grand parent depth > max depth
		if(child.getDepth() + 2 > currPedigree.maxDepth)
			return REJECT;
		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;
		

		//reject parent is a ghost and has only one child
		Node parent = child.getParents().get(0);
		if(!parent.sampled && parent.getChildren().size() == 1)
			return REJECT;
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("cutOneLinkTwo"));
		
		//cut and link
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.cutOneLinkTwo(child);
		
		
		//new to old
		int nCandidate = currPedigree.getFullSibsWithTargetSex(parent, (parent.getSex()+1)%2).size();
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nCandidate) + Math.log(moveProbs.get("cutTwoLinkOne"));
		
	

		
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

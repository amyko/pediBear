package mcmcMoves;



import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class CutTwoLinkOne extends Move{

	
	public CutTwoLinkOne(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a parent to cut
		Node parent = currPedigree.getRandomNode();

		
		//reject if parent does not have exactly 2 parents or 1 child
		if(parent.sampled || parent.getParents().size()!=2 || parent.getChildren().size()!=1)
			return REJECT;
		
		//reject if child has more than one parent
		Node child = parent.getChildren().get(0);
		if(child.getParents().size()!=1)
			return REJECT;
		
		int newParentSex = (parent.getSex() + 1) % 2;
		List<Node> newParentCandidates = currPedigree.getFullSibsWithTargetSex(parent, newParentSex);

		
		//reject there's no potential new parent
		int nCandidate = newParentCandidates.size();
		if(nCandidate==0)
			return REJECT;
		
		
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nCandidate) + Math.log(moveProbs.get("cutTwoLinkOne"));
		
		
		//cut and link
		double prevLkhd = currPedigree.getLogLikelihood();
		Node newParent = newParentCandidates.get(currPedigree.rGen.nextInt(nCandidate));
		currPedigree.cutTwoLinkOne(parent, newParent);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("cutOneLinkTwo"));
		

		
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

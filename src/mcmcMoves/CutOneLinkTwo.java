package mcmcMoves;

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
		//reject if the child doesn't have any half/full siblings
		if(child.getParents().get(0).getChildren().size() < 2)
			return REJECT;

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("cutOneLinkTwo"));
		Node parent = child.getParents().get(0);
		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.cutOneLinkTwo(child);
		
		
		//new to old
		int nCandidate = currPedigree.getFullSibs(parent, parent.getSex()).size() + 1;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nCandidate) + Math.log(moveProbs.get("cutTwoLinkOne"));
		

		
		//accept ratio
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

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}

}

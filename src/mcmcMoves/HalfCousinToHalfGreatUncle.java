
package mcmcMoves;


import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class HalfCousinToHalfGreatUncle extends Move{

	
	
	public HalfCousinToHalfGreatUncle(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size() != 1)
			return REJECT;
		
		//reject if parent doesn't have exactly 1 child or exactly one parent
		Node parent = child.getParents().get(0);
		if(parent.sampled || parent.getChildren().size()!=1 || parent.getParents().size()!=1)
			return REJECT;

		//reject if parent doesn't have any sibs
		Node gp = parent.getParents().get(0);
		if(gp.getChildren().size() < 2)
			return REJECT;
		
		
		//depth constraint
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		int maxDepth = currPedigree.getMaxDepth(child);
		if(maxDepth+3 > currPedigree.maxDepth)
			return REJECT;

					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("halfCousinToHalfGreatUncle"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.halfCousinToHalfGreatUncle(child, gp);
		
		
		//new to old
		int nSibs = 0;
		for(Node i : child.getParents().get(0).getChildren()){
			if(i!=child && i.getChildren().size() > 0) nSibs++;
		}
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nSibs) + +getLogChooseOne(2) + Math.log(moveProbs.get("halfGreatUncleToHalfCousin"));
		

		
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
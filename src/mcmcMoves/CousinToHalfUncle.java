
package mcmcMoves;

import java.util.List;

import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class CousinToHalfUncle extends Move{

	
	
	public CousinToHalfUncle(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	
		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		//reject if child doesn't have exactly one parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		Node parent = child.getParents().get(0);
		
		
		//reject if parent doesn't have two parents or exactly one child
		if(parent.sampled || parent.getChildren().size()!=1 || parent.getParents().size()!=2)
			return REJECT;
		
		//reject if parent doesn't have any full sibs
		List<Node> sibs = currPedigree.getFullSibs(parent);
		if(sibs.size()==0)
			return REJECT;
			
		
		//depth constraint
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		int maxDepth = currPedigree.getMaxDepth(child);
		if(maxDepth+1 > currPedigree.maxDepth)
			return REJECT;

		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + lnTwo + Math.log(moveProbs.get("cousinToHalfUncle"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.cousinToHalfUncle(child);
		
		
		//new to old
		int nHalfSibs = child.getParents().get(0).getChildren().size() - 1;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nHalfSibs) + Math.log(moveProbs.get("halfUncleToCousin"));
		

		
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

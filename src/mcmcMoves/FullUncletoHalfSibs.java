
package mcmcMoves;

import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class FullUncletoHalfSibs extends Move{

	
	public FullUncletoHalfSibs(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	
		//choose a child
		Node child = currPedigree.getRandomNode();

		//reject if child doesn't have exactly one parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		Node parent = child.getParents().get(0);
		
	
		//reject if parent doesn't have two parents or exactly one child
		if(parent.sampled || parent.getChildren().size()!=1 || parent.getParents().size()!=2)
			return REJECT;
		
		
		
		//reject if child doesn't have any full uncle
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
		
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("fullUncleToHalfSibs"));
		
		//cut and link
		double prevLkhd = currPedigree.getLogLikelihood();
		Node uncle = sibs.get(currPedigree.rGen.nextInt(sibs.size()));
		int targetSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		currPedigree.fullUncletoHalfSibs(child, uncle, targetSex);
		
		int nFS = currPedigree.getFullSibs(uncle).size();
		int nHS = child.getParents().get(0).getChildren().size() - 1;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nHS) + Math.log((nFS+1) * moveProbs.get("halfSibsToFullUncle"));
		
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

package mcmcMoves;

import java.util.List;

import mcmc.SimulatedAnnealing;
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
		

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		Node uncle = sibs.get(currPedigree.rGen.nextInt(sibs.size()));
		currPedigree.fullUncletoHalfSibs(child, uncle);
		

		
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
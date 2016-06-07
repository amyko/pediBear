
package mcmcMoves;

import java.util.List;

import mcmc.SimulatedAnnealing;
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
		
		//reject if child doesn't have any full cousins
		List<Node> sibs = currPedigree.getFullSibs(parent);
		int nCousinParents = 0;
		for(Node s : sibs){
			if(s.getChildren().size() > 0) nCousinParents++;
		}
		if(nCousinParents==0)
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
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("cousinToHalfUncle"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.cousinToHalfUncle(child);
		
		
		//new to old
		int nSibs = 0;
		for(Node i : child.getParents().get(0).getChildren()){
			if(i!=child && i.getChildren().size() > 0) nSibs++;
		}
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nSibs) + getLogChooseOne(2) + Math.log(nCousinParents * moveProbs.get("halfUncleToCousin"));
		

		
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

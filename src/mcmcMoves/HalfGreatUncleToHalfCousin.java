
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift down child cluster; 

public class HalfGreatUncleToHalfCousin extends Move{

	List<Node> sibs = new ArrayList<Node>();
	
	public HalfGreatUncleToHalfCousin(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;

		
		//reject if child does not have any half siblings with children
		Node parent = child.getParents().get(0);
		sibs.clear();
		for(Node sib : parent.getChildren()){
			if(sib != child && sib.getChildren().size()>0)
				sibs.add(sib);
		}
		
		if(sibs.size()==0)
			return REJECT;
		
			
		//reject if moving down child two levels goes below depth=0
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		if(currPedigree.getMinDepth(child) - 2 < 0){
			return REJECT;
		}
		
		
		//choose sib
		Node newGP = sibs.get(currPedigree.rGen.nextInt(sibs.size()));

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(sibs.size()) + getLogChooseOne(2) + Math.log(moveProbs.get("halfGreatUncleToHalfCousin"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();

		
		currPedigree.halfGreatUncleToHalfCousin(child,newGP);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("halfCousinToHalfGreatUncle"));
		

		
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


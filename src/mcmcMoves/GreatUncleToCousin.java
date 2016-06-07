
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift down child cluster; 

public class GreatUncleToCousin extends Move{

	
	private List<Node> goodFS  = new ArrayList<Node>();
	
	
	public GreatUncleToCousin(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have 2 parents
		if(child.getParents().size()!=2)
			return REJECT;

		//reject if no full sibs with children
		List<Node> fullSibs = currPedigree.getFullSibs(child);
		goodFS.clear();
		for(Node i : fullSibs){
			if(i.getChildren().size() > 0)
				goodFS.add(i);
		}
		if(goodFS.size()==0)
			return REJECT;
		

		//reject if moving down child two levels goes below depth=0
		currPedigree.clearVisit();
		for(Node k : child.getParents())
			k.setNumVisit(1);
		if(currPedigree.getMinDepth(child) - 2 < 0){
			return REJECT;
		}
		
		
		
		//choose sib
		Node sib = goodFS.get(currPedigree.rGen.nextInt(goodFS.size()));
		
		//chose sib child
		List<Node> sibChildren = sib.getChildren();
		Node sibChild = sibChildren.get(currPedigree.rGen.nextInt(sibChildren.size()));

		//get number of equivalent choices of sib child
		int nSibChoice = 1 + currPedigree.getFullSibs(sibChild).size();
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(goodFS.size()) + getLogChooseOne(sibChildren.size()) + getLogChooseOne(2) + Math.log(nSibChoice * moveProbs.get("greatUncleToCousin"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.greatUncleToCousin(child, sibChild);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("cousinToGreatUncle"));
		

		
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



package mcmcMoves;

import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class POtoHS extends Move{

	
	public POtoHS(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		//reject if no parent
		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;

		
		//get parent of targetSex
		int targetSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		Node parent = child.getParentWithSex(targetSex);
		
		if(parent==null) return REJECT;
		
		//full siblings
		List<Node> fullSibs = currPedigree.getFullSibs(child);
		
		//reject if parent would be deleted
		int leftOverChildren = parent.getChildren().size()-fullSibs.size()-1;
		if(!parent.sampled && leftOverChildren==0)
			return REJECT;


		
			
		//choose 1) shift child cluster up or 2) shift parent cluster down
		int shift = currPedigree.rGen.nextDouble() < .5 ? -1 : 1;
	
		
		
		//depth constraint
		if(shift==1){ //case 1: shift child cluster up
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int maxDepth = Math.max(currPedigree.getMaxDepth(child), child.getDepth()+1);
			if(maxDepth == currPedigree.maxDepth) return REJECT;
		}
		else{ //case 2: shift parent cluster down
			currPedigree.clearVisit();
			child.setNumVisit(1);
			for(Node x : fullSibs) x.setNumVisit(1);

			int minDepth = currPedigree.getMinDepth(parent);
			if(minDepth==0) return REJECT;
	
		} 
		
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();

		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.POtoHS(child, parent, targetSex, fullSibs, shift);
	
		
		//accept ratio
		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, heat);
		

		
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
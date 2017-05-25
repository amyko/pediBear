package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class HStoPO extends Move{

	private List<Node> halfSibs = new ArrayList<Node>();
	private List<Node> fullSibs = new ArrayList<Node>();
	
	public HStoPO(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		
		//reject if no parents
		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;
		

		//choose targetSex for parent
		int targetSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		
		
		//reject if there's no parent with targetSex
		Node parent = child.getParentWithSex(targetSex);
		if(parent==null) return REJECT;
		
		
		//get half sibs with targetSex
		halfSibs.clear();
		fullSibs.clear();
		for(Node x : parent.getChildren()){
			if(x.getIndex()==child.getIndex()) continue; //skip itself
			if(currPedigree.fullSibs(x, child)) fullSibs.add(x); //fs
			else if(x.getSex()==targetSex) halfSibs.add(x); //hs
		}

		
		if(halfSibs.size()==0) return REJECT;
		
		Node hs = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));

		
		//choose 1) shift child cluster down or 2) shift parent cluster up
		int shift = currPedigree.rGen.nextDouble() < .5 ? -1 : 1;
		


		
		//depth constraint
		if(shift==-1){ //case 1: shift child cluster down
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(child);
			if(minDepth == 0) return REJECT;
		}
		else{ //case 2: shift parent cluster up
			currPedigree.clearVisit();
			child.setNumVisit(1);
			for(Node x : fullSibs) x.setNumVisit(1);		
			
			 //don't count parent node that will be deleted
			if(!parent.sampled && (parent.getNumEdges() - fullSibs.size() -1) < 2) parent.setNumVisit(1);
			
			int maxDepth = currPedigree.getMaxDepth(hs);
			if(maxDepth==currPedigree.maxDepth) return REJECT;
			
		} 


		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new: choose targetSex, choose halfsib
		//symmetry if child has fullSibs: cancels out in reverse move
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(halfSibs.size()) + Math.log(moveProbs.get("HStoPO"));
		
		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.HStoPO(child, parent, hs, fullSibs, shift);
		
		
		//new to old: choose targetSex for parent, symmetry if child has fullSibs
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("POtoHS"));
		

		
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
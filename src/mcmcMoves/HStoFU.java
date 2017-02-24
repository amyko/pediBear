
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Pedigree;
import dataStructures.Node;
import mcmc.MCMCMC;


//cut child from parent and link to two grand parents

public class HStoFU extends Move{

	private List<Node> fullSibs = new ArrayList<Node>();
	private List<Node> halfSibs = new ArrayList<Node>();
	
	public HStoFU(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child
		Node child = currPedigree.getRandomNode();
		

		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;

		
		//choose parent to cut from
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0 : 1;
		Node parent = child.getParentWithSex(targetSex);
		if(parent==null) return REJECT;
		
		
		
		//reject if no half sibs
		halfSibs.clear();
		fullSibs.clear();
		for(Node x : parent.getChildren()){
			
			if(x.getIndex() == child.getIndex()) continue;
			
			if(currPedigree.fullSibs(x, child)) fullSibs.add(x);
			
			else halfSibs.add(x);
		}
	
		
		if(halfSibs.size()==0)
			return REJECT;
		
		
		//choose 1) shift child cluster down or 2) shift parent cluster up
		int shift = currPedigree.rGen.nextDouble() < .5 ? -1 : 1;
	
		
		//choose half sib
		Node halfSib = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));

		
		
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
			
			int maxDepth = Math.max(currPedigree.getMaxDepth(halfSib), halfSib.getDepth()+1);
			if(maxDepth==currPedigree.maxDepth) return REJECT;
			
		} 
		
		

					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		double prevLkhd = currPedigree.getLogLikelihood();
		
	
		
		//choose node, choose targetSex for parent, choose halfSib, choose shift direction
		//symmetry for fullSibs of child, fullSibs of HS; sex & direction cancel out
		int symm = (currPedigree.getFullSibs(halfSib).size()+1);
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(halfSibs.size()) +  Math.log(symm*moveProbs.get("HStoFU"));
		
	
		//modify
		currPedigree.HStoFU(child, parent, halfSib, targetSex, fullSibs, shift);
		
		//choose targetSex for parent, symmetry for fullSibs of child
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("FUtoHS"));
		
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
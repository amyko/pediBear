
package Unused;

import java.util.List;

import mcmc.MCMCMC;
import mcmcMoves.Move;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class FUtoHS extends Move{

	
	public FUtoHS(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	
		//choose a child
		Node child = currPedigree.getRandomNode();

		
		//choose parent
		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;
		

		//choose parent
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0 : 1;
		Node parent = child.getParentWithSex(targetSex);
		if(parent==null) return REJECT;
		
			
		List<Node> fullSibs = currPedigree.getFullSibs(child);
	
		//reject if parent is not a dummy node
		if(parent.sampled || parent.getParents().size()!=2 || parent.getChildren().size() != (fullSibs.size()+1))
			return REJECT;
	
		
		
		
		//reject if child doesn't have any full uncle
		List<Node> uncles = currPedigree.getFullSibs(parent);
		if(uncles.size()==0)
			return REJECT;
			
		
		//choose 1) shift child cluster down or 2) shift parent cluster up
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
			parent.setNumVisit(1);

			//don't count nodes that will be deleted
			for(Node x: parent.getParents()){
				if(!x.sampled && (x.getNumEdges() - 1) < 2) x.setNumVisit(1);	
			}
			
			int minDepth = currPedigree.getMinDepth(uncles.get(0));
			if(minDepth==0) return REJECT;
	
		} 
		
		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//modify
		double prevLkhd = currPedigree.getLogLikelihood();

		
		//choose targetSex for parent, symmetry for fullSibs of child
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("FUtoHS"));
		

		currPedigree.FUtoHS(child, parent, targetSex, fullSibs, shift);
	
		//new to old: targetSex, choose halfSib, direction
		//symmetry for fullSibs of child and fullSibs of hs
	
		int nHS = 0;
		for(Node x : child.getParentWithSex(targetSex).getChildren()){
			if(x.getIndex()==child.getIndex()) continue;
			if(currPedigree.fullSibs(x, child)) continue;
			nHS++;
		}
				
		int symm = currPedigree.getFullSibs(uncles.get(0)).size()+1;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nHS) + Math.log(symm*moveProbs.get("HStoFU"));
		
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
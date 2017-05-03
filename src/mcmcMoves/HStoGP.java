package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class HStoGP extends Move{

	private List<Node> halfSibs = new ArrayList<Node>();
	private List<Node> fullSibs = new ArrayList<Node>();
	
	public HStoGP(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		
		//reject if no parents
		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;
		
		//target sex for parent
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0 : 1;
		Node parent = child.getParentWithSex(targetSex);
		
		if(parent==null) return REJECT;
		
		
		//get half sib
		halfSibs.clear();
		fullSibs.clear();
		for(Node x : parent.getChildren()){
			if(x.getIndex()==child.getIndex()) continue;
			if(currPedigree.fullSibs(x, child)) fullSibs.add(x);
			else halfSibs.add(x);
		}
		
		if(halfSibs.size()==0) return REJECT;
		
		Node hs = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));
		
		
		//choose depth
		double p = currPedigree.rGen.nextDouble();
		int childShift = -1;
		int hsShift = 1;

		
		if(p < oneThird){
			childShift = -2;
			hsShift = 0;
		}
		else if(p < twoThirds){
			childShift = 0;
			hsShift = 2;
		}
		
		
		
		//depth constraint
		if(childShift!=0){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			if(currPedigree.getMinDepth(child) + childShift < 0)
				return REJECT;
		}
		
		if(hsShift!=0){
			currPedigree.clearVisit();
			child.setNumVisit(1);
			for(Node x : fullSibs) x.setNumVisit(1);
			//don't count if parent is to be deleted
			if(!parent.sampled && (parent.getNumEdges() - fullSibs.size() - 1) < 2)
				parent.setNumVisit(1);
			
			int maxDepth = currPedigree.getMaxDepth(hs);
			if(maxDepth + hsShift > currPedigree.maxDepth)
				return REJECT;
		}
		
		
		//old2new: choose node, choose parent targetSex, choose half sib, choose shift, symmetry for FS
		int symm = fullSibs.size() + 1;
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(halfSibs.size()) + Math.log(.5 * symm * moveProbs.get("hs2gp"));
		
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
	
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.HStoGP(child, parent, hs, fullSibs, targetSex, childShift, hsShift);
		

		//newToOld: choose node, choose targetSex for parent, shift direction, symm for nFS
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(.5 * symm * moveProbs.get("gp2hs"));
		
		
		
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
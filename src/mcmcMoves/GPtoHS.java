package mcmcMoves;

import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class GPtoHS extends Move{

	
	public GPtoHS(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		
		int parentSize = child.getParents().size();
		if(parentSize==0) return REJECT;
		
		//get parent that will be cut
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0 : 1;
		Node parent = child.getParentWithSex(targetSex);
		if(parent==null) return REJECT;
		
		
		//reject if parent is not a dummy node
		if(parent.sampled || parent.getParents().size()!=1) 
			return REJECT;
		
		Node gp = parent.getParents().get(0);

		
		//full siblings
		List<Node> fullSibs = currPedigree.getFullSibs(child);
		
		//parent not dummy node
		if(parent.getChildren().size() != (fullSibs.size()+1))
			return REJECT;
		
		
		//reject if gp is ghost and has no other children
		if(!gp.sampled && (gp.getChildren().size()) < 2) 
			return REJECT;
		
		
		//choose shift direction
		double p = currPedigree.rGen.nextDouble();
		int childShift = 1;
		int gpShift = -1;

		
		if(p < oneThird){
			childShift = 2;
			gpShift = 0;
		}
		else if(p < twoThirds){
			childShift = 0;
			gpShift = -2;
		}
		
		//depth
		if(childShift!=0){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int maxDepth = Math.max(currPedigree.getMaxDepth(child), child.getDepth()+1);
			if(maxDepth+childShift > currPedigree.maxDepth)
				return REJECT;
		}
		
		if(gpShift!=0){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(gp);
			if(minDepth + gpShift < 0)
				return REJECT;
		}

		
		
		//oldToNew: choose node, choose targetSex for parent, shift direction, symm for nFS
		int symm = fullSibs.size() + 1;
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(.5 * symm * moveProbs.get("gp2hs"));
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.GPtoHS(child, parent, gp, fullSibs, childShift, gpShift);

		
		//newToOld: choose node, choose parent targetSex, choose half sib, choose shift, symmetry for FS
		int nHS = 0;
		for(Node x : child.getParentWithSex(targetSex).getChildren()){
			if(x.getIndex()==child.getIndex()) continue;
			if(currPedigree.fullSibs(x, child)) continue;
			else nHS++;
		}
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nHS) + Math.log(.5 * symm * moveProbs.get("hs2gp"));
		
		
		
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
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class NephewtoUncle extends Move{

	
	List<Node> fullSibs = new ArrayList<Node>();
	
	public	NephewtoUncle(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a node
		Node child = currPedigree.getRandomNode();
		
		
		//exactly 1 parent
		if(child.getParents().size()!=1) return REJECT;
		
		//parent is a dummy node with exactly two parents
		Node parent = child.getParents().get(0);
		if(parent.sampled || parent.getParents().size()!=2 || parent.getChildren().size()!=1) return REJECT;
		
		//reject if no uncle
		List<Node> uncles = currPedigree.getFullSibs(parent);
		if(uncles.size()==0) return REJECT;
		
		
		//choose gp; will form full sibs with child later
		int targetSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		Node gp = parent.getParentWithSex(targetSex);

		
		//choose shift direction
		double p = currPedigree.rGen.nextDouble();
		int gpShift = -1;
		int childShift = 1;
		if(p < oneThird){
			gpShift = -2;
			childShift = 0;
		}
		else if(p < twoThirds){
			gpShift = 0;
			childShift = 2;
		}
		
		
		//depth constraint	
		if(childShift!=0){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int maxDepth = Math.max(currPedigree.getMaxDepth(child), child.getDepth()+1);
			if(maxDepth + childShift > currPedigree.maxDepth)
				return REJECT;
		}
		
		if(gpShift!=0){
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(uncles.get(0));
			if(minDepth + gpShift < 0)
				return REJECT;
		}
		
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//oldToNew: choose node, choose targetSex for gp, shift
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(.5 * moveProbs.get("uncle2nephew"));
		


		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.nephew2uncle(child, gp, childShift, gpShift);
		
		
		//newToOld: choose node, choose FS with children, choose nephew, shift, targetSex for ghostParent
		//symmetry if nephew had full siblings
		int nFS = 0;
		for(Node x : currPedigree.getFullSibs(child)){
			if(x.getChildren().size()>0) nFS++;
		}
		int nNephews = gp.getChildren().size();
		int symm = currPedigree.getFullSibs(uncles.get(0)).size() + 1;
		
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nFS) + getLogChooseOne(nNephews) + Math.log(symm * .5 * moveProbs.get("uncle2nephew"));
		
		
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
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class HStoPO extends Move{

	private List<Node> halfSibs = new ArrayList<Node>();
	
	public HStoPO(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		if(child.getParents().size()!=1)
			return REJECT;
		
		Node parent = child.getParents().get(0);
		
		//no half sibs
		int nChildren = parent.getChildren().size();
		if(nChildren < 2)
			return REJECT;
		
		//depth
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		if(currPedigree.getMinDepth(child)==0)
			return REJECT;
		
		
		//get half sib
		halfSibs.clear();
		for(Node x : parent.getChildren()){
			if(x.getIndex()==child.getIndex()) continue;
			halfSibs.add(x);
		}
		Node hs = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(halfSibs.size()) + Math.log(moveProbs.get("HStoPO"));
		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.HStoPO(child, hs);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(.5*moveProbs.get("POtoHS"));
		

		
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
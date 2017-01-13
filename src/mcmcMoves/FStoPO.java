package mcmcMoves;

import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class FStoPO extends Move{

	
	public FStoPO(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();

		//choose full sib
		List<Node> sibs = currPedigree.getFullSibs(child);

		//reject if child doesn't have any full sibs
		if(sibs.size()==0)
			return REJECT;
		
		//choose parent
		Node parent = sibs.get(currPedigree.rGen.nextInt(sibs.size()));
		
		//parent goes up or child goes down
		boolean goUp = currPedigree.rGen.nextDouble() < .5 ? true : false;
		
		//reject if going up violates depth constraint
		if(goUp){
			currPedigree.clearVisit();

			//don't count parents if they are going to be deleted
			for(Node p : parent.getParents()){
				if(!p.sampled && p.getNumEdges()<3)
					p.setNumVisit(1);
			}
			
			
			int maxDepth = currPedigree.getMaxDepth(parent);
			if(maxDepth+1 > currPedigree.maxDepth)
				return REJECT;
			
			
		}
		else{
			
			currPedigree.clearVisit();
			child.getParents().get(0).setNumVisit(1);
			child.getParents().get(1).setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(child);
			if(minDepth - 1 < 0)
				return REJECT;
			
		}

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(sibs.size()) + Math.log(moveProbs.get("FStoPO"));
		
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.FStoPO(child, parent, goUp);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("POtoFS"));
		

		
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
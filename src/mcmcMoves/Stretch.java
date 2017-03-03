package mcmcMoves;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


public class Stretch extends Move{ //WORKS; special merge not tested

	
	public Stretch(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomSampledNode();
				
		
		//reject if child has no parents
		if(child.getParents().size()==0)
			return REJECT;
		
		

		
		//determine whether this move requires HS2PO or contract as the reverse move
		if(child.getParents().size()==1){
			
			Node parent = child.getParents().get(0);
			
			if(!parent.sampled && parent.getSex()==child.getSex() && parent.getParents().size()==0 && parent.getChildren().size()>1){
				
				Node hs = null;
				for(Node x : parent.getChildren()){
					if(x.getIndex()==child.getIndex()) continue;
					hs = x;
					break;
				}
				
				if(currPedigree.getFullSibs(hs).size()+2 == parent.getChildren().size()){
					return REJECT;
				}

				
			}
			
		}
		
		

		
		//choose 1) shift child cluster down or 2) shift parent cluster up
		int shift = currPedigree.rGen.nextDouble() < .5 ? -1 : 1;

		
		//depth constraint
		if(shift==-1){ //case 1: shift child cluster down
			currPedigree.clearVisit();
			for(Node p : child.getParents()) p.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(child);
			if(minDepth == 0) return REJECT;
		}
		else{ //case 2: shift parent cluster up
			
			currPedigree.clearVisit();
			for(Node x : child.getChildren()) x.setNumVisit(1); 
			
			int maxDepth = currPedigree.getMaxDepth(child);
			if(maxDepth==currPedigree.maxDepth) return REJECT;
			
			
		}
		


		//copy pedigree
		currPedigree.copyCurrPedigree();
		

		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.stretch(child, shift);
		

		
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
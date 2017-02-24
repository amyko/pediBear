package mcmcMoves;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;

//TODO this is not exactly the reverse move of contract; need to selectively shift down descendants of child

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
		

		
		//old to new
		double oldToNew = Math.log(moveProbs.get("stretch"));


		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		

		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.stretch(child, shift);
		
		
		//newToOld
		double newToOld = Math.log(moveProbs.get("contract"));
		
		/*
		if(hsReverse){
			
			//currPedigree.printAdjMat();
			//System.out.println();
			
			int nHS = 0;
			int nFS = 0;
			int targetSex = child.getSex();
			for(Node x : hs.getParentWithSex(targetSex).getChildren()){
				if(x.getIndex()==hs.getIndex()) continue;
				if(currPedigree.fullSibs(x, hs)) nFS++;
				else if(x.getSex()==targetSex) nHS++;
			}
			


			//choose node, choose targetSex, choose halfSib, choose direction (cancels out in forward move)
			//symmetry for fullSibs of child
			newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nHS) + Math.log(.5 * (nFS+1) * moveProbs.get("HStoPO"));
		}
		
		
		else{
			newToOld = getLogChooseOne(currPedigree.numIndiv) + Math.log(moveProbs.get("contract"));
		}
		*/

		
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

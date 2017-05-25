package mcmcMoves;


import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;


public class Contract extends Move{ 

	
	public Contract(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomSampledNode();
		
		
		
		//reject if it doesn't have exactly one parent
		if(child.getParents().size()!=1) 
			return REJECT;

		
		Node parent = child.getParents().get(0);
		
		//reject if the parent is sampled or has wrong sex; must have grand parents
		if(parent.sampled || parent.getSex()!=child.getSex() || parent.getChildren().size()!=1 || parent.getParents().size()==0)
			return REJECT;
		

		
		//reject if confounds with POtoHS (i.e. both gp and parent would be deleted)
		if(parent.getParents().size()==1){
			
			Node gp = parent.getParents().get(0);
			
			//check if gp would be deleted by removing halfSib cluster
			if(!gp.sampled && gp.getSex()==child.getSex() && gp.getParents().size()==0 && gp.getChildren().size()>1){
				
				Node hs = null;
				for(Node x : gp.getChildren()){
					if(x.getIndex()==parent.getIndex()) continue;
					hs = x;
					if(currPedigree.getFullSibs(hs).size()+2 == gp.getChildren().size()){
						
						return REJECT;
						
						
					}
				}
				

						
				
			}
			
		}
		
		
		
		

		//choose 1) shift child cluster up or 2) shift parent cluster down
		int shift = currPedigree.rGen.nextDouble() < .5 ? 1 : -1;
		
		
		//depth constraint
		if(shift==1){ //case 1: shift child cluster up
			currPedigree.clearVisit();
			parent.setNumVisit(1);
			int maxDepth = currPedigree.getMaxDepth(child);
			if(maxDepth == currPedigree.maxDepth) return REJECT;
		}
		else{ //case 2: shift parent cluster down
			currPedigree.clearVisit();
			child.setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(parent);
			if(minDepth==0) return REJECT;
		}
		

		
		//oldToNew
		double oldToNew = Math.log(moveProbs.get("contract"));

		
		//record previous config
		currPedigree.copyCurrPedigree();
		
		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.contract(parent, child, shift);
		
		//newToOld
		double newToOld = Math.log(moveProbs.get("stretch"));

		
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
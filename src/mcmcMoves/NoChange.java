package mcmcMoves;


import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//chooses a ghost node and switches sex of its parent and everyone in the chain, if possible

public class NoChange extends Move {

	int prevN;
	double prevLkhd;
	double prevPrior;
	
	public NoChange(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		prevLkhd = currPedigree.getLogLikelihood();
		prevN = currPedigree.getEffectivePop();
		prevPrior = currPedigree.getPrior();
		
		currPedigree.noChange();
		

		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, 0,0, heat);
		
		
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {
			
		currPedigree.setLogLikelihood(prevLkhd);
		currPedigree.setEffectivePop(prevN);
		currPedigree.setPrior(prevPrior);

	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;	
		
	}
	
	
	//returns true if the parent's sex cannot be changed
	private boolean sexLocked(Node parent){

		parent.setNumVisit(parent.getNumVisit()+1);
		
		if(parent.sampled) 
			return true;
		
		//recurse on neighbor parents
		else{
		
			for(Node i : parent.getChildren()){
				
				if(i.getNumVisit() > 0) continue;
				i.setNumVisit(i.getNumVisit()+1);
				
				
				for(Node j : i.getParents()){
					
					if(j.getNumVisit() > 0) continue;
					
					if(sexLocked(j))
						return true;
					
				}
				
								
			}
			
			return false;
		}

		
	}
	
	

	
	
}
	


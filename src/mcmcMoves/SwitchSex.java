package mcmcMoves;


import dataStructures.Pedigree;
import dataStructures.Node;

//chooses a ghost node and switches sex of its parent and everyone in the chain, if possible

public class SwitchSex extends Move {

	
	public SwitchSex(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//get a random parent
		Node parent = currPedigree.getRandomNode();
		
		//reject if sex cannot be switched
		currPedigree.clearVisit();
		if(sexLocked(parent)) return REJECT;
			
		//switch sexes
		currPedigree.clearVisit();
		switchSex(parent);
		
		
		//always accept
		return 1;
		
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {
		
		return;

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
	
	
	//returns true if the parent's sex cannot be changed
	private void switchSex(Node parent){ //works
		
		if(parent.getNumVisit() > 0) return;
		
		//switch sex and mark visit
		parent.setSex((parent.getSex()+1) % 2);
		parent.setNumVisit(parent.getNumVisit()+1);
		
		//recurse on neighbor parents
		for(Node i : parent.getChildren()){
			
			if(i.getNumVisit() > 0) continue;
			i.setNumVisit(i.getNumVisit()+1);

			for(Node j : i.getParents()){
				switchSex(j);
			}
							
		}
	}
	
	
}
	


	
package mcmcMoves;

import mcmc.MCMCMC;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Node;
import dataStructures.Pedigree;

// 1) swap two sampled nodes

public class SwitchChildAncWithOppositeSex extends Move {

	
	List<Node> anc = new ArrayList<Node>();
	List<Node> candidates = new ArrayList<Node>();
	
	
	public SwitchChildAncWithOppositeSex(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
			
		// randomly choose a sampled node
		Node child = currPedigree.getRandomSampledNode();
		
		//reject if child is sex locked
		currPedigree.clearVisit();
		if(sexLocked(child, child.getIndex())) 
			return REJECT;
		
		
		//add to candidate list if the ancestor is sampled, is opposite sex from the child, and not sex locked
		anc.clear();
		candidates.clear();
		child.getAncestors(anc);
		int targetSex = (child.getSex()+1)%2;
		
		for(Node x : anc) {
			
			if(x.sampled && x.getSex()==targetSex) {
				
				currPedigree.clearVisit();
				if(!sexLocked(x, x.getIndex()))
					candidates.add(x);
				
			}
			
		}
		
		//reject if no candidates ancestors are left
		if(candidates.size()==0) return REJECT;
		
		
		//choose ancestor to switch with
		Node ancestor = candidates.get(currPedigree.rGen.nextInt(candidates.size()));
		

		//copy pedigree
		currPedigree.copyCurrPedigree();
	
		
		
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.switchChildAncWithOppositeSex(child, ancestor);
		

		//this move is symmetric
		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, 0d, 0d, heat);
		
	
	}
	
	
	
	//returns true if the parent's sex cannot be changed
	private boolean sexLocked(Node parent, int exception){

		//if sampled node is visited, then it's sex locked. EXCEPTION: the starting node doesn't count
		parent.setNumVisit(parent.getNumVisit()+1);
		
		if(parent.getIndex()!=exception && parent.sampled) 
			return true;
		
		//recurse on neighbor parents
		else{
		
			for(Node i : parent.getChildren()){
				
				if(i.getNumVisit() > 0) continue;
				i.setNumVisit(i.getNumVisit()+1);
				
				
				for(Node j : i.getParents()){
					
					if(j.getNumVisit() > 0) continue;
					
					if(sexLocked(j, exception))
						return true;
					
				}
				
								
			}
			
			return false;
		}

		
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
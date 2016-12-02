package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class SwapDescAnc extends Move {


	
	public SwapDescAnc(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose a sampled node
		Node child = currPedigree.getRandomSampledNode();
		
		
		//reject if the sampled node does not have any sampled ancestors
		List<Node> testing = new ArrayList<Node>();
		child.getAncestors(testing);
		List<Node> ancestors = new ArrayList<Node>();
		for(Node parent : testing){
			if(parent.sampled) ancestors.add(parent);
		}

		
		if(ancestors.size()==0)
			return REJECT;

		//choose ancestor
		Node anc = ancestors.get(currPedigree.rGen.nextInt(ancestors.size()));

		
		//reject if age incompatible 
		//TODO fix this!!!
		if(anc.getAge()!=-1 && child.getAge()!=-1)
			return REJECT;
		
		//reject if parent-offspring
		if(child.getParents().contains(anc))
			return REJECT;
		
		//check sex compatibility
		if(child.getSex() != anc.getSex()){
			
			for(Node c : child.getChildren()){
				if(c.getParents().size()==2) return REJECT;
			}
			
			for(Node c : anc.getChildren()){
				if(c.getParents().size()==2) return REJECT;
			}
			
		}
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
	

		
		//swap
		double prevLogLikelihood = currPedigree.getLogLikelihood();		
		currPedigree.swapAncDesc(anc, child);
		

		
		//likelihood
		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLogLikelihood, heat);
		
		
	}

	@Override
	protected void reverseMove(Pedigree currPedigree) {

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){

		return;
		
	}


	

	private boolean ageIncompatible(Pedigree currPedigree, Node child, Node parent){

		//both are ghost nodes
		if(child.getAge()==-1 && parent.getAge()==-1)
			return false;
		
		//both are individuals and have different ages
		if(child.getAge()!=-1 && parent.getAge()!=-1){
			return true;
		}
		
		
		if(child.getAge()!=-1){
			
			//child is too young to be an ancestor of its would-be descendants
			Node maxAgeDesc;
			for(Node c : parent.getChildren()){
				if(c==child) continue; //skip child
				if(c.sampled){
					maxAgeDesc = c;
				}
				else{
					maxAgeDesc = currPedigree.getDescendantWithMaxAge(c);
				}
				
				if(maxAgeDesc!=null && (child.getAge() <= maxAgeDesc.getAge())){
					return true;
				}
				
			}

				
		}
		
		
		if(parent.getAge()!=-1){
			
			//parent is too old to be its current ancestors's younger relative
			Node minAgeAnc = currPedigree.getAncestorWithMinAge(parent);
			if(minAgeAnc!=null && (minAgeAnc.getAge() <= parent.getAge())){
				return true;
			}
			
			
		}
		
		
		return false;
		
	}

	
}
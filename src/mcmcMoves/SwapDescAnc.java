package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) swap two sampled nodes

public class SwapDescAnc extends Move {

	List<Node> testing = new ArrayList<Node>();
	List<Node> ancestors = new ArrayList<Node>();

	
	public SwapDescAnc(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose a sampled node
		Node child = currPedigree.getRandomSampledNode();
		
		
		//reject if the sampled node does not have any sampled ancestors
		testing.clear();
		child.getAncestors(testing);
		ancestors.clear();
		
		int childDepth = child.getDepth();
		int childSex = child.getSex();
		for(Node parent : testing){
			//older than parents, sampled, and sex compatible
			if(parent.sampled && parent.getSex()==childSex && parent.getDepth() > childDepth+1)
				ancestors.add(parent);
		}

		
		if(ancestors.size()==0)
			return REJECT;

		//choose ancestor
		Node anc = ancestors.get(currPedigree.rGen.nextInt(ancestors.size()));


		
		/*
		//check sex compatibility
		if(child.getSex() != anc.getSex()){
			
			for(Node c : child.getChildren()){
				if(c.getParents().size()==2) return REJECT;
			}
			
			for(Node c : anc.getChildren()){
				if(c.getParents().size()==2) return REJECT;
			}
			
		}
		*/
		

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
	

		
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.swapAncDesc(anc, child);
		

		//this move is symmetric
		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, 0d, 0d, heat);
		
		
		
		
	}

	@Override
	protected void reverseMove(Pedigree currPedigree) {

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){

		return;
		
	}


	
	/*
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
	*/

	
}
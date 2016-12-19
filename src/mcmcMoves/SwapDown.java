package mcmcMoves;


import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class SwapDown extends Move {//works for 3 sampled nodes (2 parents, 1 child); need to test more thorougly
	
	
	public SwapDown(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {

		// randomly choose a sampled node
		Node parent = currPedigree.getRandomSampledNode();

		//reject if the sampled node is not connected to any other node
		if(parent.getNumEdges()==0)
			return REJECT;
		

		//bad cases
		if(parent.getDepth()==0) 
			return REJECT;
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		int sex = parent.getSex();
		
		//choose child
		double oldToNew = moveProbs.get("swapDown");
		Node child;
		List<Node> children = parent.getChildrenWithSex(sex);
		int nChildren = children.size();
		
		double makeNewProb = 1d/(1+nChildren);
		if(currPedigree.rGen.nextDouble() < makeNewProb){ //choose invisible child
			
			child = currPedigree.makeNewNode(parent.getDepth()-1, sex);
			child.addParent(parent);
			parent.addChild(child);
			
			oldToNew *= makeNewProb;
		}
		else{ //choose existing
			
			child = children.get(currPedigree.rGen.nextInt(children.size()));
			
			oldToNew *= (1 - makeNewProb) * 1d/nChildren;
			
		}

		
		
		
		//bad case
		/*
		if(ageIncompatible(currPedigree, child, parent)){
			
			if(nBefore < currPedigree.getNActiveNodes()){ //if a ghost node was created, destroy it
				currPedigree.deleteNode(child);

			}
			
			return REJECT;
		}
		*/

		
		//old to new
		if(child.sampled){ //different way of swapping if both child and parent are sampled
			oldToNew += moveProbs.get("swapUp");
		}
		oldToNew = Math.log(oldToNew);

		
		

		
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.swap(child, parent);
		


		//new to old
		double newToOld = moveProbs.get("swapUp");


		if(child.sampled){ //make child go down
			nChildren = child.getChildrenWithSex(sex).size();
			newToOld += moveProbs.get("swapDown") * (1-1d/(nChildren+1)) * 1d/nChildren;
		}
		
		newToOld = Math.log(newToOld);

		
	

		
		//likelihood
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
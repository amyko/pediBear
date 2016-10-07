package mcmcMoves;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class SwapUp extends Move {


	
	public SwapUp(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose an individual
		Node child = currPedigree.getRandomSampledNode();

		//reject if the sampled node is not connected to any other node
		if(child.getNumEdges()==0)
			return REJECT;
		

		//bad cases
		if(child.getDepth()==currPedigree.maxDepth || splitsPedigree(currPedigree, child)) 
			return REJECT;
		
		

		int nBefore = currPedigree.getNActiveNodes();
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		
		//choose parent
		Node parent = child.getParentWithSex(child.getSex());
		if(parent==null){
			parent = currPedigree.makeNewNode(child.getDepth()+1, child.getSex());
			parent.addChild(child);
			child.addParent(parent);
		}
		
		//bad case
		if(ageIncompatible(currPedigree, child, parent)){
			
			if(nBefore < currPedigree.getNActiveNodes()){ //if a ghost node was created, destroy it
				currPedigree.deleteNode(parent);
			}
			
			return REJECT;
		}

		
		//old to new
		double oldToNew = moveProbs.get("swapUp");
		if(parent.sampled){ //different way of swapping if both child and parent are sampled
			int nChildren = parent.getChildrenWithSex(parent.getSex()).size();
			oldToNew += (1 - 1d/(nChildren+1)) * 1d/nChildren * moveProbs.get("swapDown");
		}
		oldToNew = Math.log(oldToNew);

		


		
		//swap
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.swap(child, parent);
		


		//new to old
		double newToOld =  moveProbs.get("swapDown");
		
		//make child go down
		int nChildren = child.getChildrenWithSex(child.getSex()).size();

		//if parent is still there
		if(parent.getIndex() < currPedigree.getNActiveNodes()){
			newToOld *= (1 - 1d/(nChildren+1)) * 1d/nChildren;
		}
		else{
			newToOld *= 1d/(nChildren+1);	
		}

		if(parent.sampled){ //make parent go up
			newToOld += moveProbs.get("swapUp");
		}
		
		newToOld = Math.log(newToOld);

		
	

		
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

	
	private boolean splitsPedigree(Pedigree currPedigree, Node child){
		
		if(child==null) return false;
		if(child.getChildren().size()>0) return false;
		
		//parent to switch with
		Node parent = child.getParentWithSex(child.getSex());
		
		
		//has two parents and the parent to switch with is a ghost		
		if(child.getParents().size()==1 && parent==null) return true;

		if(child.getParents().size()==2 && !parent.sampled) return true;

		
		return false;
		
		
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
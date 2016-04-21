

package mcmcMoves;

import java.util.List;

import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class Swap extends Move {//works for 3 sampled nodes (2 parents, 1 child); need to test more thorougly

	private Node child;
	private Node parent;
	private boolean FStoPOswap; //true if child with a FS is swapped with a ghost parent AND the child does not have any descendants
	private int[] jDepthToCount = new int[maxDepth];
	private Node otherParent;
	
	
	public Swap(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//reset
		child = null;
		parent = null;
		FStoPOswap = false;
		
		// randomly choose an individual
		Node node = currPedigree.getRandomSampledNode();

		//reject if the sampled node is not connected to any other node
		if(node.getNumEdges()==0)
			return REJECT;
		
		
		// randomly choose up or down
		double oldToNew = 0d;
		if(node.getDepth()==currPedigree.maxDepth){//go down
			parent = node;
			oldToNew = 1d;
		}
		else if(node.getDepth()==0){//go up
			child = node;
			oldToNew= 1d;
		}
		else{//choose up or down randomly
			if(currPedigree.rGen.nextDouble() < .5){
				child = node;
			}
			else{
				parent = node;
			}
			oldToNew = .5;
		}
		

		
		
		//bad cases
		otherParent = null;
		if((child!=null && child.getDepth() == currPedigree.maxDepthForSamples) || splitsPedigree(currPedigree, child)) 
			return REJECT;
		
		
		//TODO reject fs to po swap
		//if(FStoPOswap)
			//return REJECT;


		int nBefore = currPedigree.getNActiveNodes();
		

		if(parent==null){//choose parent
			
			parent = child.getParentWithSex(child.getSex());
			if(parent==null){
				parent = currPedigree.makeNewNode(child.getDepth()+1, child.getSex());
				parent.addChild(child);
				child.addParent(parent);
			}

			//alternative way of swapping (i.e. swap parent with child)
			if(parent.sampled){
				double c = parent.getDepth()==currPedigree.maxDepth? 1 : .5;
				int nChildren = parent.getChildrenWithSex(parent.getSex()).size();
				oldToNew += c * (1-1d/(nChildren+1)) *1d/nChildren;
			}
			
			
		}
		else{//choose child
			List<Node> children = parent.getChildrenWithSex(parent.getSex());
			
			//switch with "invisible" child
			if(currPedigree.rGen.nextDouble() < 1d/(1+children.size())){
				child = currPedigree.makeNewNode(parent.getDepth()-1, parent.getSex());
				child.addParent(parent);
				parent.addChild(child);
				
				oldToNew *= 1d/(1+children.size());
			}
			
			//switch with existing child
			else{
				child = children.get(currPedigree.rGen.nextInt(children.size()));
				
				oldToNew *= (1 - 1d/(1+children.size())) * 1d/children.size();
			}
			
			//alternative way of swapping (i.e. swap child with parent)
			if(child.sampled){
				double c = child.getDepth()==0? 1 : .5;
				oldToNew += c;
			}
			
		}
		//oldToNew *= 1d/currPedigree.numIndiv * moveProbs[3];//these will cancel with denom
		oldToNew = Math.log(oldToNew);		
		
		//bad case
		if(ageIncompatible(currPedigree, child, parent)){
			
			if(nBefore < currPedigree.getNActiveNodes()){ //if a ghost node was created, destroy it
				if(!child.sampled) currPedigree.deleteNode(child);
				else currPedigree.deleteNode(parent);
			}
			
			return REJECT;
		}
		
		
		

		
		
		
		//swap
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.swap(child, parent);
		
		//determine if the parent is to be deleted
		boolean deleteParent = false;
		if(!parent.sampled){
			if(parent.getNumEdges()==2 && parent.getParents().size()==2)
				deleteParent = true;
			else if(parent.getNumEdges() < 2)
				deleteParent = true;
		}
		

		//new to old
		double newToOld = 0d;
		if(child.sampled){ //make child go down
			double c = child.getDepth()==currPedigree.maxDepth? 1 : .5;
			
			int nChildren = child.getChildrenWithSex(child.getSex()).size();

			//if parent is still there
			if(!deleteParent){
				newToOld = c * (1 - 1d/(nChildren+1)) * 1d/nChildren;
			}
			else{
				newToOld = c * 1d/nChildren;	
			}
			
		}
		
		if(parent.sampled){ //make parent go up
			double c = parent.getDepth()==0? 1 : .5;
			newToOld += c;
		}
		
		newToOld = Math.log(newToOld);
		

		//link probability for FStoPOswap
		if(FStoPOswap && parent.getChildren().size()==0 && !parent.sampled){
			

			Node iPrime = parent;
			double linkProb = 0d;
			
			//for every child of otherParent
			for(Node jPrime : otherParent.getChildren()){
				
				if(jPrime==iPrime) continue;
				
				currPedigree.getDepthToCount(jPrime, jDepthToCount);
				
				int targetDepth = otherParent.getDepth();
				int iDepth = iPrime.getDepth();
	
				for(int jDepth=0; jDepth <jPrime.getDepth()+1; jDepth++){
					if(iDepth==targetDepth && jDepth==targetDepth) continue;
					linkProb += jDepthToCount[jDepth] * getPowersOfHalf(3*targetDepth-Math.max(iDepth,jDepth)-iDepth-jDepth);
				}			
				
			}
			
			int nDeleted = otherParent.getNumEdges() > 2? 0 : 1;

			linkProb = Math.log(linkProb * moveProbs.get("link")) + getLogChooseTwo(currPedigree.getNActiveNodes() - nDeleted);

			newToOld += linkProb;
			
		}
		
	

		
		//likelihood
		double acceptRatio = heat * (currPedigree.getLogLikelihood() - prevLogLikelihood) + newToOld - oldToNew;

		
		if(acceptRatio > 0){
			return 1;
		}
		else{
			return Math.exp(acceptRatio);
		}
		
		
	}

	@Override
	protected void reverseMove(Pedigree currPedigree) {
		
		currPedigree.swap(parent, child);
		currPedigree.clean(child);
		currPedigree.clean(parent);
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		currPedigree.clean(child);
		currPedigree.clean(parent);
		
	}

	
	private boolean splitsPedigree(Pedigree currPedigree, Node child){
		
		if(child==null) return false;
		if(child.getChildren().size()>0) return false;
		
		//parent to switch with
		Node parent = child.getParentWithSex(child.getSex());
		
		//has two parents and the parent to switch with is a ghost
		if(child.getParents().size()==1 && parent==null) return true;
		
	
		if(child.getParents().size()==2 && !parent.sampled){
			
			
			//if child has a full sib, allow switch
			for(Node i : child.getParents().get(0).getChildren()){
				if(currPedigree.fullSibs(child, i)){
					
					FStoPOswap = true;				
					otherParent = child.getParents().get(0)==parent ? child.getParents().get(1) : child.getParents().get(0);
	
					return false;
				}
			}
			

			
			return true;

			
		}
		
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
				
				if(maxAgeDesc!=null && (child.getAge() - maxAgeDesc.getAge() < (child.getDepth() - maxAgeDesc.getDepth() + 1) * currPedigree.genTime)){
					return true;
				}
				
			}

				
		}
		
		
		if(parent.getAge()!=-1){
			
			//parent is too old to be its current ancestors's younger relative
			Node minAgeAnc = currPedigree.getAncestorWithMinAge(parent);
			if(minAgeAnc!=null && (minAgeAnc.getAge() - parent.getAge() < (minAgeAnc.getDepth() - parent.getDepth() + 1) * currPedigree.genTime)){
				return true;
			}
			
			
		}
		
		
		return false;
		
	}

	
}

package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;


public class Split extends Move {

	
	public Split(String name, double moveProb) {
		super(name, moveProb);
	}

	//this is stuff to store each move so that if it gets accepted it can then be performed
	private boolean hasFullSib;
	private int[] iDepthToCount = new int[maxDepth];
	private int[] jDepthToCount = new int[maxDepth];
	private List<Node> splitChildren = new ArrayList<Node>();
	private List<Node> stayChildren = new ArrayList<Node>();
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
				
		//get a random node to split
		Node parent = currPedigree.getRandomNode();
		

		//reject if node has < 2 children
		int nChildren = parent.getChildren().size();
		if(nChildren < 2)
			return REJECT;
		
	
		//randomly assign children to a clone (1 to 2^n-1); choose between [0,2^n - 2] and add 1
		int powerSetInd = currPedigree.rGen.nextInt((int) getPowersOfTwo(nChildren)-1) + 1;
		
		//no change if all children are assigned to clone
		if(!parent.sampled && powerSetInd==(int) getPowersOfTwo(nChildren)-1) 
			return REJECT;
		
		
		//assign children to split or stay
		splitChildren.clear();
		stayChildren.clear();
		List<Node> children = parent.getChildren();
	    for (int i = 0; i < nChildren; i++) {
	        if((1 << i & powerSetInd) != 0){ //if corresponding bit=1, add child
	        	splitChildren.add(children.get(i));
	        }  
	        else{
	        	stayChildren.add(children.get(i));
	        }
	    }
		

	    //check if any of the children between split and stay form full sibs
	    hasFullSib = false;
	    for(Node c1 : splitChildren){
	    	if(c1.getParents().size()!=2) continue;
	    	for(Node c2 : stayChildren){
	    		hasFullSib = currPedigree.fullSibs(c1, c2);
	    		if(hasFullSib) break;
	    	}
	    	if(hasFullSib) break;
	    }
	    
	    
	    //reject if creates illegal cycle
	    if(createsIllegalCycle(currPedigree, splitChildren, stayChildren, parent))
	    	return REJECT;
	    
	    
		//old to new
		//via split
	    int symmetric = (!parent.sampled && parent.getParents().size()==0) ? 1 : 0;
		double oldToNew = (1+symmetric) * getPowersOfHalf2(nChildren) * moveProbs.get("split");
	    
	    
		//save current pedigree
		currPedigree.copyCurrPedigree();
	    double prevLogLikelihood = currPedigree.getLogLikelihood();
	    int nBefore = currPedigree.getNActiveNodes();
	    int targetDepth = parent.getDepth();
	    
	    //split
	    Node splitParent = currPedigree.makeNewNode(parent.getDepth(), parent.getSex());
	    currPedigree.split(parent, splitParent, splitChildren, hasFullSib);
	    Node iPrime = currPedigree.clean(parent);
	    Node jPrime = currPedigree.clean(splitParent);	    
	    int nAfter = currPedigree.getNActiveNodes();


		//via cut
		oldToNew += .5 * (nBefore + 1 - nAfter) * moveProbs.get("cut");
		oldToNew = getLogChooseOne(nBefore) + Math.log(oldToNew);

		
		
		//new to old	
		double newToOld = 0d;
		double innerSum = 0d;
		double outerSum = 0d;
		
		currPedigree.getDepthToCount(iPrime, iDepthToCount);
		currPedigree.getDepthToCount(jPrime, jDepthToCount);
		//targetDepth = Math.max(targetDepth, iPrime.getDepth());

		
		for(int l1=0; l1 <= iPrime.getDepth(); l1++){
			innerSum = 0;
			for(int l2=0; l2 <= jPrime.getDepth(); l2++){
				
				//if(l1==targetDepth && l2==targetDepth) continue;
				
				innerSum += jDepthToCount[l2] * getPowersOfHalf(3*targetDepth - Math.max(l1,l2)-l1-l2); 
			}
			outerSum += iDepthToCount[l1] * innerSum;
		}
		newToOld = getLogChooseTwo(nAfter) + Math.log(moveProbs.get("link") * outerSum);
		

		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLogLikelihood, heat);
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {

		/*
		//node pointers
		for(int i=0; i<nActiveNodes; i++){
			Node model = nodes[i];
			Node toFix = currPedigree.getNode(i);
			
			toFix.setSex(model.getSex());
			toFix.setDepth(model.getDepth());

			toFix.getParents().clear();
			for(int idx : model.getParentIdx()){
				toFix.addParent(currPedigree.getNode(idx));
			}
			
			toFix.getChildren().clear();
			for(int idx : model.getChildrenIdx()){
				toFix.addChild(currPedigree.getNode(idx));
			}
			
			
		}
		
		//delete unecessary
		for(int i=nActiveNodes; i<currPedigree.getNActiveNodes(); i++){
			Node toDelete = currPedigree.getNode(i);
			toDelete.getParents().clear();
			toDelete.getChildren().clear();
			toDelete.setSex(-1);
			toDelete.setDepth(-1);
		}
		
		
		//relationships
		for(int i=0; i<currPedigree.numIndiv; i++){
			for(int j=0; j<currPedigree.numIndiv; j++){
				
				if(i==j) continue;
				
				Path modelRel = rel[i][j];
				
				Path toFix = currPedigree.getRelationships()[i][j];
				toFix.updatePath(modelRel.getUp(), modelRel.getDown(), modelRel.getNumVisit());

			}
		}
		
		
		currPedigree.setLogLikelihood(logLkhd);
		currPedigree.setNActiveNodes(nActiveNodes);
		*/
		
		currPedigree.reverse();

		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;

		
	}
	
	private boolean createsIllegalCycle(Pedigree currPedigree, List<Node> splitChildren, List<Node> stayChildren, Node parent){		

	    for(Node x : splitChildren){
	    	
	    	if(currPedigree.getFullSibs(x).size()==0) continue;
	    	
	    	for(Node y : splitChildren){
	    		
	    		if(x.getIndex()==y.getIndex()) continue;
	    		
	    		if(currPedigree.getFullSibs(y).size()==0) continue;
	    		
	    		if(!currPedigree.fullSibs(x,y)) return true;
	    		
	    	}
	    		
	    		
	    	
	    }
	    
	    return false;

	}
	
	

}

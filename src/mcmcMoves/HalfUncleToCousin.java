
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift down child cluster; 

public class HalfUncleToCousin extends Move{

	
	private List<Node> children  = new ArrayList<Node>();
	
	
	public HalfUncleToCousin(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		//reject if child doesn't have any siblings
		Node parent = child.getParents().get(0);
		if(parent.getChildren().size() < 2)
			return REJECT;
		
		
		//reject if moving down child two levels goes below depth=0
		currPedigree.clearVisit();
		for(Node k : child.getParents())
			k.setNumVisit(1);
		if(getMinDepth(child) - 1 < 0){
			return REJECT;
		}
		
		
		//get sibs
		children.clear();
		for(Node c : parent.getChildren()){
			if(c!=child) children.add(c);
		}
		
		Node sib = children.get(currPedigree.rGen.nextInt(children.size()));
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(children.size()) + Math.log(moveProbs.get("halfUncleToCousin"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.halfUncleToCousin(child, sib);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + lnTwo + Math.log(moveProbs.get("cousinToHalfUncle"));
		

		
		//accept ratio
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

		currPedigree.reverse();
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;
		
	}
	
	
	
	private int getMinDepth(Node node){

		int minDepth = node.getDepth();
		node.setNumVisit(1);
		
		for(Node c : node.getChildren()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMinDepth(c);
			if(currDepth < minDepth)
				minDepth = currDepth;
		}
		for(Node c : node.getParents()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMinDepth(c);
			if(currDepth < minDepth)
				minDepth = currDepth;
		}

		
		return minDepth;
		
		
	}
	
	

}

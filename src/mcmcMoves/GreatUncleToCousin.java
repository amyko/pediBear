
package mcmcMoves;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift down child cluster; 

public class GreatUncleToCousin extends Move{

	
	private Set<Node> children  = new HashSet<Node>();
	
	
	public GreatUncleToCousin(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child has no parent
		if(child.getParents().size()==0)
			return REJECT;
		//reject if not in cherry configuration
		children.clear();
		for(Node parent : child.getParents()){
			
			if(parent.sampled || parent.getParents().size()!=0 || parent.getChildren().size()!=2)
				return REJECT;
			
			children.addAll(parent.getChildren());
		}
		if(children.size()!=2) 
			return REJECT;
		
		//reject if moving down child two levels goes below depth=0
		currPedigree.clearVisit();
		for(Node k : child.getParents())
			k.setNumVisit(1);
		if(currPedigree.getMinDepth(child) - 2 < 0){
			return REJECT;
		}
		
		
		
		//get sib
		Node parent = child.getParents().get(0);
		Node sib = parent.getChildren().get(0) == child ? parent.getChildren().get(1) : parent.getChildren().get(0);
		List<Node> sibChildren = sib.getChildren();
		
		//reject if there is no sibChild 
		if(sibChildren.size()==0)
			return REJECT;
		
		
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(sibChildren.size()) + Math.log(moveProbs.get("greatUncleToCousin"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		Node sibChild = sibChildren.get(currPedigree.rGen.nextInt(sibChildren.size()));
		
		currPedigree.greatUncleToCousin(child, sib, sibChild);
		
		
		//new to old
		int nGrandParents = child.getParents().get(0).getParents().size();
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nGrandParents) + Math.log(moveProbs.get("cousinToGreatUncle"));
		

		
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
	
	
	

	
	
	
	
	

}



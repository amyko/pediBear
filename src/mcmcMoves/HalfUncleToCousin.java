
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

		
		//reject if child doesn't have any nephews
		Node parent = child.getParents().get(0);
		children.clear();
		for(Node c : parent.getChildren()){
			if(c!=child && c.getChildren().size() > 0)
				children.add(c);
		}
		
		if(children.size()==0)
			return REJECT;
		
		
		
		//reject if moving down child a level goes below depth=0
		currPedigree.clearVisit();
		for(Node k : child.getParents())
			k.setNumVisit(1);
		if(getMinDepth(child) - 1 < 0){
			return REJECT;
		}
		

		
		Node sib = children.get(currPedigree.rGen.nextInt(children.size()));
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		int nBefore = currPedigree.getNActiveNodes();

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.halfUncleToCousin(child, sib);
		
		//old to new
		//count number of equivalent sibs
		int nSibs = 0;
		for(Node i : currPedigree.getFullSibs(child.getParents().get(0))){
			if(i.getChildren().size()>0) nSibs++;
		}
		
		double oldToNew = getLogChooseOne(nBefore) + getLogChooseOne(children.size()) +  getLogChooseOne(2) + Math.log(nSibs * moveProbs.get("halfUncleToCousin"));
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("cousinToHalfUncle"));
		

		
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


package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift down child cluster; 

public class GreatUncleToCousin extends Move{

	
	private List<Node> goodFS  = new ArrayList<Node>();
	
	
	public GreatUncleToCousin(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have 2 parents
		if(child.getParents().size()!=2)
			return REJECT;

		//reject if no full sibs
		List<Node> fullSibs = currPedigree.getFullSibs(child);
		goodFS.clear();
		for(Node i : fullSibs){
			if(i.getChildren().size() > 0)
				goodFS.add(i);
		}
		if(goodFS.size()==0)
			return REJECT;
		

		//reject if moving down child two levels goes below depth=0
		currPedigree.clearVisit();
		for(Node k : child.getParents())
			k.setNumVisit(1);
		if(currPedigree.getMinDepth(child) - 2 < 0){
			return REJECT;
		}
		
		
		
		//get sib
		Node sib = goodFS.get(currPedigree.rGen.nextInt(goodFS.size()));
		List<Node> sibChildren = sib.getChildren();
		Node sibChild = sibChildren.get(currPedigree.rGen.nextInt(sibChildren.size()));

		//get number of equivalent choices of sib child
		int nSibChoice = 1;
		int nParents1 = sibChild.getParents().size();
		for(Node i : sibChildren){
			
			int nParents2 = i.getParents().size();
			
			if(i==sibChild || nParents2 != nParents1) continue;
			
			if(nParents1==1){//1 parent scenario
				if(sibChild.getParents().get(0) == i.getParents().get(0))
					nSibChoice++;
			}
			
			else{//two parent scenario
				if(sibChild.getParents().get(0)==i.getParents().get(0) && sibChild.getParents().get(1)==i.getParents().get(1))
					nSibChoice++;
				else if(sibChild.getParents().get(0)==i.getParents().get(1) && sibChild.getParents().get(1)==i.getParents().get(0))
					nSibChoice++;
			}
			
		}
		
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(goodFS.size()) + getLogChooseOne(sibChildren.size()) + getLogChooseOne(2) + Math.log(nSibChoice * moveProbs.get("greatUncleToCousin"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.greatUncleToCousin(child, sibChild);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("cousinToGreatUncle"));
		

		
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



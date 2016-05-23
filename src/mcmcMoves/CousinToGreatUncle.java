
package mcmcMoves;

import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent and link to two grand parents

public class CousinToGreatUncle extends Move{

	
	
	public CousinToGreatUncle(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child to cut
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		//reject if parent doesn't have exactly 1 child
		Node parent = child.getParents().get(0);
		if(parent.sampled || parent.getChildren().size()!=1 || parent.getParents().size()!=2)
			return REJECT;
		
		//reject if parent doesn't have full siblings
		int nSibChoice = currPedigree.getFullSibs(parent).size();
		if(nSibChoice==0)
			return REJECT;

		
		
		//depth constraint
		currPedigree.clearVisit();
		for(Node k : child.getParents()){
			k.setNumVisit(1);
		}
		int maxDepth = currPedigree.getMaxDepth(child);
		if(maxDepth+3 > currPedigree.maxDepth)
			return REJECT;

		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(2) + Math.log(moveProbs.get("cousinToGreatUncle"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		Node newSib = currPedigree.rGen.nextDouble() < .5 ? parent.getParents().get(0) : parent.getParents().get(1);
		currPedigree.cousinToGreatUncle(child, newSib);
		
		
		//new to old
		int nFS = 0;
		for(Node i : currPedigree.getFullSibs(child)){
			if(i.getChildren().size() > 0) nFS++;
		}
		int nFSChild = newSib.getChildren().size();
		

		
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nFS) + getLogChooseOne(nFSChild) + getLogChooseOne(2) + Math.log(nSibChoice* moveProbs.get("greatUncleToCousin"));
		

		
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


package mcmcMoves;

import java.util.List;

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
		int nParents = child.getParents().size();
		
		//depth constraint
		currPedigree.clearVisit();
		for(Node k : child.getParents()){
			k.setNumVisit(1);
		}
		int maxDepth = getMaxDepth(child);
		if(maxDepth+3 > currPedigree.maxDepth)
			return REJECT;
		//reject if child doesn't have exactly 1 parent
		if(nParents!=1)
			return REJECT;
		//reject if parent doesn't have exactly 1 child
		Node parent = child.getParents().get(0);
		if(parent.sampled || parent.getChildren().size()!=1 || parent.getParents().size()==0)
			return REJECT;

		
		//reject if gp has parents, or it would be deleted by removing parent
		for(Node gp : parent.getParents()){
			
			if(gp.getParents().size()!=0 || gp.getChildren().size() < 2)
				return REJECT;

		}
		
		//reject if split pedigree
		if(isSplitNode(parent)){
			return REJECT;
		}

		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//old to new
		List<Node> gp = parent.getParents();
		int nParent = gp.size();
		Node sib = gp.get(currPedigree.rGen.nextInt(nParent));
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(nParent) + Math.log(moveProbs.get("cousinToGreatUncle"));

		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();
		currPedigree.cousinToGreatUncle(child, sib);
		
		
		//new to old
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(sib.getChildren().size()) + Math.log(moveProbs.get("greatUncleToCousin"));
		

		
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

		
	private int getMaxDepth(Node node){

		int maxDepth = node.getDepth();
		node.setNumVisit(1);
		
		for(Node c : node.getChildren()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMaxDepth(c);
			if(currDepth > maxDepth)
				maxDepth = currDepth;
		}
		for(Node c : node.getParents()){
			if(c.getNumVisit() > 0) continue;
			int currDepth = getMaxDepth(c);
			if(currDepth > maxDepth)
				maxDepth = currDepth;
		}
		
		
		return maxDepth;
		
		
	}
	
	//returns true if removing this node splits the pedigree into two
	private boolean isSplitNode(Node node){
		
		if(node.sampled || node.getParents().size()==0 || node.getChildren().size() > 1) return false;
		if(node.getParents().size()==2) return true;
		
		//ghost and has 1 parent; recurse on parent
		return isSplitNode(node.getParents().get(0));
		
	}
		

	}

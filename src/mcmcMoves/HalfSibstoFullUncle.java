
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.Node;
import mcmc.MCMCMC;


//cut child from parent and link to two grand parents

public class HalfSibstoFullUncle extends Move{

	
	
	public HalfSibstoFullUncle(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a child
		Node child = currPedigree.getRandomNode();
		
		//reject if child doesn't have exactly 1 parent
		if(child.getParents().size()!=1)
			return REJECT;
		
		//reject if parent doesn't have other children
		Node parent = child.getParents().get(0);
		List<Node> halfSibs = new ArrayList<Node>();
		for(Node hs : parent.getChildren()){
			if(hs.getIndex() != child.getIndex())
				halfSibs.add(hs);
		}
		
		if(halfSibs.size()==0)
			return REJECT;
		
		
		//depth constraint
		currPedigree.clearVisit();
		parent.setNumVisit(1);
		int minDepth = currPedigree.getMinDepth(child);
		if(minDepth == 0) //can't shift down child cluster
			return REJECT;
		

		
		//choose half sib
		Node halfSib = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));
		
		/*
		//TODO test
		if(child.getIndex()==6 && halfSib.getIndex()==1){
			
			Path rel = currPedigree.getRelationships()[5][8];
			
			if(rel.getUp()==4 && rel.getDown()==3 && rel.getNumVisit()==1){
				
				rel = currPedigree.getRelationships()[0][1];
				
				if(rel.getUp()==1 && rel.getDown()==0 && rel.getNumVisit()==1){
					System.out.println();	
				}
				
			}

		}
		*/
		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		int nFS = currPedigree.getFullSibs(halfSib).size();
		
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(halfSibs.size()) +  Math.log((nFS+1) * moveProbs.get("halfSibsToFullUncle"));
		
		
		//cut and link
		double prevLkhd = currPedigree.getLogLikelihood();
		int targetSex = currPedigree.rGen.nextDouble() < .5? 0: 1;

		currPedigree.halfSibstoFullUncle(child, halfSib, targetSex);
		
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + Math.log(moveProbs.get("fullUncleToHalfSibs"));

		//accept ratio
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
	

}
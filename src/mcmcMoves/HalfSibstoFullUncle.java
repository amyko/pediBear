
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import dataStructures.Pedigree;
import dataStructures.Node;
import mcmc.SimulatedAnnealing;

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
		
		if(halfSibs.size()<2)
			return REJECT;
		
		
		//depth constraint
		currPedigree.clearVisit();
		for(Node k : child.getParents()){
			k.setNumVisit(1);
		}
		int minDepth = currPedigree.getMinDepth(child);
		if(minDepth == 0) //can't shift down child cluster
			return REJECT;

		
		//choose half sib
		Node halfSib = halfSibs.get(currPedigree.rGen.nextInt(halfSibs.size()));
		
		
					
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
	
		
		//cut and link
		double prevLogLikelihood = currPedigree.getLogLikelihood();

		currPedigree.halfSibstoFullUncle(child, halfSib);
		
		

		//accept ratio
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
	

}
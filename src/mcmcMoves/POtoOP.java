package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) swap two sampled nodes

public class POtoOP extends Move {

	List<Node> children = new ArrayList<Node>();

	
	public POtoOP(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose a sampled node
		Node upperNode = currPedigree.getRandomSampledNode();
		
		if(upperNode.getDepth() - 2 < 0)
			return REJECT;

		//overlap with shift cluster
		if(upperNode.getNumEdges()==0)
			return REJECT;
		
		int numChildren = upperNode.getChildren().size();

		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//choose middleNode
		double oldToNew = 0d;
		Node middleNode = null;
		if(currPedigree.rGen.nextDouble() < 1d/(1+numChildren)){
			
			int middleSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
			middleNode = currPedigree.makeNewNode(upperNode.getDepth()-1, middleSex);
			currPedigree.connect(upperNode, middleNode);
			
			oldToNew = 1d/(1+numChildren) * .5;
			
		}
		else{
			middleNode = upperNode.getChildren().get(currPedigree.rGen.nextInt(numChildren));
			
			oldToNew = numChildren/(1.0+numChildren) * (1d/numChildren);
		}
		

		
		Node lowerNode = null;

		oldToNew = Math.log(oldToNew * moveProbs.get("POtoOP"));
		
	
	
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.op2po(lowerNode, middleNode, upperNode);
	
		double newToOld = Math.log(moveProbs.get("OPtoPO"));

		//this move is symmetric
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
package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
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
		Node middleNode = null;
		if(currPedigree.rGen.nextDouble() < 1d/(1+numChildren)){
			
			int middleSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
			middleNode = currPedigree.makeNewNode(upperNode.getDepth()-1, middleSex);
			currPedigree.connect(upperNode, middleNode);

			
		}
		else{
			middleNode = upperNode.getChildren().get(currPedigree.rGen.nextInt(numChildren));

		}
		

		//reject if confounds with GP2HS
		if(!middleNode.sampled && middleNode.getChildren().size()>0 && upperNode.getNumEdges()==1 && middleNode.getParents().size()==1){

			
			Node child = middleNode.getChildren().get(0);
			
			if(currPedigree.getFullSibs(child).size()+1 == middleNode.getChildren().size())
				return REJECT;
			
			
		}

		
		
		
		Node lowerNode = null;

	
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.op2po(lowerNode, middleNode, upperNode);

		//this move is symmetric
		return SimulatedAnnealing.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, heat);

		
		
		
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
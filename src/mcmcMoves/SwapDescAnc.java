package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class SwapDescAnc extends Move {

	List<Node> testing = new ArrayList<Node>();
	List<Node> relatives = new ArrayList<Node>();
	
	public SwapDescAnc(String name, double moveProb) {
		super(name, moveProb);
	}

	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		
		// randomly choose a sampled node
		Node child = currPedigree.getRandomSampledNode();
		
		
		//reject if the sampled node does not have any sampled ancestors
		testing.clear();
		child.getAncestors(testing);
		child.getDescendants(testing);
		
		
		relatives.clear();
		
		for(Node parent : testing){
			//sex compatible
			if(parent.getSex()==child.getSex())
				relatives.add(parent);
		}

		
		if(relatives.size()==0)
			return REJECT;

		
		//choose ancestor
		Node secondNode = relatives.get(currPedigree.rGen.nextInt(relatives.size()));
		

		//reject if ancestor would be deleted
		if(!secondNode.sampled){
			if(child.getChildren().size()==0 || child.getNumEdges() < 2)
				return REJECT;
	
		}
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
	
		Node anc = secondNode;
		Node desc = child;
		if(secondNode.getDepth() < child.getDepth()){
			anc = child;
			desc = secondNode;
		}


		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.swapAncDesc(anc, desc);
		

		//likelihood
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
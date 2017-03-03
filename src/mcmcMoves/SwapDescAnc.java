package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


// 1) pick a random individual, 2) switch places with either its parent or child

public class SwapDescAnc extends Move {

	List<Node> testing = new ArrayList<Node>();
	List<Node> ancestors = new ArrayList<Node>();
	
	
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
		ancestors.clear();


		for(Node parent : testing){
			//sampled, and sex compatible
			if(parent.sampled && parent.getSex()==child.getSex())
				ancestors.add(parent);
		}

		
		if(ancestors.size()==0)
			return REJECT;

		//choose ancestor
		Node anc = ancestors.get(currPedigree.rGen.nextInt(ancestors.size()));

		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
	

		
		//swap
		double prevLkhd = currPedigree.getLogLikelihood();		
		currPedigree.swapAncDesc(anc, child);
		

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
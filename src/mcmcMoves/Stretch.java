package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;

//TODO this is not exactly the reverse move of contract; need to selectively shift down descendants of child

public class Stretch extends Move{ //WORKS; special merge not tested
	
	//this is stuff to store each move so that if it gets accepted it can then be performed
	private Node donor;
	private Node recipient;
	private List<Node> donorChildren = new ArrayList<Node>();
	private int nCuttableNode;
	int[] iDepthToCount = new int[maxDepth]; //keeps track of how many nodes there are in iCluster before any changes
	int[] jDepthToCount = new int[maxDepth]; //same for j
	private boolean specialMerge;
	private boolean mergingFormsFullSibs;
	
	private List<Node> halfSibs = new ArrayList<Node>();
	
	
	public Stretch(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomNode();
				
		//reject if child is a ghost AND doesn't have any descendants
		if(!child.sampled && child.getChildren().size()==0)
			return REJECT;
		
		//reject if child has no parents
		if(child.getParents().size()==0)
			return REJECT;
		
		//reject if shifting down the child cluster violates depth constraint
		currPedigree.clearVisit();
		int minCluDepth = currPedigree.getMinDepth(child);
		if(minCluDepth < 1) return REJECT;
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.stretch(child);

		
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

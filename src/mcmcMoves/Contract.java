package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Node;
import dataStructures.Pedigree;


public class Contract extends Move{ //WORKS; special merge not tested
	
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
	
	
	public Contract(String name, double moveProb) {
		super(name, moveProb);
	}


	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//choose child
		Node child = currPedigree.getRandomNode();
		
		//reject if it doesn't have exactly one parent
		if(child.getParents().size()!=1) return REJECT;
		Node parent = child.getParents().get(0);
		
		//reject if the parent is sampled or has wrong sex
		if(parent.sampled || parent.getSex()!=child.getSex())
			return REJECT;
		
		//reject depth violations
		currPedigree.clearVisit();
		int maxDepth = currPedigree.getMaxDepth(child);
		if(maxDepth == currPedigree.maxDepth)
			return REJECT;

		//reject bad cases
		if(violatesAgeConstraints(currPedigree, child, parent))
			return REJECT;


		//record previous config
		currPedigree.copyCurrPedigree();
		
		//merge
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.contract(parent, child);

		
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
	


	private boolean violatesAgeConstraints(Pedigree currPedigree, Node child, Node parent){

		if(parent.getChildren().size()==1 || child.getAge()==-1)
			return false;
		
		for(Node leaf : parent.getChildren()){
			
			if(leaf.getIndex()==child.getIndex()) continue;
			
			double maxClusterAge = currPedigree.getDescendantWithMaxAge(leaf).getAge();
			
			if(maxClusterAge > child.getAge())
				return true;
			
		}
		
		return false;
		
		
	}

	
	
	
}
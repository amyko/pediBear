package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//chooses a ghost node and switches sex of its parent and everyone in the chain, if possible

public class ShiftClusterLevel extends Move {

	
	private final List<Node> cluster = new ArrayList<Node>();
	private int offset;
	
	
	public ShiftClusterLevel(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	@Override
	protected double tryMove(Pedigree currPedigree, double heat) {
		
		//pick a random node
		Node node = currPedigree.getRandomNode();

		//get the cluster the node belongs to
		cluster.clear();
		currPedigree.clearVisit();
		node.getConnectedNodes(cluster);
		
		
		//get the lowest level of the cluster
		int oldLowestLevel = getLowestLevel(currPedigree);
		
		
		//pick the new lowest level
		int k = geometricDist(currPedigree.rGen);
		int newLowestLevel = k - 1;
		offset = newLowestLevel - oldLowestLevel;
		
		//reject if highest node goes over max depth or no change
		int highestLevel = getHighestLevel(currPedigree);
		if(highestLevel + offset > currPedigree.maxDepth)
			return REJECT;
		

		
		//shift cluster
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.shiftCluster(cluster, offset);
		
		
		//hastings ratio
		double oldToNew = Math.log(getPowersOfHalf(k));
		double newToOld = Math.log(getPowersOfHalf(oldLowestLevel+1));
		
		
		//accept ratio
		return MCMCMC.acceptanceRatio(currPedigree.getLogLikelihood(), prevLkhd, oldToNew, newToOld, heat);
		
	}
	
	
	@Override
	protected void reverseMove(Pedigree currPedigree) {		
		
		currPedigree.shiftCluster(cluster, -offset);
		
	}
	
	
	@Override
	protected void clean(Pedigree currPedigree){
		
		return;	
		
	}
	
	
	
	
	private int getLowestLevel(Pedigree currPedigree){
		
		int lowest = currPedigree.maxDepth;
		
		for(Node i : cluster){
			if(i.getDepth() < lowest){
				lowest = i.getDepth();
			}
			if(lowest==0) break;
		}
		
		
		return lowest;
		
		
	}

	
	private int getHighestLevel(Pedigree currPedigree){
		
		int highest = 0;
		
		for(Node i : cluster){
			if(i.getDepth() > highest){
				highest = i.getDepth();
			}
			if(highest==currPedigree.maxDepth) break;
		}
		
		
		return highest;
		
		
	}
	

	
}
	
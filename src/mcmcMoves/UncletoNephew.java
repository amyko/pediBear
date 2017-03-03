package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.SimulatedAnnealing;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class UncletoNephew extends Move{

	
	List<Node> fullSibs = new ArrayList<Node>();
	
	public UncletoNephew(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose uncle
		Node uncle = currPedigree.getRandomNode();
		
		//2 parents
		int nParents = uncle.getParents().size();
		if(nParents!=2) return REJECT;
		
		//full sibs
		List<Node> sibs = currPedigree.getFullSibs(uncle);
		
		//fulls sibs with children
		fullSibs.clear();
		for(Node x : sibs){
			if(x.getChildren().size()>0) fullSibs.add(x);
		}
		
		int nFS = fullSibs.size();
		if(nFS==0) return REJECT;
	
		//choose full sib
		Node fs = fullSibs.get(currPedigree.rGen.nextInt(nFS));
	
		//choose nephew
		List<Node> nephews = fs.getChildren();
		Node nephew = nephews.get(currPedigree.rGen.nextInt(nephews.size()));
		
		
		//choose shift direction
		double p = currPedigree.rGen.nextDouble();
		int uncleShift = -1;
		int nephewShift = 1;
		if(p < oneThird){
			uncleShift = -2;
			nephewShift = 0;
		}
		else if(p < twoThirds){
			uncleShift = 0;
			nephewShift = 2;
		}
		
		
		//depth constraint	
		if(uncleShift!=0){
			currPedigree.clearVisit();
			uncle.getParents().get(0).setNumVisit(1);
			uncle.getParents().get(1).setNumVisit(1);
			int minDepth = currPedigree.getMinDepth(uncle);
			if(minDepth + uncleShift < 0)
				return REJECT;
		}
		
		if(nephewShift!=0){
			currPedigree.clearVisit();
			uncle.setNumVisit(1);
			for(Node x : fs.getParents()){//don't count parents if it will be deleted
				if(!x.sampled && x.getNumEdges()==2) x.setNumVisit(1); 
			}
			int maxDepth = currPedigree.getMaxDepth(fs);
			if(maxDepth + nephewShift > currPedigree.maxDepth)
				return REJECT;
		}

		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
	
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		int targetSex = currPedigree.rGen.nextDouble() < .5 ? 0 : 1;
		currPedigree.uncle2nephew(uncle, nephew, fs, targetSex, uncleShift, nephewShift);

		
		//accept ratio
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
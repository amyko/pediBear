package mcmcMoves;

import java.util.ArrayList;
import java.util.List;

import mcmc.MCMCMC;
import dataStructures.Pedigree;
import dataStructures.Node;

//cut child from parent; shift parent cluster down; make full siblings

public class FStoSelf extends Move{

	
	List<Node> fullSibs = new ArrayList<Node>();
	
	public FStoSelf(String name, double moveProb) {
		super(name, moveProb);
	}
	
	
	protected double tryMove(Pedigree currPedigree, double heat) {
	

		//choose a node merge
		Node child = currPedigree.getRandomNode();

		//choose full sib
		fullSibs.clear();
		List<Node> temp = currPedigree.getFullSibs(child);
		for(Node x : temp){
			if(x.getSex()==child.getSex() && !(x.sampled && child.sampled))
				fullSibs.add(x);
		}
		
		
		//reject if child doesn't have any full sibs of same sex
		if(fullSibs.size()==0) return REJECT;
		

		//choose full sib
		Node fs = fullSibs.get(currPedigree.rGen.nextInt(fullSibs.size()));
		

		//donor
		Node donor = null;
		Node recipient = null;
		if(!fs.sampled){
			donor = fs;
			recipient = child;
		}
		else{
			donor = child;
			recipient = fs;
		}
		
		
		//oldToNew: choose node, choose fs to merge with; symmetry
		double oldToNew = getLogChooseOne(currPedigree.getNActiveNodes()) + getLogChooseOne(fullSibs.size()) + Math.log(2 * moveProbs.get("fs2self"));
		
		
		//copy pedigree
		currPedigree.copyCurrPedigree();
		
	
		//modify pedigree
		double prevLkhd = currPedigree.getLogLikelihood();
		currPedigree.FStoSelf(donor, recipient);
		
		
		//newToOld: choose node, choose children set, symmetery
		int symm = recipient.sampled ? 1 : 2;
		double newToOld = getLogChooseOne(currPedigree.getNActiveNodes()) + getPowersOfHalf2(recipient.getChildren().size()) + Math.log(symm * moveProbs.get("self2fs"));
		

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
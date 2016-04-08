package Unused;


import java.util.ArrayList;
import java.util.List;

import statistic.Statistic;
import dataStructures.Pedigree;
import mcmcMoves.Move;


public class MCMC {

	final Pedigree currPedigree;
	final int burnIn;
	final int runLength;
	final int sampleRate;
	public final Move[] moves; //TODO remove public
	final Statistic[] statsOfInterest;
	final List<Object[]> sampledPedigreeStats;
	
	public MCMC(Pedigree pedigree, int burnIn, int runLength, int sampleRate, Move[] moves, Statistic[] statsOfInterest){

		this.currPedigree = pedigree;
		this.burnIn = burnIn;
		this.runLength = runLength;
		this.sampleRate = sampleRate;
		this.moves = moves;
		this.statsOfInterest = statsOfInterest;
		this.sampledPedigreeStats = new ArrayList<Object[]>();
		
	}
	
	
	private Object[] sampleStats(Pedigree currPedigree){
		
		Object[] toReturn = new Object[statsOfInterest.length];
		
		for(int i=0; i<statsOfInterest.length; i++)
			toReturn[i] = statsOfInterest[i].getVal(currPedigree);
		
		return toReturn;
	}
	
	
	public List<Object[]> getSampledPedigreeStats(){
		return sampledPedigreeStats;
	}
	
	
	public void run(){

		//run burn-in
		for(int i = 0; i < burnIn; i++){
		
			double prevLkhd = currPedigree.getLogLikelihood();
			
			
			if(i>500 && i%moves.length==4){
				System.out.println("INVESTIGATE");
				currPedigree.printAdjMat();
			}
			
			
			
			
			
			//cycle through moves
			moves[i % moves.length].mcmcMove(currPedigree);

			
			System.out.println(i);
			System.out.println(i%moves.length);
			System.out.println(currPedigree.getNActiveNodes());
			currPedigree.printAdjMat();
			System.out.println(prevLkhd);
			System.out.println(currPedigree.getLogLikelihood());
			System.out.println();
			

		}
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			if(i % sampleRate == 0){
				sampledPedigreeStats.add(sampleStats(currPedigree));
			}
			moves[i % moves.length].mcmcMove(currPedigree);
		}
	}
	
	

	
	
	
}

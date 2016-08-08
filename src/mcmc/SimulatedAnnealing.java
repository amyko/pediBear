package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import mcmcMoves.Move;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import utility.DataParser;

public class SimulatedAnnealing {
	
	final Pedigree ped;
	final int runLength;
	final Move[] moves; 
	final PrintWriter writer;
	final PrintWriter convWriter;
	final Random rGen;
	
	//heating parameters
	private double[] heat;
	private int[] coolingSchedule;
	
	
	
	public double bestLkhd = Double.NEGATIVE_INFINITY;
	

	public SimulatedAnnealing(Pedigree ped, double[] heat, int[] coolingSchedule, Move[] moves, int runLength, Random rGen, String outPath, PrintWriter convWriter) throws IOException{

		this.ped = ped;
		this.runLength = runLength;
		this.moves = moves;		
		this.writer = DataParser.openWriter(outPath);
		this.convWriter = convWriter;
		
		this.rGen = rGen;
		this.heat = heat;
		this.coolingSchedule = coolingSchedule;
		
		
	}
	
	
	
	public void run(){

		runBurnIn();
		
		runSample();
	
	}

	
	
	
	private void runBurnIn(){
		
		System.out.println("Burn in...");
		
		convWriter.write(">\n");

		
		for(int t=0; t<heat.length-1; t++){

			//run burn-in
			//while(moves[0].nTried < 100 && moves[0].nAccept <50){
			for(int i = 0; i < coolingSchedule[t]; i++){
			

				Move move = chooseMove();
					
				/*
				//TESTING			
				if(!ped.sanityCheck() || i==9724){
					System.out.println(String.format("(%s,%d)", move.name, i));
				
					for(int k=0; k< ped.getNActiveNodes(); k++){
						ped.getNode(k).print();
					}
					//System.out.println();
						
					System.out.println(ped.getNActiveNodes());
					ped.printAdjMat();
					System.out.println();
						
				}
				*/

				move.mcmcMove(ped, heat[t]);
				
					
			
			}
			
			
			//count num success
			//System.out.print(String.format("%d\t%d\t%d\t%f\n", t, moves[0].nAccept, moves[0].nTried, ped.getLogLikelihood()));
			moves[0].nAccept = 0;
			moves[0].nTried = 0;
			
			//record likelihood
			convWriter.write(String.format("%d\t%f\n", coolingSchedule[t]*(t+1), ped.getLogLikelihood()));
			
			
		}
		
		

		
		
	}
	
	
	
	
	
	
	private void runSample(){
		
		System.out.println("Sampling...");
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			//record best likelihood
			double currLkhd = ped.getLogLikelihood();
			if(currLkhd > this.bestLkhd){
				this.bestLkhd = currLkhd;
				//System.out.println(bestLkhd);
			}
			
			

			sample(ped);
			//sample2(ped);

			
			//update	
			Move move = chooseMove();
			move.mcmcMove(ped, heat[heat.length-1]);

			

	
		}
		
		
		//close outfile
		writer.close();
		
	}
	
	

	
	private Move chooseMove(){
		
		double u = rGen.nextDouble();
		double cumProb = 0;
		
		for(Move move : moves){
			
			cumProb += move.getProb();
			
			if(u < cumProb)
				return move;
			
		}
		
		throw new RuntimeException("No move chosen");
		
		
	}
	
	

	//write relationship to file
	private void sample(Pedigree currPedigree){
		
		//header for this sample
		writer.write(String.format(">\t%.5f\n", currPedigree.getLogLikelihood()));
		
		for(int i=0; i<currPedigree.numIndiv; i++){
			for(int j=i+1; j<currPedigree.numIndiv; j++){
				
				Path rel = currPedigree.getRelationships()[i][j];
				
				writer.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
			}
		}
		
	}
	
	
	
	
	public static double acceptanceRatio(double proposalLkhd, double currLkhd, double heat){
		
		//accept if uphill
		if(proposalLkhd > currLkhd){
			return 1d;
		}
		
		//reject if frozen
		if(heat==0) return 0d;
		
		//accept with prob
		double prob = Math.exp(heat * (proposalLkhd - currLkhd));
		
		return prob;
	}
	
	

	
	
}

package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import mcmcMoves.Move;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import utility.ArrayUtility;
import utility.DataParser;

public class MCMCMC {
	
	//TODO these should come in as input
	final static double swapRateLow = .2;
	final static double swapRateHigh = .4;
	final static int tuneIter = 10000;
	final static int maxTuneTrialTune = 10; //max number of tune trials inside Tune
	final static int maxTuneTrialBurnIn = 10; //max number of tune trials inside BurnIn
	final static int tuneInterval = 50000;
	
	final List<Pedigree> chains;
	final int burnIn;
	final int runLength;
	final int sampleRate;
	final Move[] moves; 
	final PrintWriter writer;
	final PrintWriter writer2;
	final Random rGen;
	final int nChain;
	public int nSwapSuccess;
	public int nSwapAttempt;
	public int coldChain;
	public int swapInterval;
	boolean tuned;
	
	private double deltaT;
	private double[] heat;
	
	
	
	public double bestLkhd = Double.NEGATIVE_INFINITY;
	

	//TODO parallelize 
	public MCMCMC(List<Pedigree> chains, double deltaT, Move[] moves, int burnIn, int runLength, int sampleRate, int swapInterval, Random rGen, String outPath) throws IOException{

		this.chains = chains;
		this.deltaT = deltaT;
		this.burnIn = burnIn;
		this.runLength = runLength;
		this.sampleRate = sampleRate;
		this.moves = moves;		
		this.writer = DataParser.openWriter(outPath);
		this.writer2 = DataParser.openWriter(outPath+".2");
		this.rGen = rGen;
		this.nChain = chains.size();
		this.nSwapAttempt = 0;
		this.nSwapSuccess = 0;
		this.coldChain = 0;
		this.swapInterval = swapInterval;
		this.tuned = false;
		
		this.heat = new double[nChain];
		for(int i=0; i<nChain; i++) 
			heat[i] = 0;
			//heat[i] = 1 / (1 + deltaT*i);
		
		
	}
	
	
	
	public void run(){

		tune();

		runBurnIn();
		
		runSample();
	
	}
	
	
	
	//tune delta to achieve the desired swap rate
	private void tune(){
		
		if(nChain < 2) return;
		
		System.out.println("Tuning...");
		
		//initialize variables
		int t = 0;
		nSwapAttempt = 0;
		nSwapSuccess = 0;
		double currSwapRate = 0;
		
		//tune 
		while(!tuned && t < maxTuneTrialTune){
		
			//run mcmc
			for(int i=0; i<tuneIter; i++){
				
				//for every chain, update
				for(int j = 0; j < nChain; j++){
					
					Move move = chooseMove();				
					move.mcmcMove(chains.get(j), heat[j]);
						
					
					
									
				}
				
				if(i%swapInterval==0){
					swapStates();
				}
				
			}
			
			
			//check swap rate
			currSwapRate = (double) nSwapSuccess/ nSwapAttempt / swapInterval;
			
			if(currSwapRate > swapRateLow && currSwapRate < swapRateHigh){
				tuned = true;
			}
			
			else{
				double multiplier = currSwapRate < swapRateLow? .5 : 2;
				this.deltaT *= multiplier;
				for(int i=0; i<nChain; i++) heat[i] = 1 / (1 + deltaT*i);
			}

			//update counts
			nSwapAttempt = 0;
			nSwapSuccess = 0;
			t++;
			
		}
		
		//System.out.println(String.format("swap rate: %.2f", currSwapRate));
		//if(t==maxTuneTrialTune && !tuned)
			//System.out.println("Tuning failed inside tune()");
		
		
		
		
	}
	
	
	
	
	private void runBurnIn(){
		
		System.out.println("Burn in...");
		
		//init variables
		int tuneNum = 0;
		
		//run burn-in
		for(int i = 0; i < burnIn; i++){
		
			//for every chain, update
			for(int j = 0; j < nChain; j++){
				

				Move move = chooseMove();
				
				/*
				//TESTING			
				if(!chains.get(j).sanityCheck() || true){
					System.out.println(String.format("(%s,%d,%d)", move.name, i, j));
				
					for(int k=0; k< chains.get(j).getNActiveNodes(); k++){
						chains.get(j).getNode(k).print();
					}
					//System.out.println();
						
					System.out.println(chains.get(j).getNActiveNodes());
					chains.get(j).printAdjMat();
					System.out.println();
					

						
				}
				*/
				
				
				
				move.mcmcMove(chains.get(j), heat[j]);
				
				
		
			}
			

			
			//swap
			if(i%swapInterval==0){
				swapStates();
			}
			
			
			//retune, if necessary
			if(i%tuneInterval==0 && tuneNum < maxTuneTrialBurnIn){
				
				double currSwapRate = (double) nSwapSuccess / nSwapAttempt / swapInterval;

				if(currSwapRate < swapRateLow || currSwapRate > swapRateHigh){
					tuned = false;
					tune();
					tuneNum++;
					//reset burnIn counter
					i = 0;
				}
				
				if(!tuned && tuneNum==maxTuneTrialBurnIn){
					System.out.println("WARNING: Auto-tune failed");
					nSwapAttempt = 0;
					nSwapSuccess = 0;
				}
				

				
			}
			
			


		}
		
		
	}
	
	
	
	
	
	
	private void runSample(){
		
		System.out.println("Sampling...");
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			//record best likelihood
			double currLkhd = chains.get(this.coldChain).getLogLikelihood();
			if(currLkhd > this.bestLkhd){
				this.bestLkhd = currLkhd;
				//System.out.println(bestLkhd);
			}
			
			
			//sample from cold chain
			if(i % sampleRate == 0){
				sample(chains.get(this.coldChain));
				sample2(chains.get(this.coldChain));
			}
			
			
			//for every chain, update
			for(int j = 0; j < nChain; j++){				
				Move move = chooseMove();
				move.mcmcMove(chains.get(j), heat[j]);
			}
			
			if(i%swapInterval==0){
				swapStates();
			}
			
	
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
	

	
	
	
	private void swapStates(){
		
		if(nChain<2) return;
		
		nSwapAttempt++;
		
		//randomly choose two chains to swap
		int[] twoChains = ArrayUtility.getNRandomIdx(nChain, 2, rGen);
		int j = twoChains[0];
		int k = twoChains[1];
		
		
		//compute probability of swapping states
		double acceptRatio = heat[j] * chains.get(k).getLogLikelihood() + heat[k] * chains.get(j).getLogLikelihood() - heat[j] * chains.get(j).getLogLikelihood() - heat[k] *chains.get(k).getLogLikelihood();
		double acceptProb = 0d;
		if(acceptRatio > 0){
			acceptProb = 1;
		}
		else{
			acceptProb = Math.exp(acceptRatio);
		}
		
		
		//swap states, which is equivalent to swapping heat parameters
		if(rGen.nextDouble() < acceptProb){
			
			double temp = heat[j];
			heat[j] = heat[k];
			heat[k] = temp;
			
			nSwapSuccess++;
			
		}
		
		
		//update index of cold chain
		if(heat[j]==1){
			this.coldChain = j;
		}
		else if(heat[k]==1){
			this.coldChain = k;
		}
		
		
		
		
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
	
	
	
	private void sample2(Pedigree currPedigree){
		
		//header for this sample
		writer2.write(String.format(">\t%.5f\t%d\n", currPedigree.getLogLikelihood(), currPedigree.getNActiveNodes()));
		
		for(int i=0; i<currPedigree.getNActiveNodes(); i++){

			Node node = currPedigree.getNode(i);
			List<Node> parents = node.getParents();
			int p1 = -1;
			int p2 = -1;
			if(parents.size()==1){
				p1 = parents.get(0).getIndex();
			}
			else if(parents.size()==2){
				p1 = parents.get(0).getIndex();
				p2 = parents.get(1).getIndex();
			}
			
			int sampled = node.sampled ? 1 : 0;
			
			writer2.write(String.format("%d\t%d\t%d\t%d\t%d\n", node.getIndex(), p1, p2, node.getDepth(), sampled));

		}
		
		
	}
	
	

	
	
}

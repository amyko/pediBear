package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.List;
import java.util.Random;

import mcmcMoves.Move;
import dataStructures.Chain;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import utility.ArrayUtility;
import utility.DataParser;

public class MCMCMC {
	
	//TODO these should come in as input
	final static double swapRateLow = .2;
	final static double swapRateHigh = .6;
	final static int tuneIter = 10000;
	final static int maxTuneTrialTune = 10; //max number of tune trials inside Tune
	final static int maxTuneTrialBurnIn = 10; //max number of tune trials inside BurnIn
	final static int tuneInterval = 100000;

	
	final List<Chain> chains;
	final int burnIn;
	final int runLength;
	final int sampleRate;
	final Move[] moves; 
	final PrintWriter writer;
	final PrintWriter convWriter;
	final PrintWriter cranefootFamWriter;
	final Random rGen;
	final int nChain;
	public int nSwapSuccess;
	public int nSwapAttempt;
	public int coldChain;
	public int swapInterval;
	boolean tuned;
	int nSwaps;
	
	private double deltaT;
	
	private int missingParentCounter = 0;

	
	
	
	public double bestLkhd = Double.NEGATIVE_INFINITY;
	

	//TODO parallelize 
	public MCMCMC(List<Chain> chains, double deltaT, Move[] moves, int burnIn, int runLength, int sampleRate, int swapInterval, int nSwaps, Random rGen, String outPath) throws IOException{

		this.chains = chains;
		this.deltaT = deltaT;
		this.burnIn = burnIn;
		this.runLength = runLength;
		this.sampleRate = sampleRate;
		this.moves = moves;		
		this.writer = DataParser.openWriter(outPath+".pair");
		this.convWriter = DataParser.openWriter(outPath+".lkhd");
		this.cranefootFamWriter = DataParser.openWriter(outPath+".fam");
		this.rGen = rGen;
		this.nChain = chains.size();
		this.nSwapAttempt = 0;
		this.nSwapSuccess = 0;
		this.coldChain = nChain - 1;
		this.swapInterval = swapInterval;
		this.tuned = false;
		this.nSwaps = nSwaps;
		


		
	}
	
	
	
	public void run(){

		runBurnIn();
		
		runSample();
		
		//record last fam
		writeFamFile();
	
	}
	
	
	
	//tune delta to achieve the desired swap rate
	private void tune(){
		
		//if(true) return;
		
		//System.out.println("Tuning...");
		
		//initialize variables
		int t = 0;

		
		//tune 
		while(!tuned && t < maxTuneTrialTune){
		

			//check swap rate
			double currSwapRate = (double) nSwapSuccess/ nSwapAttempt / swapInterval;
			
			if(currSwapRate > swapRateLow && currSwapRate < swapRateHigh){
				tuned = true;
			}
			
			else{
				double multiplier = currSwapRate < swapRateLow? .8 : 1.2; 
				
				//update delta
				this.deltaT *= multiplier;
				
				for(int i=0; i<nChain; i++){
					chains.get(i).setHeat(deltaT);
					
				}
				
				
			}

			//update counts
			nSwapAttempt = 0;
			nSwapSuccess = 0;
			t++;
			
			
			//run mcmc
			for(int i=0; i<tuneIter; i++){
				
				//for every chain, update
				for(int j = 0; j < nChain; j++){
					
					Move move = chooseMove();				
					move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat());
						
					
					
									
				}
				
				if(i%swapInterval==0){
					swapStates();
				}
				
			}
			
			
			
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
			
			//System.out.println(i);
		
			//for every chain, update
			for(int j = 0; j < nChain; j++){
				

				Move move = chooseMove();
				
				/*
				//TESTING			
				if(!chains.get(j).sanityCheck()){
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


				
				move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat());
				
				
		
			}
			

			
			//swap
			if(i%swapInterval==0){
				swapStates();
			}
			
			
			
			//retune, if necessary
			if((i+1)%tuneInterval==0 && tuneNum < maxTuneTrialBurnIn){
				

				double currSwapRate = (double) nSwapSuccess / nSwapAttempt / swapInterval;
				
				//TODO testing
				//get minimum heat
				double minHeat = 1;
				for(Chain x : chains){
					if(x.getHeat() < minHeat) minHeat = x.getHeat();
				}
						
						
				System.out.println(String.format("%f, %f", currSwapRate, minHeat));
				
				
				if(currSwapRate < swapRateLow || currSwapRate > swapRateHigh){
					tuned = false;
					tune();
					tuneNum++;

				}
				
				if(!tuned && tuneNum==maxTuneTrialBurnIn){
					System.out.println("WARNING: Auto-tune failed");
				}
				
				
				

				
			}
			
			
			
			
			


		}

		
		
	}
	
	
	
	
	
	
	private void runSample(){
		
		System.out.println("Sampling...");
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			//record best likelihood
			double currLkhd = chains.get(this.coldChain).getPedigree().getLogLikelihood();
			if(currLkhd > this.bestLkhd){
				this.bestLkhd = currLkhd;
				//System.out.println(bestLkhd);
			}
			
			
			//sample from cold chain
			if(i % sampleRate == 0){
				sample(chains.get(this.coldChain).getPedigree());
				convergence(chains.get(this.coldChain).getPedigree());
			}
			
			
			//for every chain, update
			for(int j = 0; j < nChain; j++){				
				Move move = chooseMove();
				move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat());
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
		
		for(int q=0; q < nSwaps; q++){
			
			

			int j = rGen.nextInt(nChain-1);
			int k = chains.get(j).getRandomNeighborIndex();
			
			
			//compute probability of swapping states
			double acceptRatio = chains.get(j).getHeat() * chains.get(k).getLikelihood() + chains.get(k).getHeat() * chains.get(j).getLikelihood() - chains.get(j).getHeat()  * chains.get(j).getLikelihood() - chains.get(k).getHeat()  * chains.get(k).getLikelihood();
			double acceptProb = 0d;
			if(acceptRatio > 0){
				acceptProb = 1;
			}
			else{
				acceptProb = Math.exp(acceptRatio);
			}
			
		
			//swap states
			if(rGen.nextDouble() < acceptProb){
				
				Pedigree jped = chains.get(j).getPedigree();
				chains.get(j).setPedigree(chains.get(k).getPedigree());
				chains.get(k).setPedigree(jped);

				
				nSwapSuccess++;
				
			}
			
			nSwapAttempt++;
			

		}
		
		
		
		
	}
	

	//write relationship to file
	private void sample(Pedigree currPedigree){
		
		//pairwise relatioship
		//header for this sample
		writer.write(String.format(">\t%f\n", currPedigree.getLogLikelihood()));

		
		for(int i=0; i<currPedigree.numIndiv; i++){
			
			for(int j=i+1; j<currPedigree.numIndiv; j++){
				
				Path rel = currPedigree.getRelationships()[i][j];
				
				
				//writer.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
				//TODO for hastings test
				writer.write(String.format("%d\t%d\t%d\t", rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
			}
			
			
		}
		
		//TODO testing
		writer.write("\n");
	

	}
	
	private void writeFamFile(){
		
		Pedigree currPedigree = chains.get(coldChain).getPedigree();
		
		//write family relationship
		cranefootFamWriter.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		for(int i=0; i<currPedigree.getNActiveNodes(); i++){
			recordCranefootFam(currPedigree.getNode(i), currPedigree);
		}
		
		cranefootFamWriter.close();
		
	}
	
	
	private void convergence(Pedigree currPedigree){
		
		convWriter.write(String.format("%f\n", currPedigree.getLogLikelihood()));
		
	}
	

	private void recordCranefootFam(Node ind, Pedigree currPedigree){
		
		String name = ind.fid + "_" + ind.iid;
		String pa = "0";
		String ma = "0";
		String sampleStatus = ind.sampled ? "000000" : "999999";
		String sex = ind.getSex()==1 ? "1" : "7"; 
		
		//if missing individual and sex not constrained
		currPedigree.clearVisit();
		if(currPedigree.sexLocked(ind)==false) sex = "4";
			
		//get parent ids
		for(Node parent : ind.getParents()){
			
			//recordFam(parent);
		
			if(parent.getSex()==0)
				ma = parent.fid + "_" + parent.iid;
			else if(parent.getSex()==1)
				pa = parent.fid + "_" + parent.iid;
			else
				throw new RuntimeException("Parent with unknown sex");
			
		}
		
		//if only one parent is present
		if(ind.getParents().size()==1){
			
			//make missing parent
			int missingParentSex = ind.getParents().get(0).getSex()==1 ? 0 : 1;
			Node missingParent = new Node("missingParent", missingParentCounter+"", missingParentSex, false, -1);
			
			//connect temporarily
			currPedigree.connect(missingParent, ind);
			recordCranefootFam(missingParent, currPedigree);
			currPedigree.disconnect(missingParent, ind);
			
			if(missingParentSex==0) ma = String.format("missingParent_%d", missingParentCounter);
			else pa = String.format("missingParent_%d", missingParentCounter);
			
			missingParentCounter++;
			
		}

		
		
		
		
		//write to file
		cranefootFamWriter.write(String.format("%s\t%s\t%s\t%s\t%s\n", name, pa, ma, sex, sampleStatus));
		
		
	}
	
	
	
	public static double acceptanceRatio(double newLkhd, double oldLkhd, double oldToNew, double newToOld, double heat){
		
		
		if(newLkhd > 0){
			System.out.println("Positive lkhd");
		}
		
		
		
		if(oldToNew==Double.NEGATIVE_INFINITY || newToOld==Double.NEGATIVE_INFINITY){
			System.out.println("proposal prob inifinity");
		}
		
		
		//TODO testing
		//heat = 0;
		
		
		
		double acceptRatio = heat * (newLkhd - oldLkhd) + newToOld - oldToNew;
		

		
		if(acceptRatio > 0){
			return 1;
		}
		else{
			return Math.exp(acceptRatio);
		}

	}
	
	

	
	
}

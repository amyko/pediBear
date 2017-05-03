package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import mcmcMoves.Move;
import dataStructures.Chain;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.PedInfo;
import dataStructures.Pedigree;
import utility.DataParser;

public class MCMCMC {
	
	//TODO these should come in as input
	final static double swapRateLow = .2;
	final static double swapRateHigh = .6;
	final static int tuneIter = 10000;
	final static int maxTuneTrialTune = 10; //max number of tune trials inside Tune
	final static int maxTuneTrialBurnIn = 10; //max number of tune trials inside BurnIn
	final static int tuneInterval = 100000;
	final static int maxDelta = 10;
	final static int minDelta = 0;

	String outPath; 
	final List<Chain> chains;
	final int burnIn;
	final int runLength;
	final int sampleRate;
	final Move[] moves; 
	final PrintWriter writer;
	final PrintWriter convWriter;
	PrintWriter cranefootFamWriter;
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
	public Map<String, PedInfo> ped2info = new HashMap<String, PedInfo>();


	
	
	
	public double bestLkhd = Double.NEGATIVE_INFINITY;
	

	//TODO parallelize 
	public MCMCMC(List<Chain> chains, double deltaT, Move[] moves, int burnIn, int runLength, int sampleRate, int swapInterval, int nSwaps, Random rGen, String outPath) throws IOException{

		this.outPath = outPath;
		this.chains = chains;
		this.deltaT = deltaT;
		this.burnIn = burnIn;
		this.runLength = runLength;
		this.sampleRate = sampleRate;
		this.moves = moves;		
		this.writer = DataParser.openWriter(outPath+".pair");
		this.convWriter = DataParser.openWriter(outPath+".lkhd");
		//this.cranefootFamWriter = DataParser.openWriter(outPath+".fam");
		this.rGen = rGen;
		this.nChain = chains.size();
		this.nSwapAttempt = 0;
		this.nSwapSuccess = 0;
		this.coldChain = nChain - 1;
		this.swapInterval = swapInterval;
		this.tuned = false;
		this.nSwaps = nSwaps;
		


		
	}
	
	
	
	public void run() throws IOException{

		runBurnIn();
		
		runSample();

	
	}
	
	
	
	//tune delta to achieve the desired swap rate
	private void tune(){
		
		//if(true) return;
		
		//System.out.println("Tuning...");
		
		//initialize variables
		int t = 0;

		boolean deltaOutOfBounds = false;
		
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
				
				//failed auto tune
				if(deltaT > maxDelta || deltaT < minDelta){
					deltaOutOfBounds = true;
					break;
				}
				
				for(int i=0; i<nChain; i++){
					chains.get(i).setHeat(deltaT);
					
				}
				
				
			}
			
			//failed to auto-tune
			if(deltaOutOfBounds) break;

			//update counts
			nSwapAttempt = 0;
			nSwapSuccess = 0;
			t++;
			
			
			//run mcmc
			for(int i=0; i<tuneIter; i++){
				
				//for every chain, update
				for(int j = 0; j < nChain; j++){
					
					Move move = chooseMove();				
					move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat(),0);
						
					
					
									
				}
				
				if(i%swapInterval==0){
					swapStates();
				}
				
			}
			
			
			
		}
		

		//System.out.println(String.format("swap rate: %.2f", currSwapRate));
		if(t==maxTuneTrialTune && !tuned)
			System.out.println("Tuning failed");
		if(deltaOutOfBounds){
			System.out.println("Tuning paramater delta out of bounds");
		}
		
		
		
		
	}
	
	
	
	
	private void runBurnIn() throws IOException{
		
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
				if(!chains.get(j).getPedigree().sanityCheck()){
					System.out.println(String.format("(%s,%d,%d)", move.name, i, j));

						
				}
				
				 
				
				if(i==1838){
					//Node myNode = chains.get(j).getPedigree().getNode(6);
					//String toWrite = "";
					//for(Node x : myNode.getParents()) toWrite += x.iid+" ";
					//System.out.println(toWrite);
					//cranefootFamWriter = DataParser.openWriter(outPath+".fam");
					//writeFamFile(chains.get(j).getPedigree());
					//System.out.println(String.format("%d %d", i, j));
					System.out.println(move.name);
					//System.out.println(chains.get(j).getLikelihood());
					
					chains.get(j).getPedigree().printAdjMat();
					System.out.println();
					
					
				}
				*/
				
				
				
				
				//chains.get(j).getPedigree().printAdjMat();
				//System.out.println();
				

				move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat(), i);
				
				
		
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
	
	
	
	
	
	
	private void runSample() throws IOException{
		
		System.out.println("Sampling...");
		
		boolean written = false;
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			//record best likelihood
			double currLkhd = chains.get(this.coldChain).getPedigree().getLogLikelihood();
			if(currLkhd > this.bestLkhd){
				this.bestLkhd = currLkhd;
				
				cranefootFamWriter = DataParser.openWriter(outPath+".fam");
				missingParentCounter = 0;
				writeFamFile(chains.get(coldChain).getPedigree());
				
			}
			
			
			if(Math.abs(currLkhd+47393.39) < 1e-1 && written==false){

				cranefootFamWriter = DataParser.openWriter(outPath+".2.fam");
				missingParentCounter = 0;
				writeFamFile(chains.get(coldChain).getPedigree());	
				written = true;
			}
			
			
			
			
			
	
			
			
			//sample from cold chain
			if(i % sampleRate == 0){
				sample(chains.get(this.coldChain).getPedigree());
				convergence(chains.get(this.coldChain).getPedigree());
			}
			
			
			//for every chain, update
			for(int j = 0; j < nChain; j++){				
				Move move = chooseMove();
				move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat(),0);
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
			double acceptRatio = chains.get(j).getHeat()*chains.get(k).getLikelihood() + chains.get(k).getHeat()*chains.get(j).getLikelihood() - chains.get(j).getHeat()*chains.get(j).getLikelihood() - chains.get(k).getHeat()*chains.get(k).getLikelihood();

			
			
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

		//number of ancestors for each sampled node
		String numAncString = "";
		
		String toWrite = "";
		
		for(int i=0; i<currPedigree.numIndiv; i++){
			
			//num ancestors
			List<Node> anc = new ArrayList<Node>();
			currPedigree.getNode(i).getAncestors(anc);
			int nAnc = anc.size();
			
			numAncString += nAnc;
			
			for(int j=i+1; j<currPedigree.numIndiv; j++){
				
				Path rel = currPedigree.getRelationships()[i][j];
				 
				
				//writer.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, rel.getUp(), rel.getDown(), rel.getNumVisit()));
				
				//TODO for hastings test
				toWrite += String.format("%d%d%d", rel.getUp(), rel.getDown(), rel.getNumVisit());
				
			}
			
			
		}
		
		
		
		toWrite += numAncString;
		
		//TODO testing
		writer.write(toWrite+"\n");
		
		
		//multiplier and likelihood
		if(!ped2info.containsKey(toWrite)){			
			
			PedInfo info = new PedInfo();
			info.lkhd = currPedigree.getLogLikelihood();
			info.multiplier = computeMultiplier(currPedigree);
			
			ped2info.put(toWrite, info);		
		}

		
		
		//count
		if(!ped2info.containsKey(toWrite)){
			PedInfo info = new PedInfo();
			info.count = 1;
			ped2info.put(toWrite, info);
			
		}
		else{
			PedInfo info = ped2info.get(toWrite);
			info.count = info.count + 1;
			ped2info.put(toWrite, info);
		}
		
		
	

	}
	
	

	private void writeFamFile(Pedigree currPedigree){
		

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
		//currPedigree.clearVisit();
		//if(currPedigree.sexLocked(ind)==false) sex = "4";
			
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

	

	private int computeMultiplier(Pedigree currPedigree){
		
		
		int toReturn = 1;
		
		currPedigree.clearVisit();
		
		for(int i=0; i<currPedigree.getNActiveNodes(); i++){
			
			//skip visited cluster
			Node x = currPedigree.getNode(i);
			if(x.getNumVisit()>0) continue;
			
			
			//get cluster
			List<Node> cluster = new ArrayList<Node>();
			x.getConnectedNodes(cluster);
			
			
			//depth
			int minDepth = currPedigree.maxDepth;
			int maxDepth = 0;
			int maxSampleDepth = 0;
	
			//clear cluster visits
			for(Node y : cluster) y.setNumVisit(0);
			


			for(Node y : cluster){
				
				//depth
				int yDepth = y.getDepth();
				
				if(yDepth < minDepth) minDepth = yDepth;
				if(yDepth > maxDepth) maxDepth = yDepth;
				
				if(y.sampled && maxSampleDepth < y.getDepth())
					maxSampleDepth = y.getDepth();
				
				//only 1 parent: not sex locked
				if(y.getParents().size()==1){
					
					Node p = y.getParents().get(0);
					
					if(p.getNumVisit()>0) continue;
					
					if(!p.sampled)
						toReturn *= 2;
					
					p.setNumVisit(1);
					
				}
				
				//2 parents
				else if(y.getParents().size()==2){
					
					Node p1 = y.getParents().get(0);
					Node p2 = y.getParents().get(1);
					
					
					if(!p1.sampled && !p2.sampled && p1.getNumVisit()==0 && p2.getNumVisit()==0){
						
						int nCommonChildren = 1 + currPedigree.getFullSibs(y).size();
						
						if(p1.getParents().size()!=0 || p2.getParents().size()!=0 || p1.getChildren().size()!=nCommonChildren || p2.getChildren().size()!=nCommonChildren)
							toReturn *= 2;
						
					}
					
					
					
					p1.setNumVisit(1);
					p2.setNumVisit(1);
					
					
				}

				
			}
			
			for(Node y : cluster) y.setNumVisit(1);
			
			

			
			int depthFactor = 1;
			int k=1;
			while(k+maxDepth < currPedigree.maxDepth+1){ //going up
				if(maxSampleDepth+k <= currPedigree.maxSampleDepth){
					depthFactor++;
					k++;
				}
				else
					break;
			}
			
			k=1;
			while(minDepth-k >= 0){
				depthFactor++;
				k++;
			}

			
			toReturn = toReturn * depthFactor;
			
			
			
		}
		
		
		
		return toReturn;
			
		
		
	}
	
	
	public static double acceptanceRatio(double newLkhd, double oldLkhd, double oldToNew, double newToOld, double heat){
		
		
		if(newLkhd > 0){
			System.out.println("Positive lkhd");
		}
		
		
		
		if(oldToNew==Double.NEGATIVE_INFINITY || newToOld==Double.NEGATIVE_INFINITY || oldToNew==Double.POSITIVE_INFINITY || newToOld==Double.POSITIVE_INFINITY ){
			System.out.println("proposal prob inifinity");
		}
		
		
		//TODO testing
		heat = 0;
		
		
		
		double acceptRatio = heat * (newLkhd - oldLkhd) + newToOld - oldToNew;
		

		
		if(acceptRatio > 0){
			return 1;
		}
		else{
			return Math.exp(acceptRatio);
		}

	}
	
	

	
	
}

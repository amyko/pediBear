package mcmc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
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
	final static int tuneInterval = 50000;
	final static int maxDelta = 1;
	final static int minDelta = 0;

	String outPath; 
	final List<Chain> chains;
	final int burnIn;
	final int runLength;
	final int sampleRate;
	final Move[] moves; 
	final PrintWriter lkhdWriter;
	final PrintWriter countWriter; // label + count
	final PrintWriter NeWriter;
	final PrintWriter pairWriter;//count of each relationship (FS, HS, UR, FC, hC) for each pair (i,j)
	//final PrintWriter famWriter;// label + fam
	//final PrintWriter thetaWriter;
	final Random rGen;
	final int nChain;
	public int nSwapSuccess;
	public int nSwapAttempt;
	public int coldChain;
	public int swapInterval;
	boolean tuned;
	int nSwaps;
	
	private double deltaT;	
	
	//FOR SAMPLING
	private int missingParentCounter = 0;
	public Map<String, PedInfo> ped2info = new HashMap<String, PedInfo>();
	private final List<Node> anc = new ArrayList<Node>();
	public int[][][] pair_results; // (id1, id2) -> count(FS, HS, UR, FC, HC)
	public Map<Integer, Integer> Ne_count;

	


	//TODO parallelize 
	public MCMCMC(List<Chain> chains, double deltaT, Move[] moves, int burnIn, int runLength, int sampleRate, int swapInterval, int nSwaps, Random rGen, String outPath) throws IOException{

		this.outPath = outPath;
		this.chains = chains;
		this.deltaT = deltaT;
		this.burnIn = burnIn;
		this.runLength = runLength;
		this.sampleRate = sampleRate;
		this.moves = moves;		
		this.lkhdWriter = DataParser.openWriter(outPath+".lkhd");
		//TODO testing
		countWriter = DataParser.openWriter("/Users/amy/eclipse-workspace/mcmc/simulations2/sample.count");
		//this.countWriter = DataParser.openWriter(outPath+".count");
		this.NeWriter = DataParser.openWriter(outPath+".Ne");
		//this.thetaWriter = DataParser.openWriter(outPath+".theta");
		//this.famWriter = DataParser.openWriter(outPath+".fam");
		this.pairWriter = DataParser.openWriter(outPath+".pairAssignment");
		this.rGen = rGen;
		this.nChain = chains.size();
		this.nSwapAttempt = 0;
		this.nSwapSuccess = 0;
		this.coldChain = nChain - 1;
		this.swapInterval = swapInterval;
		this.tuned = false;
		this.nSwaps = nSwaps;
		
		int n = chains.get(0).getPedigree().numIndiv;
		pair_results = new int[n][n][5];
		Ne_count = new HashMap<Integer, Integer>();
		
		//header for .theta file
		//thetaWriter.print("N\talpha\tbeta\tNe\n");
		
		
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
					move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat());
						
					
					
									
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
	
	
	
	
	
	
	private void runSample() throws IOException{
		
		System.out.println("Sampling...");
		
		//now start sampling
		for(int i = 0; i < runLength; i++){
			
			if((i%100000) == 0) {
				System.out.println(i);
			}
			
			
			//sample from cold chain
			if(i % sampleRate == 0){
				
				sample(chains.get(this.coldChain).getPedigree());
				convergence(chains.get(this.coldChain).getPedigree());
				//writeTheta(thetaWriter, chains.get(this.coldChain).getPedigree());
			
			}
			
			
			//for every chain, update
			for(int j = 0; j < nChain; j++){				
		
				
				Move move = chooseMove();
				move.mcmcMove(chains.get(j).getPedigree(), chains.get(j).getHeat());
				
			}
			
			
			//swap chains if necessary
			if(i%swapInterval==0)
				swapStates();
			
			
	
		}
		
		//TODO
		//write counts
		System.out.println(String.format("%.4f", chains.get(coldChain).getPedigree().getAlpha()));
		System.out.println(String.format("%.4f", chains.get(coldChain).getPedigree().getBeta()));
		System.out.println(String.format("%d", chains.get(coldChain).getPedigree().getN()));
		writeCounts();

		//close outfile
		lkhdWriter.close();
		countWriter.close();
		NeWriter.close();
		//thetaWriter.close();
		//famWriter.close();
		pairWriter.close();

		
		
	}
	
	
	private void writeCounts(){
		
		
		//write pair counts
		for(int i=0; i<pair_results.length; i++) {
			
			for (int j = i+1; j<pair_results.length; j++) {
				
				pairWriter.write(String.format("%d %d %d %d %d %d %d\n", i, j, pair_results[i][j][0], pair_results[i][j][1], pair_results[i][j][2] ,pair_results[i][j][3], pair_results[i][j][4]));
				
			}
			
		}
		
		
		//write counts of Ne
		List<Integer> sortedKeys = new ArrayList<Integer>(Ne_count.keySet());
		Collections.sort(sortedKeys);
		for(int pop : sortedKeys) {
			NeWriter.write(String.format("%d %d\n", pop, Ne_count.get(pop)));
		}
		
		
		//write pedigree label
		for(String key : ped2info.keySet()){

			countWriter.write(String.format("%s\t%d\n", key, ped2info.get(key).count));	
	
		}
		
		
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
		
		//pop count
		int N = currPedigree.getNe();
		N = (int) Math.round(N / 10.0) * 10;// round to nearest 10
		int N_count = Ne_count.containsKey(N)? Ne_count.get(N) : 0;
		Ne_count.put(N, N_count+1);
		
		

		//number of ancestors for each sampled node
		String numAncString = "";
		StringBuilder sb = new StringBuilder(currPedigree.numIndiv * (currPedigree.numIndiv) / 2 * 3);
		
		//ind 1
		for(int i=0; i<currPedigree.numIndiv; i++){
			
			//num ancestors
			anc.clear();
			currPedigree.getNode(i).getAncestors(anc);
			numAncString += anc.size();
			
			//ind 2
			for(int j=i+1; j<currPedigree.numIndiv; j++){
				
				Path rel = currPedigree.getRelationships()[i][j];			 
				//toWrite += String.format("%d%d%d", rel.getUp(), rel.getDown(), rel.getNumVisit());
				sb.append(rel.getUp());
				sb.append(rel.getDown());
				sb.append(rel.getNumVisit());
				
				
				//tally pair
				int r = 2; // unrelated
				if(rel.getUp()==1) {
					r = rel.getNumVisit() == 2 ? 0 : 1;
				}
				else if(rel.getUp()==2) {
					r = rel.getNumVisit() == 2 ? 3 : 4;
				}
				pair_results[i][j][r]++;
				
				
			}
			
			
		}
		
		//TODO testing
		//write
		//pairWriter.write(String.format("%s\n", toWrite));
		//toWrite += numAncString;
		
		
		//if new pedigree, record fam, multiplier, and likelihood
		PedInfo info; 
		String toWrite = sb.toString();
		if(!ped2info.containsKey(toWrite)){			
			
			//lkhd & multiplier & initialize count
			info = new PedInfo();
			info.multiplier = computeMultiplier(currPedigree);
			info.count = 1;
			
			//fam
			missingParentCounter = 0;
			//writeFamFile(famWriter, chains.get(coldChain).getPedigree(), toWrite);
		
			//record pedigree
			ped2info.put(toWrite, info);
		
		}

	
		//lkhd and count
		info = ped2info.get(toWrite);
		info.count++;
		


	}
	
	

	private void writeFamFile(PrintWriter famWriter, Pedigree currPedigree, String pedLabel){
		
		//header 
		famWriter.write(String.format(">\t%s\t%f\n", pedLabel, currPedigree.getLogLikelihood()));
		
		//write family relationship
		famWriter.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		for(int i=0; i<currPedigree.getNActiveNodes(); i++){
			recordCranefootFam(currPedigree.getNode(i), currPedigree);
		}
		
		
	}
	
	
	private void convergence(Pedigree currPedigree){
		
		lkhdWriter.write(String.format("%f\n", currPedigree.getLogLikelihood()));
		
	}
	
	
	// write current parameters out to a file
	private void writeTheta(PrintWriter thetaWriter, Pedigree currPedigree) {
		
		thetaWriter.write(String.format("%d\t%f\t%f\t%d\n", currPedigree.getN(), currPedigree.getAlpha(), currPedigree.getBeta(), currPedigree.getNe()));

		
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
		//famWriter.write(String.format("%s\t%s\t%s\t%s\t%s\n", name, pa, ma, sex, sampleStatus));
		
		
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
			throw new RuntimeException("Positive lkhd");
		}
		
		
		
		if(oldToNew==Double.NEGATIVE_INFINITY || newToOld==Double.NEGATIVE_INFINITY || oldToNew==Double.POSITIVE_INFINITY || newToOld==Double.POSITIVE_INFINITY ){
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
	
	
	public static int getIndex(double d) {
		
		if(Math.abs(d - .001) < 1e-6) {
			return 0;
		}
		else if(Math.abs(d-.01) < 1e-6) {
			return 1;
		}
		else if(Math.abs(d-.1) < 1e-6) {
			return 2;
		}
		else return 3;
	}

	
	
}

package executables;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import dataStructures.Chain;
import dataStructures.Path;
import dataStructures.PedInfo;
import dataStructures.Pedigree;
import dataStructures.Relationship;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import likelihood.PreProcess;
import mcmc.MCMCMC;
import mcmcMoves.Contract;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.FStoPO;
import mcmcMoves.FStoSelf;
import mcmcMoves.FUtoHS;
import mcmcMoves.GPtoHS;
import mcmcMoves.HStoGP;
import mcmcMoves.HStoPO;
import mcmcMoves.HStoFU;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.NephewToUncle;
import mcmcMoves.OPtoPO;
import mcmcMoves.POtoFS;
import mcmcMoves.POtoHS;
import mcmcMoves.POtoOP;
import mcmcMoves.SelftoFS;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.Stretch;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import mcmcMoves.UncletoNephew;
import statistic.Accuracy;
import utility.DataParser;

public class RunMCMC{
	
	//MCMC parameter
	public static String fileName = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/simPed4/";
	public static String refPopFileName = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/mcmcTest/mcmcTest";
	public static String ageFileName = "";
	public static double maf = 0.01;
	public static double errorRate = 0.01;
	public static int maxDepth = 4;
	public static int sampleDepth = maxDepth;
	public static double back = 0.04;
	public static double startTemp = 100;
	public static double tempFact = 1.01;
	public static int iterPerTemp = 40000;
	public static int maxIter = 10000000;
	public static double conv = 1;
	public static int numIndiv = 18;
	public static double poissonMean = numIndiv;
	public static boolean conditional = true;
	public static int numRun = 1;
	public static int runLength = 1;
	public static int numThreads = 1;
	
	//misc
	public static int maxNumNodes = 200;
	public static Map<String, Double> name2age = null;
	public static Random rGen = new Random(102762);
	

	
	
	//relationships for likelihood computation
	private static List<Relationship> relationships = new ArrayList<Relationship>();


	public static void computeLikelihoods() throws IOException{
		
		//init relationships
		initRelationships();
		
		//compute info
		LDStreamPedMissing.writeLdOutfile(refPopFileName, fileName + ".info", back, conditional);
		
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(errorRate, back, numIndiv);
		
		//compute marginal likelihood
		System.out.println("Computing marginal likelihoods...");
		int[] indCols = new int[numIndiv];
		for(int i=0; i<numIndiv; i++) indCols[i] = i;
		computeMarginals(core, fileName, indCols);
		
		//compute pairwise likelihood
		System.out.println("Computing pairwise likelihoods...");
		computePairwise(core, fileName, indCols, relationships);
		

		
	}
	
	
	
	public static void initRelationships(){
		
		relationships.add(new Relationship(15d, new double[] {1d, 0d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
		
		// depth = 1 relationship
		relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
		relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib	
		relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
		
		//depth = 2 relationships
		if(maxDepth >= 2){
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
		}
		
		//depth = 3 relationships
		if(maxDepth >= 3){
			relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1)})); 
			relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1)}));
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
			relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2)})); 
		}
		
		//depth = 4 relationships 
		if(maxDepth >= 4){
			relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1)}));
			relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1)})); 
			relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2)}));
			relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2)}));
		
		}
		
	}
	

	public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids) throws IOException{
		
		//compute
		double[] marginals = core.computeMarginal(fileName+".tped", fileName+".info", ids);

		//open outfile
		PrintWriter writer = DataParser.openWriter(fileName+".marginal");				
		
		//write marginals to file
		for (double item : marginals){
			writer.write(String.format("%f\n",item));
		}
		writer.close();
	}
	

	public static void computePairwise(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, List<Relationship> relationships) throws IOException{
		
		int numIndiv = ids.length;
		int numRel = relationships.size();
		
		//likelihood
		PrintWriter writer = DataParser.openWriter(fileName+".pairwise");
			
		double[][][] lkhd = core.forwardAlgorithm(fileName+".tped", fileName+".info", ids, relationships);
		

		//write to file
		for(int k=0; k<numRel; k++){
			
			for(Path path : relationships.get(k).getAllPaths()){

				String header = String.format(">\t%d\t%d\t%d\n",path.getUp(), path.getDown(), path.getNumVisit());
				writer.write(header);
				//DataParser.writeMatrixToFile(writer, lkhd[k]);
				
				for(int i=0; i<numIndiv; i++){
					for(int j=i+1; j<numIndiv; j++){
						
						writer.write(String.format("%d\t%d\t%f\n", i,j,lkhd[k][i][j]));
						
					}
				}
				
				
				writer.flush();					
			}

		
		}
	
		//write to file
		writer.close();
		
	}
	
	//returns fid+iid --> age
	public static void setName2age() throws NumberFormatException, IOException{
		
	
		if(ageFileName.equals(""))
			name2age = null;
			
		else{
			//read age information
			BufferedReader reader = DataParser.openReader(ageFileName);
			
			name2age = new HashMap<String, Double>();
			
			String line;
			while((line=reader.readLine()) != null){
				
				String[] fields = line.split("\\s+");
				if(fields.length < 3) continue; //skip if age was not given for this fii+iid
				
				name2age.put(fields[0]+fields[1], Double.parseDouble(fields[2]));
				
			}
		}
		
		
	}
	

	public static void runThreads(String myFile, String outfile) throws IOException{
		
		//Random rGen = new Random(102762);
		
		//arguments
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .02), new Split("split", .02),  
				new CutLink("cutLink", .05), new SplitLink("splitLink", .05), 
				new ShiftClusterLevel("shiftClusterLevel", .02),  new SwitchSex("switchSex", .02),  
				new FStoSelf("FStoSelf", .07), new SelftoFS("selfToFS", .07),
				new HStoGP("HStoGP", .06), new GPtoHS("GPtoHS", .06),	
				new UncletoNephew("uncleToNephew", .07), new NephewToUncle("nephewToUncle", .07),
				new SwapDescAnc("swapDescAnc", .02),
				new OPtoPO("OPtoPO", .02), new POtoOP("POtoOP", .02),
				new FStoPO("FStoPO", .05), new POtoFS("POtoFS", .05), //confounds with HS2FU
				new HStoPO("HStoPO", .05), new POtoHS("POtoHS", .06),
				new Contract("contract", .05), new Stretch("stretch", .05), //confounds with HS2PO, HS2FU?
				
				new HStoFU("HStoFU",.0), new FUtoHS("FUtoHS", .0), //confounds with POtoFU
				new Split2("split2", 0), new SwapDescAnc("swapDescAnc", .0),
				new SwapUp("swapUp", .0), new SwapDown("swapDown", .0),
				new CutOneLinkTwo("cutOneLinkTwo", .0), new CutTwoLinkOne("cutTwoLinkOne", .0)};
		
		
		double prob = 0d;
		for(Move m : moves){
			prob += m.getProb();
		}
		System.out.println(prob);
		
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(errorRate, back, numIndiv);

		//mcmc parameters
		int nChain = 7;
		int nBranch = 1;
		int burnIn = 10;
		int runLength = 100000;
		int sampleRate = 25;
		double deltaT = .5;
		int swapInterval = 1;
		int nSwaps = 1;
		

		
		/*
		//testing age info
		name2age = new HashMap<String, Double>();
		name2age.put("1_1", 1d); 
		name2age.put("1_2", 3d);
		name2age.put("1_3", 2d);
		name2age.put("1_4", 1d);
		name2age.put("1_5", 2d);
		name2age.put("1_6", 2d);
		name2age.put("1_7", 1d);
		name2age.put("1_8", 2d);
		name2age.put("1_9", 3d);
		name2age.put("1_10", 2d);
		*/
		
		
		
		
		
		
		//init chains
		List<Chain> chains = new ArrayList<Chain>(nChain);

		if(nChain>1){
		
			int currIdx = 0;
			Pedigree ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age);
			Chain superHeatedChain = new Chain(nChain-1, ped);
			superHeatedChain.setHeat(deltaT);
			chains.add(superHeatedChain);
			currIdx++;
			
			
			for(int branch=0; branch < nBranch; branch++){
			
				for(int chain=nChain-2; chain >= 0; chain--){
					
					ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age);
					Chain myChain = new Chain(chain, ped);
					
					//temp
					myChain.setHeat(deltaT);
					
					
					//neighbors				
					if(chain==nChain-2){ //next to super heated chain
						superHeatedChain.addNeighbor(currIdx);
						myChain.addNeighbor(0);
					}
					else{
						myChain.addNeighbor(currIdx-1);
						chains.get(currIdx-1).addNeighbor(currIdx);
					}
					
					
					//add chain
					chains.add(myChain);
					currIdx++;
		
				}
			}
			
			
			
			/*
			//connect cold chains together
			for(int i=0; i<nBranch; i++){
				
				int iIdx = (i+1) * (nChain-1);
				Chain coldChain1 = chains.get(iIdx);
				
				for(int j=i+1; j<nBranch; j++){
					
					int jIdx = (j+1) * (nChain-1);
					Chain coldChain2 = chains.get(jIdx);
					
					coldChain1.addNeighbor(jIdx);
					coldChain2.addNeighbor(iIdx);
					
				}
				
			}
			*/
			
		}
		
		else{
			Pedigree ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age);
			Chain myChain = new Chain(0, ped);
			
			//temp
			myChain.setHeat(deltaT);
			chains.add(myChain);
		}
		
		
		MCMCMC mcmcmc = new MCMCMC(chains, deltaT, moves, burnIn, runLength, sampleRate, swapInterval, nSwaps, rGen, outfile);
		mcmcmc.run();
		
		
		
		/*
		//print counts
		for(PedInfo info : mcmcmc.ped2info.values()){
			System.out.println(String.format("%f %d",info.lkhd, info.count));
		}
		
		//check relative occupancy
		String[] peds = getTwoPeds(mcmcmc);
		//System.out.println(peds[0]);
		//System.out.println(peds[1]);
		if(peds[0]!=null && peds[1]!=null)
			checkRelativeOccupancy(outfile, mcmcmc, peds[0], peds[1]);
		System.out.println(mcmcmc.bestLkhd);
		*/
		
		

		
	}
	
	
	//validate output
	public static void validate(String outfile) throws IOException{
		
		//read target counts
		int numTarget = DataParser.countLines(outfile+".target");
		String[] target = new String[numTarget];
		int[] expectedCount = new int[numTarget];
		
		BufferedReader targetReader = DataParser.openReader(outfile+".target");
		String line;
		int lineNum=0;
		while((line=targetReader.readLine())!=null){
			
			String[] fields = line.split("#");
			target[lineNum] = fields[0]+"\t";
			expectedCount[lineNum++] = Integer.parseInt(fields[1].split("\t")[1]);
			
			
		}
		targetReader.close();
		
		

		//count
		int[] counts = new int[target.length];
		
		BufferedReader reader = DataParser.openReader(outfile+".pair");
		
		while((line=reader.readLine())!=null){
			
			if(line.charAt(0)=='>') continue;
			
			String myLine = line.split("\n")[0]+"\t";
			
			
			for(int i=0; i<target.length; i++){
				if(target[i].equals(myLine)){
					counts[i]++;
					break;
				}
			}

		}
		
		
		PrintWriter writer = DataParser.openWriter(outfile+".stat");
		
		for(int i=0; i<counts.length; i++){
			writer.write(String.format("%d\t%.2f\n", expectedCount[i], counts[i]/1.0));
		}
		
		writer.close();
		
		
	}
	
	//return two most likely
	public static String[] getTwoPeds(MCMCMC mcmcmc){
		
		
		double lkhd1 = Double.NEGATIVE_INFINITY;
		double lkhd2 = Double.NEGATIVE_INFINITY;
		String ped1 = null;
		String ped2 = null;
		
		for(String key : mcmcmc.ped2info.keySet()){
			
			PedInfo info = mcmcmc.ped2info.get(key);
			
			
			if(info.lkhd > lkhd1){
				lkhd2 = lkhd1;
				ped2 = ped1;
				lkhd1 = info.lkhd;
				ped1 = key;
			}
			else if(info.lkhd > lkhd2){
				lkhd2 = info.lkhd;
				ped2 = key;
			}
			
		}
		
		return new String[]{ped1, ped2};
		
	
		
	}
	
	
	public static void checkRelativeOccupancy(String outfile, MCMCMC mcmcmc, String ped1, String ped2) throws IOException{
		
		//compute multiplier
		PedInfo info1 = mcmcmc.ped2info.get(ped1);
		PedInfo info2 = mcmcmc.ped2info.get(ped2);
		
		System.out.println(String.format("%d %f %d", info1.count, info1.lkhd, info1.multiplier));
		System.out.println(String.format("%d %f %d", info2.count, info2.lkhd, info2.multiplier));
		
		double factor = Math.exp(info1.lkhd - info2.lkhd + Math.log(info1.multiplier) - Math.log(info2.multiplier));

		//TODO test
		//factor = Math.exp(Math.log(info1.multiplier) - Math.log(info2.multiplier));
		
		//record observed/expected
		BufferedReader reader = DataParser.openReader(outfile+".pair");
		PrintWriter writer = DataParser.openWriter(outfile+".conv");
		int lineCounter = 0;
		int count1 = 0;
		int count2 = 0;
		
		String line = "";
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			if(fields[0].equals(">")) continue;

			if(fields[0].equals(ped1)){
				count1++;
			}
			else if(fields[0].equals(ped2)){
				count2++;
			}
			
			lineCounter++;
			
			
			//recompute ratio
			if(lineCounter%1000==0){
				double expected = count2 * factor;
				double ratio = expected / count1;
				writer.write(String.format("%d\t%.3f\n", lineCounter, ratio));
			}
			
			
		}
		
		writer.close();
		
		
	}
	
	
	public static void writeMap(PrintWriter writer, int t) throws NumberFormatException, IOException{
			
		//get mapAcc
		String truePath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim4.true";
		String outPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim5.0.pair";
		String pathToOmega = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/pathToOmega.txt";
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		double[][] mapAcc = Accuracy.mapAccuracy(outPath, truePath, numIndiv, numIndiv, pathToKinship);
		
		
		//header
		writer.write(String.format(">\t%d\n", t));
		
		int count = 0;
		
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				writer.write(String.format("%d\t%d\t%d\n", i, j, (int)mapAcc[i][j]));
				
				count += (int)mapAcc[i][j];
				
				
			}
		}
		
		
		//TODO testing
		System.out.println(count);
		

		//flush
		writer.flush();

		
		
		
	}
	
	
	public static void main(String[] args) throws IOException{
	
		
		//open output file
		//PrintWriter writer = DataParser.openWriter("/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim4.mcmc.5chains.mapAcc");
		PrintWriter writer = DataParser.openWriter("/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/testing");

		//run
		for(int i=0; i<100; i++){
			
			System.out.println(i);
			
			for(int j=0; j<1; j++){
			
				String myFile = fileName + "sim4." + i;
				String outfile = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim5."+j;

				
				runThreads(myFile, outfile);
				
				//validate(outfile);
				
				writeMap(writer, i);
			}
			
		}

		
		
		writer.close();


		
		System.out.println("DONE");
		
		
		
		
		
		
		
	}



	

}
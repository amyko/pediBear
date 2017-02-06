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
import dataStructures.Pedigree;
import dataStructures.Relationship;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import likelihood.PreProcess;
import mcmc.MCMCMC;
import mcmcMoves.Contract;
import mcmcMoves.CousinToGreatUncle;
import mcmcMoves.CousinToHalfUncle;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.FStoPO;
import mcmcMoves.FullUncletoHalfSibs;
import mcmcMoves.GreatUncleToCousin;
import mcmcMoves.HStoPO;
import mcmcMoves.HalfCousinToHalfGreatUncle;
import mcmcMoves.HalfGreatUncleToHalfCousin;
import mcmcMoves.HalfSibstoFullUncle;
import mcmcMoves.HalfUncleToCousin;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.OPtoPO;
import mcmcMoves.POtoFS;
import mcmcMoves.POtoHS;
import mcmcMoves.POtoOP;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.Stretch;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import statistic.Accuracy;
import utility.DataParser;

public class RunMCMC{
	
	//MCMC parameter
	public static String fileName = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/simPed2/";
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
	public static int numIndiv = 3;
	public static double poissonMean = numIndiv;
	public static boolean conditional = true;
	public static int numRun = 1;
	public static int runLength = 1;
	public static int numThreads = 1;
	
	//misc
	public static int maxNumNodes = 200;
	public static Map<String, Double> name2age = null;

	
	
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
		
		//arguments
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .05), new Split("split", .05), new Split2("split2", .02), 
				new CutLink("cutLink", .05), new SplitLink("splitLink", .05), 
				new ShiftClusterLevel("shiftClusterLevel", .02),  new SwitchSex("switchSex", .03), 
				new SwapUp("swapUp", .03), new SwapDown("swapDown", .03), new SwapDescAnc("swapDescAnc", .02), 
				new CutOneLinkTwo("cutOneLinkTwo", .05), new CutTwoLinkOne("cutTwoLinkOne", .05),
				new Contract("contract", .05), new Stretch("stretch", .05),
				new FStoPO("FStoPO", .05), new POtoFS("POtoFS", .05),  
				new OPtoPO("OPtoPO", .05), new POtoOP("POtoOP", .05),
				new HStoPO("HStoPO", .05), new POtoHS("POtoHS", .05),
				new HalfSibstoFullUncle("halfSibsToFullUncle", .05), new FullUncletoHalfSibs("fullUncleToHalfSibs", .05),
				new HalfUncleToCousin("halfUncleToCousin_REDUNDANT", 0), new CousinToHalfUncle("cousinToHalfUncle_REDUNDANT", 0),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0),  
				new CousinToGreatUncle("cousinToGreatUncle", 0), new GreatUncleToCousin("greatUncleToCousin", 0)};
		
		double prob = 0d;
		for(Move m : moves){
			prob += m.getProb();
		}
		System.out.println(prob);
		
		Random rGen = new Random(123);
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(errorRate, back, numIndiv);

		//mcmc parameters
		int nChain = 5;
		int nBranch = 1;
		int burnIn = 10000;
		int runLength = 2000000;
		int sampleRate = 50;
		double deltaT = .2;
		int swapInterval = 1;
		int nSwaps = 1;
		

		
		/*
		//testing age info
		name2age = new HashMap<String, Double>();
		name2age.put("1_1", 5d); 
		name2age.put("1_2", 4d);
		name2age.put("1_3", 4d);
		name2age.put("1_4", 3d);
		name2age.put("1_5", 2d);
		name2age.put("1_6", 1d);
		name2age.put("1_7", 3d);
		name2age.put("1_8", 2d);
		name2age.put("1_9", 1d);
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
			
			String myLine = line.split("\t\n")[0];
			
			
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
	
	
	
	
	
	public static void writeMap(PrintWriter writer, int t) throws NumberFormatException, IOException{
			
		//get mapAcc
		String truePath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim4.true";
		String outPath = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim4.pair";
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
		//PrintWriter writer = DataParser.openWriter("/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim4.mcmc.3chains.mapAcc");
		//PrintWriter writer = DataParser.openWriter("/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/testing");

		//run
		for(int i=0; i<1; i++){
			
			System.out.println(i);
			
			for(int j=0; j<1; j++){
			
				String myFile = fileName + "sim2." + i;
				String outfile = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/results/sim2."+j;

				
				runThreads(myFile, outfile);
				
				
				validate(outfile);
				
				//writeMap(writer, i);
			}
			
		}

		
		
		//writer.close();


		
		System.out.println("DONE");
		
		
		
		
		
		
		
	}



	

}
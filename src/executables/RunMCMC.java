package executables;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import Unused.CutOneLinkTwo;
import Unused.CutTwoLinkOne;
import Unused.FUtoHS;
import Unused.HStoFU;
import dataStructures.Chain;
import dataStructures.Path;
import dataStructures.PedInfo;
import dataStructures.Pedigree;
import dataStructures.Relationship;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import likelihood.PreProcess;
import likelihood.Prior;
import likelihood.PriorNoLoop;
import mcmc.MCMCMC;
import mcmcMoves.Contract;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.FStoPO;
import mcmcMoves.FStoSelf;
import mcmcMoves.GPtoHS;
import mcmcMoves.HStoGP;
import mcmcMoves.HStoPO;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.NephewtoUncle;
import mcmcMoves.OPtoPO;
import mcmcMoves.POtoFS;
import mcmcMoves.POtoHS;
import mcmcMoves.POtoOP;
import mcmcMoves.SelftoFS;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.Split;
import mcmcMoves.SplitLink;
import mcmcMoves.Stretch;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchChildAncWithOppositeSex;
import mcmcMoves.SwitchSex;
import mcmcMoves.UncletoNephew;
import mcmcMoves.UpdateAlpha;
import mcmcMoves.UpdateBeta;
import mcmcMoves.UpdateNe;
import simulator.SimulatePedigreeUnderPrior;
import statistic.Accuracy;
import statistic.CompareWithColony;
import utility.ArrayUtility;
import utility.DataParser;
import utility.MapUtility;

public class RunMCMC{
	
	//MCMC parameter
	public static String fileName = "/Users/amy/eclipse-workspace/mcmc/simulations2/";
	public static String refPopFileName = "/Users/amy/eclipse-workspace/mcmc/simulations/";
	public static String ageFileName = "";
	public static double maf = 0.01;
	public static double errorRate = 0.01;
	public static int maxDepth = 2; //only simulated up to depth=3
	public static int sampleDepth = 0;
	public static double back = 0.04;
	public static double startTemp = 100;
	public static double tempFact = 1.01;
	public static int iterPerTemp = 40000;
	public static int maxIter = 10000000;
	public static double conv = .1;
	public static int numIndiv = 20;
	public static double poissonMean = numIndiv;
	public static boolean conditional = true;
	public static int numRun = 1;
	//public static int runLength = 1;
	public static int numThreads = 1;
	public static double credibleInterval = .95;
	
	
	//MCMC parameters
	public static int nChain = 1;
	public static int nBranch = 1;
	public static int burnIn = 1000000;
	public static int runLength = 2000000;
	public static int sampleRate = 50;
	public static double deltaT = .5;
	public static int swapInterval = 1;
	public static int nSwaps = 1;
	
	//misc
	public static int maxNumNodes = 1000;
	public static Map<String, Double> name2age = null;
	public static Random rGen = new Random(102574);
	public static int minN = 5;
	public static int maxN = 3000;
	public static int stepSize = 50;
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
	

	public static MCMCMC runThreads(String myFile, String outfile, double alpha_sd, double alpha_min, double alpha_max, double beta_sd, double beta_min, double beta_max) throws IOException{
		
		Random rGen = new Random(1025742);
		
		/*
		//arguments
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .05), new Split("split", .05),  
				new ShiftClusterLevel("shiftClusterLevel", .05),  new SwitchSex("switchSex", .05),  
				new FStoSelf("fs2self", .05), new SelftoFS("self2fs", .05),
				new HStoGP("hs2gp", .05), new GPtoHS("gp2hs", .05),	
				new UncletoNephew("uncle2nephew", .05), new NephewtoUncle("nephew2uncle", .05),
				new SwapDescAnc("swapDescAnc", .05),
				new OPtoPO("OPtoPO", .02), new POtoOP("POtoOP", .02),
				new FStoPO("FStoPO", .05), new POtoFS("POtoFS", .05), //confounds with HS2FU
				new HStoPO("HStoPO", .03), new POtoHS("POtoHS", .03),
				new Contract("contract", .05), new Stretch("stretch", .05), //confounds with HS2PO, HS2FU?
				
				new SwitchChildAncWithOppositeSex("switchChildParentWithOppositeSex", .04),
				new NoChange("noChange", .06),
				
				new CutLink("cutLink", .0), new SplitLink("splitLink", .0), //these don't work if donor node is deleted in link
				new HStoFU("HStoFU",.0), new FUtoHS("FUtoHS", .0), //confounds with POtoFU
				new SwapUp("swapUp", .0), new SwapDown("swapDown", .0),
				new CutOneLinkTwo("cutOneLinkTwo", .0), new CutTwoLinkOne("cutTwoLinkOne", .0)};
		
		
		*/
		

		//single generation inference
		Move[] moves = new Move[]{new Link("link", .2), new Cut("cut", .2), 
				new Split("split", .1),  new SwitchSex("switchSex", .1),  
				new FStoSelf("fs2self", .1), new SelftoFS("self2fs", .1),
				new UpdateAlpha("updateAlpha", .05, alpha_sd, alpha_min, alpha_max), new UpdateBeta("updateBeta", .05, beta_sd, beta_min, beta_max), new UpdateNe("updateNe", .1, minN, maxN, stepSize)};
		
		
		double prob = 0d;
		for(Move m : moves){
			prob += m.getProb();
		}
		System.out.println(prob);
		
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(errorRate, back, numIndiv);
		
		// TODO testing prior no loop
		PriorNoLoop priorCalculatorNoLoop = new PriorNoLoop(rGen, maxDepth);
		Prior priorCalculator = new Prior(rGen, maxDepth);
		

		
		/*
		//testing age info
		name2age = new HashMap<String, Double>();
		name2age.put("1_1", 10d); 
		name2age.put("1_2", 9d);
		name2age.put("1_3", 9d);
		name2age.put("1_4", 8d);
		name2age.put("1_5", 7d);
		name2age.put("1_6", 6d);
		name2age.put("1_7", 5d);
		name2age.put("1_8", 4d);
		name2age.put("1_9", 3d);
		//name2age.put("1_10", 2d);
		*/
		
		
		
		
		
		
		//init chains
		List<Chain> chains = new ArrayList<Chain>(nChain);

		if(nChain>1){
		
			int currIdx = 0;
			Pedigree ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age, 
					minN, maxN, (alpha_min+alpha_max)/2, (beta_min+beta_max)/2, priorCalculator, priorCalculatorNoLoop);
			Chain superHeatedChain = new Chain(nChain-1, ped);
			superHeatedChain.setHeat(deltaT);
			chains.add(superHeatedChain);
			currIdx++;
			
			
			for(int branch=0; branch < nBranch; branch++){
			
				for(int chain=nChain-2; chain >= 0; chain--){
					
					ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age, 
							minN, maxN, (alpha_min+alpha_max)/2, (beta_min+beta_max)/2, priorCalculator, priorCalculatorNoLoop);
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
			Pedigree ped = new Pedigree(myFile, core, maxDepth, sampleDepth, rGen, maxNumNodes, poissonMean, numIndiv, name2age, 
					minN, maxN, (alpha_min+alpha_max)/2, (beta_min+beta_max)/2, priorCalculator, priorCalculatorNoLoop);
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
		
		return mcmcmc;

		
	}
	
	
	//validate output
	public static void validate(String outfile, String targetFile) throws IOException{
		
		//read target counts
		int numTarget = DataParser.countLines(targetFile);
		String[] target = new String[numTarget];
		int[] expectedCount = new int[numTarget];
		
		BufferedReader targetReader = DataParser.openReader(targetFile);
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
	
	/*
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
	*/
	
	
	public static void writeAcc(PrintWriter writer, String truePed, String countPath, int t) throws NumberFormatException, IOException{
			

		double denom = runLength / sampleRate;
		
		
		//in MAP?
		String[] mapResult = Accuracy.inMAP(countPath, truePed);
		//double percentMap = Integer.parseInt(mapResult[1]) / denom;
		
		
		//in credible interval?
		String[] ciResult = Accuracy.inCIPair(countPath, truePed, credibleInterval, denom);
		
		//in credible interval except false positives?
		
		writer.write(String.format("%d\t%s\t%s\t%s\n", t, ciResult[0], ciResult[1], ciResult[2]));
		
		//writer.write(String.format("%d\t%s\t%.3f\t%s\t%s\n", t, mapResult[0], percentMap, ciResult[0], ciResult[1]));
		
		
	
		//flush
		writer.flush();

		
		
		
	}
	
	public static String getTruePed(String truePath) throws IOException{
		
		BufferedReader reader = DataParser.openReader(truePath);
		
		String toReturn = reader.readLine().split("\\s")[0];
		reader.close();
		
		return toReturn;
		
	}
	
	
	public static String getMapPed(String inPath) throws NumberFormatException, IOException{
		//get MAP ped
		double bestLkhd = Double.NEGATIVE_INFINITY;
		String bestPed = "";
		
		
		BufferedReader reader = DataParser.openReader(inPath);
		
		String line;
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			double currLkhd = Double.parseDouble(fields[1]);
			
			if(currLkhd > bestLkhd){
				bestLkhd = currLkhd;
				bestPed = fields[0];
			}
			
			
		}
		reader.close();
		
		return bestPed;
		
	}
	
	
	//write cranefoot fam file for MAP estimate
	public static void writeMapFam(String outfile) throws NumberFormatException, IOException{
		
		String bestPed = getMapPed(outfile+".count");
		
		
		//write fam
		BufferedReader reader = DataParser.openReader(outfile+".fam");
		PrintWriter writer = DataParser.openWriter(outfile+".mapFam");
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			if(!fields[0].equals(">")) 
				continue;
			
			if(fields[1].equals(bestPed))
				break;
			
		}
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			if(fields[0].equals(">")) 
				break;
			
			writer.write(line+"\n");
			
		}
		
		reader.close();
		writer.close();
		
	}
	
	
	
	//get MAP estimate and compare pairwise
	public static void writePairAcc(PrintWriter writer, String fileName, int n, Map<Path,double[]> path2omega, int t) throws NumberFormatException, IOException{
		
		//get MAP ped
		String bestPed = getMapPed(fileName+".count");
		
		//name to index
		Map<String, Integer> name2idx = new HashMap<String, Integer>();
		BufferedReader reader = DataParser.openReader(String.format("%s.%d.tfam", fileName, t));
		
		String line;
		int idx = 0;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			name2idx.put(fields[1], idx);
			
			idx++;	
		}
		
		reader.close();
		
		
		//get true relationship matrix
		Path[][] trueRel = new Path[n][n];
		reader = DataParser.openReader(String.format("%s.%d.true", fileName, t));
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			int i = name2idx.get(fields[0]);
			int j = name2idx.get(fields[1]);
			Path rel = new Path(Integer.parseInt(fields[2]), Integer.parseInt(fields[3]), Integer.parseInt(fields[4]));
			
			trueRel[i][j] = rel;
			trueRel[j][i] = rel;
			
		}
		
		reader.close();
		
		
		//evaluate pairwise accuracy
		//header
		writer.write(String.format(">\t%d\n", t));
		
		idx = 0;
		for(int i=0; i<n; i++){
			for(int j=i+1; j<n; j++){
				
				//true rel
				Path truth = trueRel[i][j];
				
				int up = Integer.parseInt(bestPed.charAt(idx)+"");
				int down = Integer.parseInt(bestPed.charAt(idx+1)+"");
				int nVisit = Integer.parseInt(bestPed.charAt(idx+2)+"");
				
				double[] inferredOmega = path2omega.get(new Path(up,down,nVisit));
				double[] trueOmega = path2omega.get(truth);
				
				int acc = 0;
				if(inferredOmega[0]==trueOmega[0] && inferredOmega[1]==trueOmega[1] && inferredOmega[2]==trueOmega[2]){
					acc = 1;
				}
				
				writer.write(String.format("%d\t%d\t%d\n", i,j,acc));
				
				idx+=3;
				
				//System.out.print(String.format("truth: %d %d %d\t", truth.getUp(), truth.getDown(), truth.getNumVisit()));
				//System.out.println(String.format("inferred: %d %d %d", up, down, nVisit));
				
				
				
			}
		}
		
		
		
	}
	
	
	/*
	public static void writePrior(PrintWriter writer, MCMCMC mcmcmc, int t){
		
		//header
		writer.write(String.format(">\t%d\n", t));
		
		//print prior
		int numSamples = runLength /sampleRate;
		for(int i=0; i<mcmcmc.logLkhdOfNe.length; i++){
			writer.write(String.format("%d %f\n", i*stepSize + minN, mcmcmc.logLkhdOfNe[i] - numSamples));
		}
		
		writer.flush();
		
	}
	*/
	
	public static void writePosteriorForParameter(PrintWriter writer, MCMCMC mcmcmc, int t){
		

		//header
		writer.write(String.format(">\t%d\n", t));
		
		//print 
		for(int i=minN; i<=maxN; i++){
			
			int index = i - minN;
			double totalPosterior = 0d;
			
			//for every tree sampled, get likelihood P(G, i)
			for(String key : mcmcmc.ped2info.keySet()){
				
				//double posterior = mcmcmc.ped2info.get(key).lkhd[index];
				double posterior = mcmcmc.ped2info.get(key).counts[index];
				
				if(posterior!=Double.NaN && !Double.isInfinite(posterior)){
					totalPosterior += posterior;
				}
					
				
			}
			
			//if this theta was never sampled, skip
			if(totalPosterior==0)
				totalPosterior = Double.NEGATIVE_INFINITY;
			
			writer.write(String.format("%d %f\n", i, totalPosterior));
			
			
		}
		
		writer.flush();
		
	}
	
	
	public static void writePosteriorForPedigree(PrintWriter writer, MCMCMC mcmcmc, int t){
		
		//header
		writer.write(String.format(">\t%d\n", t));
		

		//for every tree sampled, get likelihood P(G, i)
		for(String key : mcmcmc.ped2info.keySet()){
			
			double totalPosterior = 0;
			
			for(int i=0; i<mcmcmc.ped2info.get(key).lkhd.length; i++){
				
				double posterior = mcmcmc.ped2info.get(key).lkhd[i];
				if(posterior!=Double.NaN){
					totalPosterior += posterior;
				}
				
			}
			
			//if this theta was never sampled, skip
			if(totalPosterior==0)
				totalPosterior = Double.NEGATIVE_INFINITY;
			
			writer.write(String.format("%s %f\n", key, totalPosterior));
			
		}
		

			
			
		
	
		writer.flush();
		
	}
	
	
	
	public static void writePairCounts(String outPath, MCMCMC mcmcmc, String[] ids) throws IOException {
		
		
		//open file
		PrintWriter writer = DataParser.openWriter(outPath);
		
		for(int i = 0; i< mcmcmc.pair_results.length; i++) {
			
			for(int j = i+1; j < mcmcmc.pair_results.length; j++) {
				
				writer.write(String.format("%s %s ", ids[i], ids[j]));
				
				
				for(int k = 0; k < 5; k++) {
					
					writer.write(mcmcmc.pair_results[i][j][k] + " ");
					
				}
				
				writer.write("\n");
				
			}
			
		}
		
		
		writer.close();
		
		
	}
	
	
	
	public static void assign(String outPath, String inPath, double thresh) throws NumberFormatException, IOException {
		
		
		PrintWriter writer = DataParser.openWriter(outPath);
		BufferedReader reader = DataParser.openReader(inPath);
		
		String line;
		while((line = reader.readLine()) != null) {
			
			String[] fields = line.split("\\s");
			
			//find max vote idx
			int maxIdx = -1;
			int max = -1;
			int unrelCount = 0;
			for(int i=0; i<5; i++) { 
				int curr = Integer.parseInt(fields[i + 2]);
				if(curr > max) {
					maxIdx = i;
					max = curr;
				}
				if(i == 2)
					unrelCount = curr;
			}

			
			//if count / UR_count < c, assign as UR
			if(maxIdx !=2 && unrelCount > 0 && (double) max / unrelCount < thresh) {
				maxIdx = 2;
			}
			
			
			String assign = "";
			if(maxIdx == 0)
				assign = "FS";
			else if(maxIdx == 1)
				assign  = "HS";
			else if(maxIdx == 2)
				assign  = "UR";
			else if(maxIdx == 3)
				assign  = "FC";
			else
				assign  = "HC";

			//write results
			writer.write(String.format("%s %s %s\n", fields[0], fields[1], assign));
			
			
			
		}
		
		
		writer.close();
		
		
	}
	
	
	public static void writePop(String outPath, MCMCMC mcmcmc) throws IOException {
		
		//open file
		PrintWriter writer = DataParser.openWriter(outPath);
		
		for(int i = 0; i< mcmcmc.pop_results.length; i++) {
			
			writer.write(String.format("%d %d\n", 100*(i+1), mcmcmc.pop_results[i]));
			
		}
		
		
		writer.close();
		
		
	}
	
	
	
	public static void main(String[] args) throws IOException{
	
	
		/*
		/////// for frog data ////////
		
		//run mcmc
		System.out.println("Running MCMC");
		String outFile = String.format("/Users/amy/eclipse-workspace/mcmc/frogs/frogs.juv");
		String myFile = outFile;

		
		MCMCMC mcmcmc = runThreads(myFile, outFile);
		
		//write pairwise results to file
		System.out.println("Analyzing pairwise results");
		String [] ids = new String[numIndiv];
		BufferedReader reader = DataParser.openReader(myFile + ".tfam");
		for(int i=0; i<ids.length; i++) {
			ids[i] = reader.readLine().split("\\s")[1];
		}
		
		writePairCounts(outFile + ".pairCounts", mcmcmc, ids);
		writePop(outFile + ".Ne", mcmcmc);
			
		
		// assign relationship
		assign(outFile + ".pairAssignment", outFile + ".pairCounts", 1);
		 */
		
		
		
		//////// FOR SIMULATIONS /////////
		//file paths
		String fileName = "/Users/amy/eclipse-workspace/mcmc/simulations2/";
		String sim = "sample";
		
		//get truePed
		//Map<Path,double[]> path2omega = Accuracy.getPathToOmega(outPath + "pathToOmega.txt");
		
		double[][] totalAcc = new double[5][5];
		int T = 10;
		double thresh = 1; // threshold for calling a relative
		int[] pops = new int[T];
		double alpha_sd = .005; //.005
		double alpha_min = .001; //.001
		double alpha_max = 2; //2
		double beta_sd = .005; //.005
		double beta_min = .001; //.001
		double beta_max = 1; //1
		
		/*
		//TODO testing prior for true pedigree
		List<String> samples = new ArrayList<String>();
		for(int i=0; i<200; i++) {
			samples.add(String.format("0_%d", i));
		}
		new Pedigree(outPath+"gasbarra.fam", outPath+"test", samples, new Prior(rGen, maxDepth));
		*/

		//run
		for(int t=0; t<T; t++){
			
			System.out.println(t);
			
			String outFile = String.format("%s%s.%d", fileName, sim, t);
			String myFile = String.format("%s%s.%d", fileName, sim, t);
			String truePed = getTruePed(String.format("%s%s.%d.true", fileName, sim, t));
			
			MCMCMC mcmcmc = runThreads(myFile, outFile, alpha_sd, alpha_min, alpha_max, beta_sd, beta_min, beta_max);
		
			//accuracy matrix
			double[][] acc = CompareWithColony.accuracyMatrix(String.format("%ssample.%d.true", fileName, t), String.format("%ssample.%d.count", fileName, t), numIndiv, thresh);
			
			for(int i=0; i<acc.length; i++) {
				for(int j=0; j<acc[0].length; j++) {
					
					totalAcc[i][j] += acc[i][j];
					
				}
				
			}
			
			
			//best estiamte of Ne
			pops[t] = MapUtility.maxValue(mcmcmc.Ne_count);
			
			
			//writeAcc(writer, truePed, outFile+".count", t);
			//writeMapFam(outFile);
			
			//testing relative occupancy without likelihood
			//String targetFile = String.format("%s.%d.target", outFile, 0);
			//validate(outFile, targetFile);
			
			
			
		}

		//writer.close();
		//priorWriter.close();
		
		
		//print accuracy matrix
		for(int i=0; i<totalAcc.length; i++) {
			
			double rowSum = 0d;
			
			for(int j=0; j<totalAcc[0].length; j++) {
				
				rowSum += totalAcc[i][j];
				
			}
			
			for(int j=0; j<totalAcc[0].length; j++) {
				
				System.out.print(String.format("%d ", (int) totalAcc[i][j]));
				
			}
			System.out.println();
			
		}
		
		//print pops
		for(int x : pops)
			System.out.print(String.format("%d, ", x));
		System.out.println("\n");
		
		System.out.println("DONE");
		
		
		
	}



	

}
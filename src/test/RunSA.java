package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Random;

import statistic.Accuracy;
import utility.DataParser;
import dataStructures.Path;
import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import mcmc.SimulatedAnnealing;
import mcmcMoves.CousinToGreatUncle;
import mcmcMoves.CousinToHalfUncle;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.FStoPO;
import mcmcMoves.FStoSelf;
import mcmcMoves.GPtoHS;
import mcmcMoves.HStoGP;
import mcmcMoves.HStoPO;
import mcmcMoves.HalfCousinToHalfGreatUncle;
import mcmcMoves.HalfGreatUncleToHalfCousin;
import mcmcMoves.HalfUncleToCousin;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.NephewToUncle;
import mcmcMoves.OPtoPO;
import mcmcMoves.POtoFS;
import mcmcMoves.POtoHS;
import mcmcMoves.POtoOP;
import mcmcMoves.SelftoFS;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.GreatUncleToCousin;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.Contract;
import mcmcMoves.Stretch;
import mcmcMoves.HalfSibstoFullUncle;
import mcmcMoves.FullUncletoHalfSibs;
import mcmcMoves.CutShiftLink;
import mcmcMoves.UncletoNephew;



public class RunSA {

	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int maxDepth = 4;
		int sampleDepth = maxDepth;
		int numIndiv = 20;
		int totalIndiv = 20;
		double seqError = 0.01;
		int back = 30000;
		int maxNumNodes = 200;
		double prior = numIndiv;
		double stopThresh = 0;
		double beta = 30;
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmega = dir + "pathToOmega.txt";


		//SA parameters
		double[] heat = new double[300]; //200
		heat[0] = .5; //.1
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*1.01;
		System.out.println(heat[heat.length-1]);
		int coolingTime = 10000;
 		int runLength = 1;
 		int numRun = 6;
		Random rGen = new Random(1942083275L);
		
					
		//cooling schedule
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = coolingTime;
		}

		
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .02), new Split("split", .02), 
				new CutLink("cutLink", .05), new SplitLink("splitLink", .05), 
				new ShiftClusterLevel("shiftClusterLevel", .02),  new SwitchSex("switchSex", .02),  
				new FStoSelf("FStoSelf", .07), new SelftoFS("selfToFS", .07),
				new HStoGP("HStoGP", .06), new GPtoHS("GPtoHS", .06),	
				new UncletoNephew("uncleToNephew", .07), new NephewToUncle("nephewToUncle", .07),
				new SwapDescAnc("swapDescAnc", .02),
				new OPtoPO("OPtoPO", .02), new POtoOP("POtoOP", .02),
				new FStoPO("FStoPO", .05), new POtoFS("POtoFS", .05),
				new HStoPO("HStoPO", .05), new POtoHS("POtoHS", .06),
				new Contract("contract", .05), new Stretch("stretch", .05)};
		

		

		double mySum = 0d;
		for(Move mov : moves) mySum += mov.getProb();
		System.out.println(mySum);
		

		
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		
		
			
		//input files
		String testName = "sim2";
		String truePath = dir + "results/sim2.true";
		String accPath = String.format(dir + "results/testing");
		
		//true path
		Path[][] trueRel = Accuracy.getTruePath(truePath, totalIndiv);
		
		//open outfiles
		PrintWriter mapWriter = DataParser.openWriter(accPath+".mapAcc");
		PrintWriter distWriter = DataParser.openWriter(accPath+".kinshipDist");
		PrintWriter nWriter = DataParser.openWriter(accPath+".nCorrect");
		PrintWriter pairWriter = DataParser.openWriter(accPath+".pair");
			


		
		
		for(int t=0; t<100; t++){

			System.out.println(t);       
			
			String fileName = dir + "simPed2/"+testName+"."+t;
			String outDir = dir + "results/mcmc";
			
			double bestLkhd = Double.NEGATIVE_INFINITY;
			int bestRun = 0;
			double bestLkhdAdj = Double.NEGATIVE_INFINITY;
			
			for(int run=0; run<numRun; run++){
				
				

				
				//System.out.println(run);
				
				//initialize pedigree
				Pedigree ped = new Pedigree(fileName, core, maxDepth, sampleDepth, rGen, maxNumNodes, prior, numIndiv, beta);

			
				//initialize SA
				SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, outDir+"."+run, stopThresh);

				
				//run MCMC
				double startTime = System.nanoTime();
				sa.run();
				double endTime = System.nanoTime();

				double duration = (endTime - startTime)/1e9; 

				//System.out.println(String.format("Number of singletons: %d", ped.nSingletons[ped.curr]));
				System.out.println(String.format("Running time: %.1f seconds", duration));
				
				
				//if currlkhd is better than best, save mcmc output
				double currLkhd = ped.getLogLikelihood();
				if(currLkhd > bestLkhd){
					bestLkhd = currLkhd;
					bestRun = run;
					bestLkhdAdj = currLkhd - ped.getSingletonProb();
				}
				
				//sanity check
				//System.out.println(ped.likelihoodAllPedigrees());
				
				
			}
			
			//System.out.println(String.format("Best MCMC lkhd: %f  %f", bestLkhdAdj, bestLkhd));

			
			//////////////////////////////////
			//likelihood for true pedigree

			//copy relationship
			Pedigree ped = new Pedigree(fileName, core, maxDepth, sampleDepth, rGen, maxNumNodes, prior, numIndiv, beta);
			Path[][] mcmcmcRel = ped.getRelationships();
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = trueRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			ped.nSingletons[ped.curr] = 11;
			
			double trueLkhd = ped.likelihoodAllPedigrees();	
			double trueLkhdAdj = ped.likelihoodAllPedigrees() - ped.getSingletonProb();		
			//System.out.println(String.format("lkhd of true pedigree: %.2f %.2f", trueLkhdAdj, trueLkhd));
			
			
			
			////////////////////////////////////////
			//kinship accuracy
			double[][] mapAcc = Accuracy.mapAccuracy(outDir+"."+bestRun+".pair", truePath, totalIndiv, numIndiv, pathToKinship);
			
			//write header for output path
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhdAdj - bestLkhdAdj));
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					mapWriter.write(String.format("%d\t%d\t%.3f\n", i, j, mapAcc[i][j]));
				}
			}

			//flush
			mapWriter.flush();
			
			
			//kinship distance
			mapAcc = Accuracy.kinshipDist(outDir+"."+bestRun+".pair", truePath, totalIndiv, numIndiv, pathToKinship);
			
			
			//write header for output path
			distWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhdAdj - bestLkhdAdj));

			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					distWriter.write(String.format("%d\t%d\t%.8f\n", i, j, mapAcc[i][j]));		
				}
			}

			//flush
			distWriter.flush();
			
			
			//num correct
			int nCorrect = Accuracy.numCorrect(outDir+"."+bestRun+".pair", truePath, totalIndiv, numIndiv, pathToKinship);
			nWriter.write(String.format("%d\n", nCorrect));
			nWriter.flush();
			
			
			//TODO
			String line;
			BufferedReader pairFile = DataParser.openReader(outDir+"."+bestRun+".pair");
			while((line=pairFile.readLine())!=null){
				pairWriter.write(line+"\n");
			}
			pairWriter.flush();
			
			
			

		}

		mapWriter.close();
		distWriter.close();
		nWriter.close();
		
		

		
		
	}
		
	
	

}
	
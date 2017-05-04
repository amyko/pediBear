package test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;
import java.util.Random;

import Unused.CousinToGreatUncle;
import Unused.CousinToHalfUncle;
import Unused.CutOneLinkTwo;
import Unused.CutShiftLink;
import Unused.CutTwoLinkOne;
import Unused.FUtoHS;
import Unused.GreatUncleToCousin;
import Unused.HStoFU;
import Unused.HalfCousinToHalfGreatUncle;
import Unused.HalfGreatUncleToHalfCousin;
import Unused.HalfUncleToCousin;
import Unused.Split2;
import statistic.Accuracy;
import utility.DataParser;
import dataStructures.Path;
import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import mcmc.SimulatedAnnealing;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.FStoPO;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.POtoFS;
import mcmcMoves.Split;
import mcmcMoves.SplitLink;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.Contract;
import mcmcMoves.Stretch;



public class RunSA {

	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int maxDepth = 4;
		int sampleDepth = maxDepth;
		int numIndiv = 18;
		int totalIndiv = 18;
		double seqError = 0.01;
		double r = 1.3e-8;
		int back = 30000;
		int maxNumNodes = 200;
		double prior = numIndiv;
		int stopThresh = 5;
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmega = dir + "pathToOmega.txt";


		//SA parameters
		double[] heat = new double[800]; //200
		heat[0] = .01; //.1
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*1.01;
		System.out.println(heat[heat.length-1]);
		int coolingTime = 100000;
 		int runLength = 1;
 		int numRun = 1;
		Random rGen = new Random(1942083275L);
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .1), new Split("split", .02), new Split2("split2", 0.02), new SwapUp("swapUp", 0.02), new SwapDown("swapDown", 0.02), new SwitchSex("switchSex", 0.02), 
				new CutLink("cutLink", 0.2), new SplitLink("splitLink", 0.04), new ShiftClusterLevel("shiftClusterLevel", .02), new CutOneLinkTwo("cutOneLinkTwo", 0.15), new CutTwoLinkOne("cutTwoLinkOne", 0.02),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.02), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.02), new FStoPO("FStoPO", 0.02), new POtoFS("POtoFS",0.02), 
				new HalfUncleToCousin("halfUncleToCousin", 0.02), new CousinToHalfUncle("cousinToHalfUncle", 0.02), new CousinToGreatUncle("cousinToGreatUncle", 0.02), new GreatUncleToCousin("greatUncleToCousin", 0.02),
				new SwapDescAnc("swapDescAnc", 0.04), new Contract("contract", 0.02), new Stretch("stretch", 0.02), new HStoFU("halfSibstoFullUncle", 0.02), new FUtoHS("fullUncleToHalfSibs", 0.02),
				new ShiftClusterLevel("shiftClusterLevel", 0.04)};
		String testName = "sim4";
		String truePath = dir + "results/sim4.true";
		String accPath = dir + "results/sim4.4gen.n";
		
		double mySum = 0d;
		for(Move mov : moves) mySum += mov.getProb();
		System.out.println(mySum);
		

		//cooling schedule
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = coolingTime;
		}
		
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		
		 
		//true path
		Path[][] trueRel = Accuracy.getTruePath(truePath, totalIndiv);
		
		//open outfiles
		PrintWriter mapWriter = DataParser.openWriter(accPath+".mapAcc");
		PrintWriter distWriter = DataParser.openWriter(accPath+".kinshipDist");
			
			
		for(int t=0; t<100; t++){

			System.out.println(t);       
			
			String fileName = dir + "simPed4/"+testName+"."+t;
			String outDir = dir + "results/mcmc";
			
			double bestLkhd = Double.NEGATIVE_INFINITY;
			int bestRun = 0;
			
			for(int run=0; run<numRun; run++){
				
				//System.out.println(run);
				
				//initialize pedigree
				Pedigree ped = new Pedigree(fileName, core, maxDepth, sampleDepth, rGen, maxNumNodes, prior, numIndiv);

			
				//initialize SA
				SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, outDir+"."+run, stopThresh);

				
				//run MCMC
				double startTime = System.nanoTime();
				sa.run();
				double endTime = System.nanoTime();

				double duration = (endTime - startTime)/1e9; 

				System.out.println(String.format("Number of singletons: %d", ped.nSingletons[ped.curr]));
				System.out.println(String.format("Running time: %.1f seconds", duration));
				
				
				//if currlkhd is better than best, save mcmc output
				if(ped.getLogLikelihood() > bestLkhd){
					bestLkhd = ped.getLogLikelihood();
					bestRun = run;
				}
				
				//sanity check
				//System.out.println(ped.likelihoodAllPedigrees());
				
				
			}
			
			System.out.println(String.format("Best MCMC lkhd: %f", bestLkhd));

			
			//////////////////////////////////
			//likelihood for true pedigree

			//copy relationship
			Pedigree ped = new Pedigree(fileName, core, maxDepth, sampleDepth, rGen, maxNumNodes, prior, numIndiv);
			Path[][] mcmcmcRel = ped.getRelationships();
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = trueRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			ped.nSingletons[ped.curr] = 10;
			
			double trueLkhd = ped.likelihoodAllPedigrees();		
			System.out.println(String.format("lkhd of true pedigree: %.2f", trueLkhd));
			
			
			
			////////////////////////////////////////
			//kinship accuracy
			double[][] mapAcc = Accuracy.mapAccuracy(outDir+"."+bestRun+".pair", truePath, totalIndiv, numIndiv, pathToKinship);
			
			//write header for output path
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - bestLkhd));
			
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
			distWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - bestLkhd));
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					distWriter.write(String.format("%d\t%d\t%.8f\n", i, j, mapAcc[i][j]));
				}
			}

			//flush
			distWriter.flush();
			
			
			
			
		}
		

		
		mapWriter.close();
		distWriter.close();
		
		
		
	}
	

}
	
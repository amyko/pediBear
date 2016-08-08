package test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import statistic.Accuracy;
import statistic.Convergence;
import utility.DataParser;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStream2;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import mcmc.SimulatedAnnealing;
import mcmcMoves.CousinToGreatUncle;
import mcmcMoves.CousinToHalfUncle;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.FStoPO;
import mcmcMoves.HalfCousinToHalfGreatUncle;
import mcmcMoves.HalfGreatUncleToHalfCousin;
import mcmcMoves.HalfUncleToCousin;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.POtoFS;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.GreatUncleToCousin;



public class RunSA {
	
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int depth = 5;
		int maxDepthForSamples = depth;
		int numIndiv = 20;
		int totalIndiv = 20;
		double seqError = 0.01;
		double r = 1.3e-8;
		int back = 30000;
		int maxNumNodes = 200;
		int genTime = 16;
		double prior = .99*numIndiv;
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, r, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmega = dir + "pathToOmega.txt";

		//String marginalPath = dir + ".marginal";
		//String lkhdPath = dir + ".pairwise";
		
		//SA parameters
		double[] heat = new double[75];
		heat[0] = .3;
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*1.05;
		int coolingTime = 20000;
 		int runLength = 1;
 		int numRun = 5;
		Random rGen = new Random(19420838275L);
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .2), new Split("split", .02), new Split2("split2", 0.02), new SwapUp("swapUp", 0.02), new SwapDown("swapDown", 0.02), new SwitchSex("switchSex", 0.02), 
				new CutLink("cutLink", 0.15), new SplitLink("splitLink", 0.15), new ShiftClusterLevel("shiftClusterLevel", .02), new CutOneLinkTwo("cutOneLinkTwo", 0.15), new CutTwoLinkOne("cutTwoLinkOne", 0.02),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.02), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.02), new FStoPO("FStoPO", 0.02), new POtoFS("POtoFS",0.02), 
				new HalfUncleToCousin("halfUncleToCousin", 0.02), new CousinToHalfUncle("cousinToHalfUncle", 0.02), new CousinToGreatUncle("cousinToGreatUncle", 0.02), new GreatUncleToCousin("greatUncleToCousin", 0.02)};
		String testName = "test12";
		String outPath = dir + "results/mcmc.sample";
		String truePath = dir + "results/" +testName + ".true";
		String mapAccPath = dir + "results/"+testName+".sa.map.acc.prior.10k";
		//String mapAccPath = dir + "results/testing.sa.map.acc";

		//cooling schedule
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = coolingTime;
		}
		
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		
		
		
		//open accuracy path
		PrintWriter mapWriter = DataParser.openWriter(mapAccPath);
		 
		//true path
		Path[][] trueRel = Accuracy.getTruePath(truePath, totalIndiv);
		
			
			
		for(int t=0; t<100; t++){

			System.out.println(t);       
			
			String marginalPath = dir + "genotypes/test12.pruned.10k."+t+".marginal";
			String lkhdPath = dir + "genotypes/test12.pruned.10k."+t+".pairwise";
			
			double bestLkhd = Double.NEGATIVE_INFINITY;
			int bestRun = 0;
			
			for(int run=0; run<numRun; run++){
				
				//System.out.println(run);
				
				//initialize pedigree
				Node[] inds = new Node[numIndiv];
				for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
					inds[i] = new Node(true, i, i%2, 0, -1);
				}
				Pedigree ped = new Pedigree(depth, maxDepthForSamples, inds, core, marginalPath, lkhdPath, rGen, maxNumNodes, genTime, prior);

			
				//initialize SA
				SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, String.format("%s.%d", outPath, run), DataParser.openWriter("dummy.txt"));

				
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
				
				
			}
			
			
			
			
			
			/*
			//likelihood of best pairwise inference
			Path[][] pairwiseRel = Accuracy.pairwiseBestPed(lkhdPath, totalIndiv, numIndiv);
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = pairwiseRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			inds2.clear();
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds2.add(new Node(true, i, i%2, 0, -1));
			}
			
			
			double pairwiseLkhd = ped.likelihoodLocalPedigree(inds2);
			System.out.println();
			System.out.println(String.format("lkhd of pairwise pedigree: %.2f", pairwiseLkhd));
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					Path bestPath = pairwiseRel[i][j];
					System.out.println(String.format("(%d, %d) : (%d, %d, %d)\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
				}
			}
			*/
			
			
			
			//////////////////////////////////
			//likelihood for true pedigree
			Node[] inds = new Node[numIndiv];
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds[i] = new Node(true, i, i%2, 0, -1);
			}
			
			//copy relationship
			Pedigree ped = new Pedigree(depth, maxDepthForSamples, inds, core, marginalPath, lkhdPath, rGen, maxNumNodes, genTime, prior);
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
			double[][] mapAcc = Accuracy.mapAccuracy(String.format("%s.%d", outPath, bestRun), truePath, totalIndiv, numIndiv, pathToKinship);
			
			//write header for output path
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - bestLkhd));
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					mapWriter.write(String.format("%d\t%d\t%.3f\n", i, j, mapAcc[i][j]));
				}
			}

			//flush
			mapWriter.flush();
			
			
			
			
		}
		
		//Convergence.distanceFromTruth(outPath+"."+run, outPath+".dist."+run, truePath, numIndiv, totalIndiv, pathToKinship);
		
		//}
		

		mapWriter.close();
		
		
		/*
		//test convergence
		double[] lkhds = Convergence.getTwoHighestLikelihoods(outPath);
		Path[][] targetPed1 = Convergence.getTargetPed(outPath, lkhds[0], numIndiv);
		Path[][] targetPed2 = Convergence.getTargetPed(outPath, lkhds[1], numIndiv);
		
		double lnPrior1 = Convergence.computePrior(outPath+".2", lkhds[0], numIndiv, depth);
		double lnPrior2 = Convergence.computePrior(outPath+".2", lkhds[1], numIndiv, depth);
		System.out.println(Math.exp(lkhds[0]-lkhds[1] + lnPrior1 - lnPrior2));
		
		
		Convergence.getSampleProportion(outPath, outPath+".conv", targetPed1, targetPed2, sampleRate, numIndiv);
		*/
		
		
		//distance from truth
		//Convergence.distanceFromTruth(outPath+".t4", outPath+".dist.t4", truePath, numIndiv, totalIndiv, pathToKinship);

		
	}
	

}
	
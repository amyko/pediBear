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
		PairwiseLikelihoodCoreStream2 core = new PairwiseLikelihoodCoreStream2(seqError, r, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmega = dir + "pathToOmega.txt";

		//String marginalPath = dir + ".marginal";
		//String lkhdPath = dir + ".pairwise";
		
		//SA parameters
		double[] heat = new double[]{.2,.4,.6,.8, 1, 2, 4, 6, 8, 10};
		int coolingTime = 100000;
 		int runLength = 100;
		Random rGen = new Random(19420838275L);
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .05), new Split("split", .05), new Split2("split2", 0.05), new SwapUp("swapUp", 0.05), new SwapDown("swapDown", 0.05), new SwitchSex("switchSex", 0.05), 
				new CutLink("cutLink", 0.05), new SplitLink("splitLink", 0.05), new ShiftClusterLevel("shiftClusterLevel", .05), new CutOneLinkTwo("cutOneLinkTwo", 0.05), new CutTwoLinkOne("cutTwoLinkOne", 0.05),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.05), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.05), new FStoPO("FStoPO", 0.05), new POtoFS("POtoFS",0.05), 
				new HalfUncleToCousin("halfUncleToCousin", 0.05), new CousinToHalfUncle("cousinToHalfUncle", 0.05), new CousinToGreatUncle("cousinToGreatUncle", 0.05), new GreatUncleToCousin("greatUncleToCousin", 0.05)};
		String testName = "test11";
		String outPath = dir + "results/mcmc.sample";
		String truePath = dir + "results/" +testName + ".true";
		String mapAccPath = dir + "results/"+testName+".sa.map.acc";
		//String mapAccPath = dir + "results/testing.mcmc.map.acc";

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
			
			
			String marginalPath = dir + "pairwiseLikelihood/"+testName+".marginal."+t;
			String lkhdPath = dir + "pairwiseLikelihood/"+testName+".pairwise."+t;
			
			//initialize pedigree
			Node[] inds = new Node[numIndiv];
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds[i] = new Node(true, i, i%2, 0, -1);
			}
			Pedigree ped = new Pedigree(depth, maxDepthForSamples, inds, core, marginalPath, lkhdPath, rGen, maxNumNodes, genTime);

		
			//initialize SA
			SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, outPath);

			
			//run MCMC
			double startTime = System.nanoTime();
			sa.run();
			double endTime = System.nanoTime();

			double duration = (endTime - startTime)/1e9; 

			System.out.println(String.format("Running time: %.1f seconds", duration));

			
			//likelihood for true pedigree
			Path[][] mcmcmcRel = ped.getRelationships();
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = trueRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			List<Node> inds2 = new ArrayList<Node>();
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds2.add(new Node(true, i, i%2, 0, -1));
			}
			
			
			double trueLkhd = ped.likelihoodLocalPedigree(inds2);
			
			System.out.println(String.format("lkhd of true pedigree: %.2f", trueLkhd));
			
			
			
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
			
		
			
			//kinship accuracy
			double[][] mapAcc = Accuracy.mapAccuracy(outPath, truePath, totalIndiv, numIndiv, pathToKinship);
			
			//write header for output path
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - sa.bestLkhd));
			
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
	
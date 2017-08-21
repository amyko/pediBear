package test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import Unused.CousinToGreatUncle;
import Unused.CousinToHalfUncle;
import Unused.CutOneLinkTwo;
import Unused.CutTwoLinkOne;
import Unused.GreatUncleToCousin;
import Unused.HalfCousinToHalfGreatUncle;
import Unused.HalfGreatUncleToHalfCousin;
import Unused.HalfUncleToCousin;
import Unused.Split2;
import statistic.Accuracy;
import statistic.Convergence;
import utility.DataParser;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStream2;
import mcmc.MCMCMC;
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



public class RunMCMCMC {
	
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int depth = 4;
		int maxDepthForSamples = depth;
		int numIndiv = 20;
		int totalIndiv = 20;
		double seqError = 0.01;
		double r = 1.3e-8;
		int back = 30000;
		int maxNumNodes = 200;
		int genTime = 16;
		PairwiseLikelihoodCoreStream2 core = new PairwiseLikelihoodCoreStream2(seqError, r, back, numIndiv);
		String dir = "/Users/amy/eclipse-workspace/mcmc/simulations/";
		String pathToOmega = dir + "pathToOmega.txt";

		//String marginalPath = dir + ".marginal";
		//String lkhdPath = dir + ".pairwise";
		
		//MCMC parameters
		int nChain = 4;
		int burnIn = 1000000;
		int runLength = 1000000;
		int sampleRate = 50;
		double deltaT = .5;
		int swapInterval = 1;
		Random rGen = new Random(19420838275L);
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .05), new Split("split", .05), new Split2("split2", 0.05), new SwapUp("swapUp", 0.05), new SwapDown("swapDown", 0.05), new SwitchSex("switchSex", 0.05), 
				new CutLink("cutLink", 0.05), new SplitLink("splitLink", 0.05), new ShiftClusterLevel("shiftClusterLevel", .05), new CutOneLinkTwo("cutOneLinkTwo", 0.05), new CutTwoLinkOne("cutTwoLinkOne", 0.05),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.05), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.05), new FStoPO("FStoPO", 0.05), new POtoFS("POtoFS",0.05), 
				new HalfUncleToCousin("halfUncleToCousin", 0.05), new CousinToHalfUncle("cousinToHalfUncle", 0.05), new CousinToGreatUncle("cousinToGreatUncle", 0.05), new GreatUncleToCousin("greatUncleToCousin", 0.05)};
		String testName = "test";
		String outPath = dir + "results/mcmc.sample";
		String truePath = dir + "results/" +testName + ".true";
		String meanAccPath = dir + "results/"+testName+".mcmc.mean.acc";
		String mapAccPath = dir + "results/"+testName+".mcmc.map.acc";
		//String meanAccPath = dir + "results/testing.mcmc.mean.acc";
		//String mapAccPath = dir + "results/testing.mcmc.map.acc";


		
		
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		
		
		
		//open accuracy path
		PrintWriter meanWriter = DataParser.openWriter(meanAccPath);
		PrintWriter mapWriter = DataParser.openWriter(mapAccPath);
		 
		//true path
		Path[][] trueRel = Accuracy.getTruePath(truePath, totalIndiv);
		

			
		for(int t=0; t<1; t++){

			System.out.println(t);         
			
			
			String marginalPath = dir + "pairwiseLikelihood/"+testName+".marginal."+t;
			String lkhdPath = dir + "pairwiseLikelihood/"+testName+".pairwise."+t;
			
			//initialize chains
			List<Pedigree> chains = new ArrayList<Pedigree>(nChain);
			for(int chain=0; chain < nChain; chain++){
				
				//data
				Node[] inds = new Node[numIndiv];
				for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
					inds[i] = new Node(true, i, i%2, 0, -1);
				}
				
				Random rGen2 = new Random(2304985 + chain);
				chains.add(new Pedigree(depth, maxDepthForSamples, inds, core, marginalPath, lkhdPath, rGen2, maxNumNodes, genTime));
			}
			
			
			//initialize MCMC
			MCMCMC mcmcmc = new MCMCMC(chains, deltaT, moves, burnIn, runLength, sampleRate, swapInterval, rGen, outPath);

			
			//run MCMC
			double startTime = System.nanoTime();
			mcmcmc.run();
			double endTime = System.nanoTime();

			double duration = (endTime - startTime)/1e9; 
			
			System.out.println(String.format("cold chain index: %d", mcmcmc.coldChain));
			System.out.println(String.format("final swap rate: %.2f", mcmcmc.nSwapSuccess/((double)(burnIn+runLength)/swapInterval)));
			System.out.println(String.format("Running time: %.1f seconds", duration));

			
			//likelihood for true pedigree
			Path[][] mcmcmcRel = chains.get(mcmcmc.coldChain).getRelationships();
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = trueRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			List<Node> inds = new ArrayList<Node>();
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds.add(new Node(true, i, i%2, 0, -1));
			}
			
			
			double trueLkhd = chains.get(mcmcmc.coldChain).likelihoodLocalPedigree(inds);
			
			System.out.println(String.format("lkhd of true pedigree: %.2f", trueLkhd));
			
			
			
			//likelihood of best pairwise inference
			Path[][] pairwiseRel = Accuracy.pairwiseBestPed(lkhdPath, totalIndiv, numIndiv);
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = pairwiseRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			inds.clear();
			for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
				inds.add(new Node(true, i, i%2, 0, -1));
			}
			
			
			double pairwiseLkhd = chains.get(mcmcmc.coldChain).likelihoodLocalPedigree(inds);
			System.out.println();
			System.out.println(String.format("lkhd of pairwise pedigree: %.2f", pairwiseLkhd));
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					Path bestPath = pairwiseRel[i][j];
					System.out.println(String.format("(%d, %d) : (%d, %d, %d)\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
				}
			}
			
			
			
			
			
			//kinship accuracy
			double[][] meanAcc = Accuracy.kinshipAccuracy(outPath, truePath, totalIndiv, pathToKinship);
			double[][] mapAcc = Accuracy.mapAccuracy(outPath, truePath, totalIndiv, numIndiv, pathToKinship);
			
			
			//write header for output path
			meanWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - mcmcmc.bestLkhd));
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - mcmcmc.bestLkhd));
			//ibdWriter.write(String.format(">\t%d\t%.2f\n", t, trueLkhd - mcmcmc.bestLkhd));
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					meanWriter.write(String.format("%d\t%d\t%.3f\n", i, j, meanAcc[i][j]));
					mapWriter.write(String.format("%d\t%d\t%.3f\n", i, j, mapAcc[i][j]));
					//ibdWriter.write(String.format("%d\t%d\t%.5f\n", i, j, ibdAcc[i][j]));
				}
			}

			//flush
			meanWriter.flush();
			mapWriter.flush();
			
			
			
			
		}
		
		//Convergence.distanceFromTruth(outPath+"."+run, outPath+".dist."+run, truePath, numIndiv, totalIndiv, pathToKinship);
		
		
		
		meanWriter.close();
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
		
		
		
		/*
		///////TESTING
		//counts of each relationship
		List<Path> relationships = new ArrayList<Path>();

		relationships.add(new Path(0,0,0)); //unrelated

		// depth = 1 relationship
		//relationships.add(new Path(1,0,1)); //parent
		//relationships.add(new Path(1,1,1)); //half-sib
		//relationships.add(new Path(1,1,2)); //full-sib
		
		
		//depth = 2 relationships
		//relationships.add(new Path(2,0,1)); //grand parents
		//relationships.add(new Path(2,1,1)); //half uncle
		//relationships.add(new Path(2,2,1)); //half cousins
		//relationships.add(new Path(2,1,2)); //uncle
		//relationships.add(new Path(2,2,2)); //first cousins
		
		//depth = 3 relationships
		//relationships.add(new Path(3,0,1)); 
		//relationships.add(new Path(3,1,1)); 
		//relationships.add(new Path(3,2,1)); 
		//relationships.add(new Path(3,3,1));
		//relationships.add(new Path(3,1,2));
		//relationships.add(new Path(3,2,2)); 
		relationships.add(new Path(3,3,2)); 
		
		

		for(Path rel : relationships){
			int count = Accuracy.numSamplesWithGivenRel(outPath, rel);
			//int count = Accuracy.numSamplesWithGivenRel(outPath, 0,1,rel, false);
			System.out.println(String.format("(%d,%d,%d) %d", rel.getUp(), rel.getDown(), rel.getNumVisit(), count));
		}
		*/
		
		
		
		
		
		
		

		
	}
	

}
	
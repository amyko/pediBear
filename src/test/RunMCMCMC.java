package test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;

import statistic.Accuracy;
import utility.DataParser;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStream2;
import mcmc.MCMCMC;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.Swap;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;



public class RunMCMCMC {
	
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int depth = 4;
		int maxDepthForSamples = 4;
		int numIndiv = 6;
		double seqError = 0.01;
		double r = 1.3e-8;
		int back = 30000;
		int maxNumNodes = 10;
		int genTime = 16;
		PairwiseLikelihoodCoreStream2 core = new PairwiseLikelihoodCoreStream2(seqError, r, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmegaPath = dir + "pathToOmega.txt";

		//String marginalPath = dir + ".marginal";
		//String lkhdPath = dir + ".pairwise";
		
		//MCMC parameters
		int nChain = 6;
		int burnIn = 2000000;
		int runLength = 100000;
		int sampleRate = 100;
		double deltaT = .5;
		int swapInterval = 1;
		Random rGen = new Random(1068580L);
		Move[] moves = new Move[]{new Link("link", .1), new Cut("cut", .1), new Split("split", .05), new Split2("split2", 0.05), new Swap("swap", 0.05), new SwitchSex("switchSex", 0.05), 
				new CutLink("cutLink", .2), new SplitLink("splitLink", .2), new ShiftClusterLevel("shiftClusterLevel", .05), new CutOneLinkTwo("cutOneLinkTwo", .1), new CutTwoLinkOne("cutTwoLinkOne", .05)};
		String testName = "test6";
		String outPath = dir + "results/test.out";
		String truePath = dir + "results/" +testName + ".true";
		String relAccPath = dir + "results/"+testName+".rel.acc";
		String kinshipAccPath = dir + "results/"+testName+".kinship.acc";
		//String relAccPath = dir + "results/testing.rel.acc";
		//String kinshipAccPath = dir + "results/testing.kinship.acc";
		
		
		//open accuracy path
		PrintWriter writer1 = DataParser.openWriter(kinshipAccPath);
		PrintWriter writer2 = DataParser.openWriter(relAccPath);
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmegaPath);
		
		//true path
		truePath = dir + "results/test3.true";
		Path[][] trueRel = Accuracy.getTruePath(truePath, numIndiv);
		
		
		for(int t=0; t<100; t++){

			System.out.println(t);
			
			//write header for output path
			writer1.write(String.format(">\t%d\n", t));
			writer2.write(String.format(">\t%d\n", t));
			
			String marginalPath = dir + "pairwiseLikelihood/"+testName+".marginal."+t;
			String lkhdPath = dir + "pairwiseLikelihood/"+testName+".pairwise."+t;
			
			//initialize chains
			List<Pedigree> chains = new ArrayList<Pedigree>(nChain);
			for(int chain=0; chain < nChain; chain++){
				
				//data
				Node[] inds = new Node[numIndiv];
				for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
					inds[i] = new Node(true, i, i%2, 0, 30);
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

			System.out.println(String.format("final swap rate: %.2f", mcmcmc.nSwapSuccess/((double)(burnIn+runLength)/swapInterval)));
			System.out.println(String.format("Running time: %.1f seconds", duration));
			
			
			//Results
			System.out.println();
			for(int j=0; j<nChain; j++){
				//System.out.println(chains.get(j).getNActiveNodes());
				System.out.println(chains.get(j).getLogLikelihood());
				chains.get(j).printAdjMat();
				System.out.println();

			}
			
			
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
				inds.add(new Node(true, i, i%2, 0, 30));
			}
			
			
			double trueLkhd = chains.get(mcmcmc.coldChain).likelihoodLocalPedigree(inds);
			
			System.out.println(String.format("lkhd of true pedigree: %.2f", trueLkhd));
			

			System.out.println(String.format("cold chain index: %d", mcmcmc.coldChain));
			System.out.println(String.format("final swap rate: %.2f", mcmcmc.nSwapSuccess/((double)(burnIn+runLength)/swapInterval)));
			System.out.println(String.format("Running time: %.1f seconds", duration));

			
			
			//Accuracy
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					Accuracy.mostLikelyPath(outPath, i, j);
					
				}
			}
			
			
			
			
			
			
			double[][] kinshipAcc = Accuracy.kinshipAccuracy(outPath, truePath, numIndiv, pathToKinship);
			double[][] relAcc = Accuracy.relAccuracy(outPath, truePath, numIndiv);
			
			//System.out.println("Accuracy:");
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					//System.out.print(String.format("%.2f\t", kinshipAcc[i][j]));
					writer1.write(String.format("%d\t%d\t%.3f\n", i, j, kinshipAcc[i][j]));
				}
				//System.out.println();
			}
			
			//System.out.println("Kinship accuracy:");
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					//System.out.print(String.format("%.1e\t", relAcc[i][j]));
					writer2.write(String.format("%d\t%d\t%.3f\n", i, j, relAcc[i][j]));
				}
				//System.out.println();
			}
			
			
			//flush
			writer1.flush();
			writer2.flush();
			
			
			
		}
		
		
		writer1.close();
		writer2.close();
		
		

		
		
		
		
		/*
		///////TESTING
		//counts of each relationship
		List<Path> relationships = new ArrayList<Path>();

		relationships.add(new Path(0,0,0)); //unrelated

		// depth = 1 relationship
		relationships.add(new Path(1,0,1)); //parent
		relationships.add(new Path(1,1,1)); //half-sib
		relationships.add(new Path(1,1,2)); //full-sib
		
		
		//depth = 2 relationships
		relationships.add(new Path(2,0,1)); //grand parents
		relationships.add(new Path(2,1,1)); //half uncle
		relationships.add(new Path(2,2,1)); //half cousins
		relationships.add(new Path(2,1,2)); //uncle
		relationships.add(new Path(2,2,2)); //first cousins
		
		//depth = 3 relationships
		relationships.add(new Path(3,0,1)); 
		relationships.add(new Path(3,1,1)); 
		relationships.add(new Path(3,2,1)); 
		relationships.add(new Path(3,3,1));
		relationships.add(new Path(3,1,2));
		relationships.add(new Path(3,2,2)); 
		relationships.add(new Path(3,3,2)); 
		
		
		for(Path rel : relationships){
			int count = Accuracy.numSamplesWithGivenRel(outPath, 0,1,rel, false);
			System.out.println(String.format("(%d,%d,%d) %d", rel.getUp(), rel.getDown(), rel.getNumVisit(), count));
		}
		*/
		
		
		
		

		
	}
	

}
	
package Unused;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import dataStructures.Node;
import dataStructures.Pedigree;
import statistic.InferredRelationships;
import statistic.Statistic;
import mcmcMoves.Cut;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.Split;
import mcmcMoves.SwitchSex;

public class TestMCMC {
	
	public static void main(String[] args) throws IOException{
		
		//MCMC parameters
		int burnIn = 1000;
		int runLength = 0;
		int sampleRate = 10;
		int genTime = 16;
		Random rGen = new Random(192580);
		Move[] moves = new Move[]{new Link(.2, rGen), new Cut(.2, rGen), new Split(.2, rGen), new Split2(.2, rGen), new Swap(.2, rGen), new SwitchSex(0, rGen)};
		Statistic[] stats = new Statistic[]{new InferredRelationships()};
		
		//pedigree parameters
		int depth = 1;
		int numIndiv = 3;
		double seqError = 0.01;
		double r = 1.3e-8;
		int back = 100;
		int maxNumNodes = 50;
		PairwiseLikelihoodCoreStream core = new PairwiseLikelihoodCoreStream(seqError, r, back, numIndiv);
		String dir = System.getProperty("user.dir") + "/data/simulations/";
		String marginalPath = dir + "sim.2fam.marginal";
		String lkhdPath = dir + "sim.2fam.pairwise";
		
		//data (sampled, index, sex, depth, age)
		Node[] inds = new Node[numIndiv];
		inds[0] = new Node(true, 0, 0, 0, 60);
		inds[1] = new Node(true, 1, 1, 0, 60);
		inds[2] = new Node(true, 2, 0, 0, 20);
		//inds[3] = new Node(true, 3, 0, 1, 20);
		//inds[4] = new Node(true, 4, 0, 0, 20);
		//inds[5] = new Node(true, 5, 0, 0, 50);
		//inds[6] = new Node(true, 6, 1, 0, 50);
		

		//initialize pedigree
		Pedigree pedigree = new Pedigree(depth, inds, core, marginalPath, lkhdPath, rGen, maxNumNodes, genTime);
		
		//initialize MCMC
		MCMC mcmc = new MCMC(pedigree, burnIn, runLength, sampleRate, moves, stats);

		/*
		//TESTING FULL SIB SWAP
		Node mom = pedigree.makeNewNode(1, 1);
		Node dad = pedigree.makeNewNode(1, 0);
		pedigree.connect(mom, inds[0]);
		pedigree.connect(mom, inds[1]);
		pedigree.connect(dad, inds[0]);
		pedigree.connect(dad, inds[1]);
		pedigree.relationships[0][1].updatePath(1, 1, 2);
		pedigree.relationships[1][0].updatePath(1, 1, 2);
		pedigree.logLikelihood = pedigree.likelihoodLocalPedigree(new ArrayList<Node>(Arrays.asList(inds)));
		
		mcmc.moves[4].mcmcMove(pedigree);
		*/
		
		long startTime = System.nanoTime();
		mcmc.run();
		long endTime = System.nanoTime();

		long duration = (endTime - startTime); 
		mcmc.moves[1].mcmcMove(pedigree, 1);
			
		System.out.println(pedigree.getNActiveNodes());
		pedigree.printAdjMat();
		System.out.println(pedigree.getLogLikelihood());
		System.out.println(pedigree.likelihoodLocalPedigree(new ArrayList<Node>(Arrays.asList(inds))));

		System.out.println(duration/1e9); //in seconds
		

		
	}
	

}

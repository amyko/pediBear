package test;

import java.io.IOException;
import java.util.Random;

import dataStructures.Pedigree;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import mcmc.SimulatedAnnealing;
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
import mcmcMoves.NephewToUncle;
import mcmcMoves.OPtoPO;
import mcmcMoves.POtoFS;
import mcmcMoves.POtoHS;
import mcmcMoves.POtoOP;
import mcmcMoves.SelftoFS;
import mcmcMoves.Split;
import mcmcMoves.SplitLink;
import mcmcMoves.Stretch;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.UncletoNephew;



public class RunInuits {
	
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int maxDepth = 4;
		int sampleDepth = 2;
		int numIndiv = 100;
		double seqError = 0.01;
		int back = 30000;
		int maxNumNodes = 200;
		double prior = numIndiv;
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/inuits/";

		
		//SA parameters		
		double[] heat = new double[300]; //200
		heat[0] = .5; //.1
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*1.01;
		System.out.println(heat[heat.length-1]);
		int coolingTime = 20000;
 		int runLength = 1;
 		int numRun = 10;
 		double stopThresh = .01;
		Random rGen = new Random(194208375L);
		
		/*
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .1), new Split("split", .02), new Split2("split2", 0.02), new SwapUp("swapUp", 0.02), new SwapDown("swapDown", 0.02), new SwitchSex("switchSex", 0.02), 
				new CutLink("cutLink", 0.17), new SplitLink("splitLink", 0.07), new ShiftClusterLevel("shiftClusterLevel", .02), new CutOneLinkTwo("cutOneLinkTwo", 0.15), new CutTwoLinkOne("cutTwoLinkOne", 0.02),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.02), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.02), new FStoPO("FStoPO", 0.02), new POtoFS("POtoFS",0.02), 
				new HalfUncleToCousin("halfUncleToCousin", 0.02), new CousinToHalfUncle("cousinToHalfUncle", 0.02), new CousinToGreatUncle("cousinToGreatUncle", 0.02), new GreatUncleToCousin("greatUncleToCousin", 0.02),
				new SwapDescAnc("swapDescAnc", 0.04), new Contract("contract", 0.02), new Stretch("stretch", 0.02), new HalfSibstoFullUncle("halfSibstoFullUncle", 0.02), new FullUncletoHalfSibs("fullUncleToHalfSibs", 0.02),
				new ShiftClusterLevel("shiftClusterLevel", 0.04)};
				*/
		
		
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
		
		
		String testName = "100tasiilaq.admixed0.5.aims1e-5.prune0.05";
		
		


		//cooling schedule
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = coolingTime;
		}


		double bestLkhd = Double.NEGATIVE_INFINITY;
		int bestRun = 0;
		
		for(int run=0; run<numRun; run++){
			
			//initialize pedigree
			Pedigree ped = new Pedigree(dir+testName, core, maxDepth, sampleDepth, rGen, maxNumNodes, prior, numIndiv);

			/*
			//FOR TESTING ONLY
			Node dad = ped.getNode(0);
			Node mom = ped.getNode(1);
			
			
			Node c1 = ped.getNode(2);
			c1.addParent(dad);
			c1.addParent(mom);
			dad.addChild(c1);
			mom.addChild(c1);
			
			
			Node c2 = ped.getNode(3);
			c2.addParent(dad);
			c2.addParent(mom);
			dad.addChild(c2);
			mom.addChild(c2);
			
			
			Node c3 = ped.getNode(4);
			c3.addParent(dad);
			c3.addParent(mom);
			dad.addChild(c3);
			mom.addChild(c3);
			
			
			Node c4 = ped.getNode(5);
			c4.addParent(dad);
			c4.addParent(mom);
			dad.addChild(c4);
			mom.addChild(c4);
			
			
			ped.updateAdjMat(c1);
			ped.updateAdjMat(c2);
			ped.updateAdjMat(c3);
			ped.updateAdjMat(c4);
			System.out.println(ped.pairwiseLkhd());
			*/
		
			//initialize SA
			SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, dir+testName+".revised."+run, stopThresh);

			
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
		
		System.out.println(bestRun);
			
					
	
	}
	

}
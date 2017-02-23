package test;

import java.io.BufferedReader;
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
import mcmcMoves.Contract;
import mcmcMoves.CousinToGreatUncle;
import mcmcMoves.CousinToHalfUncle;
import mcmcMoves.Cut;
import mcmcMoves.CutLink;
import mcmcMoves.CutOneLinkTwo;
import mcmcMoves.CutTwoLinkOne;
import mcmcMoves.FStoPO;
import mcmcMoves.FUtoHS;
import mcmcMoves.HalfCousinToHalfGreatUncle;
import mcmcMoves.HalfGreatUncleToHalfCousin;
import mcmcMoves.HStoFU;
import mcmcMoves.HalfToFullSibs;
import mcmcMoves.HalfUncleToCousin;
import mcmcMoves.Link;
import mcmcMoves.Move;
import mcmcMoves.POtoFS;
import mcmcMoves.Split;
import mcmcMoves.Split2;
import mcmcMoves.SplitLink;
import mcmcMoves.Stretch;
import mcmcMoves.SwapDescAnc;
import mcmcMoves.SwapDown;
import mcmcMoves.SwapUp;
import mcmcMoves.SwitchSex;
import mcmcMoves.ShiftClusterLevel;
import mcmcMoves.GreatUncleToCousin;



public class RunJays {
	
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	public static void main(String[] args) throws IOException{
		

		//pedigree parameters
		int depth = 4;
		int numIndiv = 75;
		int totalIndiv = 75;
		double seqError = 0.01;
		double r = 1.5e-8;
		int back = 30000;
		int maxNumNodes = 500;
		int prior = numIndiv/2;
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, r, back, numIndiv);
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/";
		String pathToOmega = dir + "pathToOmega.txt";

		//SA parameters
		double[] heat = new double[800]; //200
		heat[0] = .01; //.1
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*1.01;
		System.out.println(heat[heat.length-1]);
		int coolingTime = 100000;
 		int runLength = 1;
 		int numRun = 1;


		Random rGen = new Random(1942038275L);
		Move[] moves = new Move[]{new Link("link", .05), new Cut("cut", .1), new Split("split", .02), new Split2("split2", 0.02), new SwapUp("swapUp", 0.02), new SwapDown("swapDown", 0.02), new SwitchSex("switchSex", 0.02), 
				new CutLink("cutLink", 0.17), new SplitLink("splitLink", 0.07), new ShiftClusterLevel("shiftClusterLevel", .02), new CutOneLinkTwo("cutOneLinkTwo", 0.15), new CutTwoLinkOne("cutTwoLinkOne", 0.02),
				new HalfCousinToHalfGreatUncle("halfCousinToHalfGreatUncle", 0.02), new HalfGreatUncleToHalfCousin("halfGreatUncleToHalfCousin", 0.02), new FStoPO("FStoPO", 0.02), new POtoFS("POtoFS",0.02), 
				new HalfUncleToCousin("halfUncleToCousin", 0.02), new CousinToHalfUncle("cousinToHalfUncle", 0.02), new CousinToGreatUncle("cousinToGreatUncle", 0.02), new GreatUncleToCousin("greatUncleToCousin", 0.02),
				new SwapDescAnc("swapDescAnc", 0.04), new Contract("contract", 0.02), new Stretch("stretch", 0.02), new HStoFU("halfSibstoFullUncle", 0.02), new FUtoHS("fullUncleToHalfSibs", 0.02),
				new ShiftClusterLevel("shiftClusterLevel", 0.04)};
		

		//cooling schedule
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = coolingTime;
		}
		
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmega);
		
		
		//files
		String testName = "75jays.geno0.prune10.0.1";
		String truePath = dir + "75jays.true";
		String indexPath = dir + testName + ".index2name";
		String marginalPath = dir +testName+".marginal";
		String lkhdPath = dir + testName+".pairwise";
		String convPath = dir + testName + ".conv";
		

		//true path
		Path[][] trueRel = Accuracy.getTruePath(truePath, totalIndiv);
		
		/*
		// samples
		BufferedReader reader = DataParser.openReader(indexPath);
		reader.readLine();
		Node[] inds = new Node[numIndiv];
		for(int i=0; i<numIndiv; i++){ //(sampled, index, sex, depth ,age)
			int sex = Integer.parseInt(reader.readLine().split("\t")[2]);
			inds[i] = new Node(true, i, sex, 0, -1);
		}
		reader.close();
		*/
		
			
		//convergence
		PrintWriter convWriter = DataParser.openWriter(convPath);
		
			
		for(int t=0; t<numRun; t++){

			System.out.println(t);  
			String outPath = dir + "mcmc.sample."+t;
			
			//open accuracy path
			String mapAccPath = dir + testName+".sa.mapAcc."+t;
			PrintWriter mapWriter = DataParser.openWriter(mapAccPath);

			//initialize pedigree
			Pedigree ped = new Pedigree(dir+testName, core, depth, rGen, maxNumNodes, prior, numIndiv);

		
			//initialize SA
			SimulatedAnnealing sa = new SimulatedAnnealing(ped, heat, coolingSchedule, moves, runLength, rGen, dir+testName+"."+t);

			
			
			
			//run MCMC
			double startTime = System.nanoTime();
			sa.run();
			double endTime = System.nanoTime();
			double duration = (endTime - startTime)/1e9;
			System.out.println(String.format("Running time: %.1f seconds", duration));

			
			/*
			//sanity check
			for(int k=0; k<numIndiv; k++){
				ped.updateAdjMat(ped.getNode(k));
			}
			System.out.println(String.format("MCMC likelihood should be: %f", ped.likelihoodAllPedigrees()));
			

			
			
			//likelihood for true pedigree
			Path[][] mcmcmcRel = ped.getRelationships();
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path real = trueRel[i][j];
					mcmcmcRel[i][j].updatePath(real.getUp(), real.getDown(), real.getNumVisit());
					
				}
			}
			
			double trueLkhd = ped.likelihoodAllPedigrees();		
			System.out.println(String.format("lkhd of true pedigree: %.2f", trueLkhd));
			 */
			

			
			//kinship accuracy
			double[][] mapAcc = Accuracy.mapAccuracy(dir+testName+"."+t+".pair", truePath, totalIndiv, numIndiv, pathToKinship);
			
			//write header for output path
			mapWriter.write(String.format(">\t%d\t%.2f\n", t, sa.bestLkhd));
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					mapWriter.write(String.format("%d\t%d\t%.3f\n", i, j, mapAcc[i][j]));
				}
			}

			//flush
			mapWriter.flush();
			mapWriter.close();
			
			
			
		}
		
		convWriter.close();
		
	}
	

}
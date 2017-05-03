package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import likelihood.LDStreamPed;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Relationship;
import dataStructures.Path;
import simulator.SimulatorStream;
import simulator.SimulatorStreamPed;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class SimulatePedigree {
		
		public static double recombRate = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(524140L);

		

		public static void computeInfo(String hapmapPath, String simPath, String outPath, int back) throws IOException{
			
			//concatenate
			String tempPath = dir+"hapmap.children.txt";
			DataParser.concatFiles(new String[]{hapmapPath, simPath}, tempPath, new  int[][]{{},{3,4,5,8,9,10}});
			
			//compute info
			LDStream.writeLdOutfile(tempPath,outPath,back);
		}
		
		
		public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, int t) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginal(fileName+"."+t+".tped", dir+ "unrelated/msprime.unrel.try2.pruned.info", ids);
			
			
			//open outfile
			PrintWriter writer = DataParser.openWriter(fileName+"."+t+".marginal");				
			
			//write marginals to file
			for (double item : marginals){
				writer.write(String.format("%f\n",item));
			}
			writer.close();
		}
		
		
		public static void computePairwise(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, List<Relationship> relationships, int t) throws IOException{
			
			int numIndiv = ids.length;
			int numRel = relationships.size();
			
			//likelihood
			PrintWriter writer = DataParser.openWriter(fileName+"."+t+".pairwise");
				
			double[][][] lkhd = core.forwardAlgorithm(fileName+"."+t+".tped", dir+ "unrelated/msprime.unrel.try2.pruned.info", ids, relationships);

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
		
		
		
		public static void simPed5(SimulatorStreamPed sim, String dataPath, String outPath, int t) throws IOException{
					
			
			//make 3 full sibs
			int momID = 0;
			int dadID = 1;
			sim.makeChildren(dataPath, dataPath, outPath+".3sibs", momID, dadID, rgen, 3);
			
			//first cluster
			momID = 0;
			dadID = 2;
			sim.makeChildren(outPath+".3sibs", dataPath, outPath+".indiv2", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 3;
			sim.makeChildren(outPath+".indiv2", dataPath, outPath+".indiv3", momID, dadID, rgen, 1);
			
			
			//second cluster
			momID = 2;
			dadID = 4;
			sim.makeChildren(outPath+".3sibs", dataPath, outPath+".indiv45", momID, dadID, rgen, 2);
			
			momID = 1;
			dadID = 5;
			sim.makeChildren(outPath+".indiv45", dataPath, outPath+".indiv6", momID, dadID, rgen, 1);
			
			
			int[] unrel = new int[2*13];
			for(int i=0; i<unrel.length/2; i++){
				unrel[2*i] = 2*(i+2+3);
				unrel[2*i+1] = 2*(i+2+3)+1;
			}
			
			//concatenate individuals
			String[] paths = new String[]{outPath+".3sibs", outPath+".indiv2", outPath+".indiv3", outPath+".indiv45", outPath+".indiv6", dataPath};
			int[][] cols = new int[][]{new int[]{0,1,2,3,6,7,8,9}, new int[]{4,5}, new int[]{4,5}, new int[]{4,5,6,7}, new int[]{4,5}, unrel};
			DataParser.concatFiles(paths, outPath+".noError.tped", cols);
			
			
			//add error
			sim.addError(outPath+".noError.tped", outPath+"."+t+".tped", seqError, rgen);
			
			
		}
		
		
		public static void simPed1(SimulatorStreamPed sim, String dataPath, String outPath, int t) throws IOException{
					
			
			//make 3 full sibs
			int momID = 0;
			int dadID = 1;
			sim.makeChildren(dataPath, dataPath, outPath+".3sibs", momID, dadID, rgen, 3);
			
			//first cluster
			momID = 0;
			dadID = 2;
			sim.makeChildren(outPath+".3sibs", dataPath, outPath+".2sibs", momID, dadID, rgen, 2);
			
			momID = 0;
			dadID = 3;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 4;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv1", momID, dadID, rgen, 1);
			
			momID = 1;
			dadID = 5;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 6;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv2", momID, dadID, rgen, 1);
			
			
			//second cluster
			momID = 1;
			dadID = 7;
			sim.makeChildren(outPath+".3sibs", dataPath, outPath+".2sibs", momID, dadID, rgen, 2);
			
			momID = 0;
			dadID = 8;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv34", momID, dadID, rgen, 2);
			
			momID = 0; //half sib
			dadID = 9;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0; 
			dadID = 10;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv5", momID, dadID, rgen, 1);
			
			momID = 1; 
			dadID = 11;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0; 
			dadID = 12;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv6", momID, dadID, rgen, 1);
			
			
			//third cluster
			momID = 2;
			dadID = 13;
			sim.makeChildren(outPath+".3sibs", dataPath, outPath+".3sibsAgain", momID, dadID, rgen, 3);
			
			momID = 0;
			dadID = 14;
			sim.makeChildren(outPath+".3sibsAgain", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 15;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv7", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 14;
			sim.makeChildren(outPath+".3sibsAgain", dataPath, outPath+".indiv8", momID, dadID, rgen, 1);
			
			
			momID = 1;
			dadID = 16;
			sim.makeChildren(outPath+".3sibsAgain", dataPath, outPath+".indiv9", momID, dadID, rgen, 1);
			
			momID = 2;
			dadID = 17;
			sim.makeChildren(outPath+".3sibsAgain", dataPath, outPath+".temp", momID, dadID, rgen, 1);
			
			momID = 0;
			dadID = 18;
			sim.makeChildren(outPath+".temp", dataPath, outPath+".indiv10", momID, dadID, rgen, 1);
			
			
			dadID++;
			int[] unrel = new int[20];
			for(int i=0; i<unrel.length/2; i++){
				unrel[2*i] = 2*(i+dadID+2);
				unrel[2*i+1] = 2*(i+dadID+2)+1;
			}
			
			//concatenate individuals
			String[] paths = new String[]{outPath+".indiv1", outPath+".indiv2", outPath+".indiv34", outPath+".indiv5", 
					outPath+".indiv6", outPath+".indiv7", outPath+".indiv8", outPath+".indiv9", outPath+".indiv10", dataPath};
			int[][] cols = new int[][]{new int[]{0,1,2,3,4,5}, new int[]{4,5}, new int[]{4,5,6,7}, new int[]{4,5}, new int[]{4,5},
					new int[]{4,5},new int[]{4,5},new int[]{4,5},new int[]{4,5},unrel};
			DataParser.concatFiles(paths, outPath+".noError.tped", cols);
			
			
			//add error
			sim.addError(outPath+".noError.tped", outPath+"."+t+".tped", seqError, rgen);
			
			
		}
		

		public static void outbredFirstCousins(SimulatorStreamPed sim, String dataPath, String outPath, int t) throws IOException{
					
			
			//make 2 full sibs
			int momID = 0;
			int dadID = 1;
			sim.makeChildren(dataPath, dataPath, outPath+".2sibs", momID, dadID, rgen, 2);
			
			//first cluster
			momID = 0;
			dadID = 2;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv1", momID, dadID, rgen, 1);
			
			momID = 1;
			dadID = 3;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv2", momID, dadID, rgen, 1);

			
			//concatenate individuals
			String[] paths = new String[]{outPath+".indiv1", outPath+".indiv2"};
			int[][] cols = new int[][]{new int[]{0,1,2,3,4,5}, new int[]{4,5}};
			DataParser.concatFiles(paths, outPath+".noError.tped", cols);
			
			
			//add error
			sim.addError(outPath+".noError.tped", outPath+"."+t+".tped", seqError, rgen);
			
			
		}
		
		
		public static void inbredFirstCousins(SimulatorStreamPed sim, String dataPath, String outPath, int t) throws IOException{
					
			
			//make 2 full sibs
			int momID = 0;
			int dadID = 1;
			sim.makeChildren(dataPath, dataPath, outPath+".fs", momID, dadID, rgen, 2);
			
			/*
			momID = 0;
			dadID = 2;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv1", momID, dadID, rgen, 1);
			
			momID = 1;
			dadID = 3;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv2", momID, dadID, rgen, 1);
			
			
			//inbred 
			momID = 0;
			dadID = 0;
			sim.makeChildren(outPath+".indiv1", outPath+".indiv2", outPath+".2sibs", momID, dadID, rgen, 2);
			*/
			
			momID = 0;
			dadID = 0;
			sim.makeChildren(outPath+".fs", outPath+".fs", outPath+".2sibs", momID, dadID, rgen, 2);
			
			momID = 0;
			dadID = 4;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv1", momID, dadID, rgen, 1);
			
			momID = 1;
			dadID = 5;
			sim.makeChildren(outPath+".2sibs", dataPath, outPath+".indiv2", momID, dadID, rgen, 1);
			
			
			
			
			
			//concatenate individuals
			String[] paths = new String[]{outPath+".indiv1", outPath+".indiv2"};
			int[][] cols = new int[][]{new int[]{0,1,2,3,4,5}, new int[]{4,5}};
			DataParser.concatFiles(paths, outPath+".noError.tped", cols);
			
			
			//add error
			sim.addError(outPath+".noError.tped", outPath+"."+t+".tped", seqError, rgen);
			
			
		}
		
		
		public static void parentOffspring(SimulatorStreamPed sim, String dataPath, String outPath, int t) throws IOException{
					
			
			//make 2 full sibs
			int momID = 0;
			int dadID = 1;
			sim.makeChildren(dataPath, dataPath, outPath+".2sibs", momID, dadID, rgen, 1);
		
			
			//concatenate individuals
			String[] paths = new String[]{dataPath, outPath+".2sibs"};
			int[][] cols = new int[][]{new int[]{0,1,2,3,4,5}, new int[]{4,5}};
			DataParser.concatFiles(paths, outPath+".noError.tped", cols);
			
			
			//add error
			sim.addError(outPath+".noError.tped", outPath+"."+t+".tped", seqError, rgen);
			
			
		}
		
		
		public static void main(String[] args) throws IOException{
			
			////////////////////////////
			
			//initialize simulator
			SimulatorStreamPed sim = new SimulatorStreamPed(recombRate);
			
			//relationships
			List<Relationship> relationships = new ArrayList<Relationship>();

			relationships.add(new Relationship(13d, new double[] {1d, 0d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
			// depth = 1 relationship
			relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
			relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
			
			
			relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
			
			//depth = 2 relationships
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2)})); //first cousins
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(3,1,2)}));
			
			
			//depth = 3 relationships
			relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1)})); 
			relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1)}));
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2)})); 
			relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2)})); 
			
			//depth = 4 relationships
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(4,1,2)}));
			relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1)}));
			relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1)}));
			relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2)}));
			relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2)}));
			
			//////////////////////////////
	
			
			//dir
			String simDir = dir + "simulations/simPed7/parentOffspring";
			String dataDir = dir + "unrelated/msprime.unrel.try2.pruned.tped";

			
			//individuals
			int numIndiv = 2;
			int[] indCols = new int[numIndiv];
			for(int i=0; i<numIndiv; i++) indCols[i] = i;			
			
			
			
			//pairwise core
			PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, numIndiv);
			

	
			for(int t=0; t<1; t++){
				
				System.out.println(t);
		
				
				System.out.println("Simulating");		
				parentOffspring(sim, dataDir, simDir, t);
				
				/*
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginals(core, simDir, indCols, t);
				
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");
				long startTime = System.nanoTime();		
				computePairwise(core, simDir, indCols, relationships, t);	
				System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
				*/
				
			
			}
		
			
			
			//true path
			//sim.writeTruePathFile(simDir+"test12.ped", simDir+"10indiv.ped", ids);
			
			
			

		}
		
		


}
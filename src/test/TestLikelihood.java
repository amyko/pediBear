package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import likelihood.LDStreamPed;
import likelihood.PairwiseLikelihoodCoreStream2;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Relationship;
import dataStructures.Path;
import simulator.SimulatorStream;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class TestLikelihood {
		
		public static double recombRate = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(526564L);


	
		
		public static void simulateCousins(String inPath, String outPath, int numGen, boolean full, int howMany, int nCluster, int start) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(recombRate);

			
			String[] outNames = new String[nCluster];
			int[][] outCols = new int[nCluster][howMany];
			for(int i=0; i<nCluster; i++){
				for(int j=0; j<howMany; j++)
					outCols[i][j] = j;
			}
			
			
			for(int clu=0; clu<nCluster; clu++){
			
				String cluOutPath = dir + "cluster."+clu+".out";
				outNames[clu] = cluOutPath;
				
				
				int commonAnc = start;
				start++;
				//System.out.println(start);
				
				//first generation
				if(full){
					sim.makeChildren(inPath, inPath, inPath, cluOutPath, commonAnc, start, rgen, howMany, "p\n");
					start++;
					//System.out.println(start);
				}
				else{
					
					String[] fileNames = new String[howMany];
					int[][] cols = new int[howMany][howMany];
					
					//mate with different spouse
					for(int i=0; i<howMany; i++){
						String tempPath = dir + "temp."+i+".txt";
						fileNames[i] = tempPath;
						cols[i] = new int[1];
						
						sim.makeChildren(inPath, inPath, inPath, tempPath, commonAnc, start, rgen, 1, "p\n");
						start++;
						//System.out.println(start);
					}
					
					//concatenate
					DataParser.concatFiles(fileNames, cluOutPath, cols);
					
				}
				
				
				
				//later generations
				int currGen = numGen - 1;
				String tempIn = cluOutPath;
				while(currGen > 0){
					
					//file names
					String[] fileNames = new String[howMany];
					int[][] cols = new int[howMany][howMany];
					
					//make children for each family
					for(int i=0; i<howMany; i++){
						
						String tempOut = dir + "temp."+i+".txt";
						fileNames[i] = tempOut;
						cols[i] = new int[1];
						
						sim.makeChildren(inPath, tempIn, inPath, tempOut, i, start, rgen, 1, "p\n");
						start++;	
						//System.out.println(start);
						
					}
					
					//concatenate
					DataParser.concatFiles(fileNames, cluOutPath, cols);
					
					//increment
					currGen--;
					
					
				}
			
			}
			
			
			//concatenate
			DataParser.concatFiles(outNames, outPath, outCols);

			
			
		}
		
		

		public static int simulateChildren(String inPath, String outPath, String momPath, int numGen, int howMany, int momCol, int start) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(recombRate);
			
			//first generation
			String firstPath = dir+"simulations/parent.out";
			sim.makeChildren(inPath, momPath, inPath, firstPath, momCol, start, rgen, howMany, "p\n");
			start++;
			
			//later generations
			int currGen = numGen - 1;
			while(currGen > 0){
				
				//file names
				String[] fileNames = new String[howMany];
				int[][] cols = new int[howMany][howMany];
				
				//make children for each family
				for(int i=0; i<howMany; i++){
					
					String tempOut = dir + "simulations/ctemp."+i+".txt";
					fileNames[i] = tempOut;
					cols[i] = new int[1];
					
					sim.makeChildren(inPath, firstPath, inPath, tempOut, i, start, rgen, 1, "p\n");
					start++;	
				}
				
				//concatenate
				DataParser.concatFiles(fileNames, outPath, cols);
				
				//increment
				currGen--;
				firstPath = outPath;
				
				
			}
			
			
			return start;
		
		}


		
		public static void computeInfo(String hapmapPath, String simPath, String outPath, int back) throws IOException{
			
			//concatenate
			String tempPath = dir+"hapmap.children.txt";
			DataParser.concatFiles(new String[]{hapmapPath, simPath}, tempPath, new  int[][]{{},{3,4,5,8,9,10}});
			
			//compute info
			LDStream.writeLdOutfile(tempPath,outPath,back);
		}
		
		
		public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, int t) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginal(fileName+"."+t+".tped", dir+ "simulations/simPed2/msprime.unrel.try2.10k.info", ids);
			
			
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
				
			double[][][] lkhd = core.forwardAlgorithm(fileName+"."+t+".tped", dir+ "simulations/simPed2/msprime.unrel.try2.10k.info", ids, relationships);

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
		
		
		

		
		public static void main(String[] args) throws IOException{
			
			////////////////////////////
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(recombRate);
			
			//relationships
			List<Relationship> relationships = new ArrayList<Relationship>();

			relationships.add(new Relationship(13d, new double[] {1d, 0d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
			//relationships.add(new Relationship(13d, new double[] {511d/512d, 1d/512d, 0d}, new Path[]{new Path(0,0,0)})); //unrelated
			
			// depth = 1 relationship
			relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
			relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
			
			
			relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
			
			//depth = 2 relationships
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
			
			//depth = 3 relationships
			relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1), new Path(5,0,1)})); 
			relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1), new Path(5,1,1)}));
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
			relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2), new Path(5,1,2)})); 
			
			//depth = 4 relationships
			relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1), new Path(5,2,1)}));
			relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1), new Path(5,3,1)}));
			relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2), new Path(5,2,2)}));
			relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2), new Path(5,3,2)}));
			
			
			//depth = 5 relationships
			relationships.add(new Relationship(9d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,4,1)}));
			relationships.add(new Relationship(10d, new double[] {511d/512d, 1d/512d, 0d}, new Path[]{new Path(5,5,1)}));
			relationships.add(new Relationship(11d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(5,4,2)}));
			relationships.add(new Relationship(12d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,5,2)}));
 			
			
			//////////////////////////////
	
			
			//dir
			String simDir = dir + "simulations/simPed2/";

			
			//individuals
			int numIndiv = 20;
			int[] indCols = new int[numIndiv];
			for(int i=0; i<numIndiv; i++) indCols[i] = i;			
			
			
			int chrStart = 1;
			int chrEnd = 23;
			boolean full = true;
			int nSmallCluster = 2;
			int nBigCluster = 1;
			int nChildren = 5;
			int nGen = 4;
			int howManyUnrelated = 10;
			String testName = "sim2";
			
			
			//pairwise core
			PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, recombRate, back, numIndiv);
			

			
			
			//compute info
			//System.out.println("Computing info");
			//LDStreamPed.writeLdOutfile(dir+"unrelated/msprime.unrel.try2.10k", simDir+"msprime.unrel.try2.10k.info", back);	
			
			
	
			//int[] cols = new int[]{8, 15, 37, 46, 56, 65, 83, 93, 114, 130};
			int[] cols = new int[]{23, 24, 36, 38, 56, 65, 83, 93, 114, 130};
			int[] ids = new int[]{30, 31,47,50,74,86,111,124,152,174};
			//int[] ids = new int[]{10,19,49,61,74,86,111,124,152,174};
			
			
			for(int t=0; t<50; t++){
				
				System.out.println(t);
		
				//sim.simulatePopulation(dataDir+"msprime.geno.pruned.", simDir+"sim.test.geno.", simDir+"pairwiseLikelihood/"+testName+".ped", 4, rgen, chrStart, chrEnd, cols, howManyUnrelated);
				
				/*
				// simulation
				System.out.println("Simulating pedigrees");
				for(int i=chrStart; i<chrEnd; i++){
					
					//files
					String unrelated = dataDir+"msprime.geno.pruned."+i;
					
					simulateCousins(unrelated, simDir+"sim.test.geno."+i, nGen, full, nChildren, nSmallCluster, 8);	
					
					
					//first gen
					int start = 8;
					sim.makeChildren(unrelated, unrelated, unrelated, simDir+"first.out", start, start+1, rgen, nSmallCluster, "children\n");
					start+=2;
					
					//later gen
					String[] fileNames = new String[nSmallCluster];
					int[][] cols = new int[nSmallCluster][nChildren];
					for(int j=0; j<nSmallCluster; j++){
						
						//file names
						fileNames[j] = simDir + "sim.fam.geno."+j;
						for(int k=0; k<nChildren; k++) cols[j][k] = k;
						
						
						start = simulateChildren(unrelated, fileNames[j], simDir+"first.out", nGen-1, nChildren, j, start);

					}
					
					
					
					//concatenate families
					DataParser.concatFiles(fileNames, simDir+"sim.test.geno."+i, cols);
					
										
					
					
				}
				*/
				
				/*
				//add error
				System.out.println("Adding error");
				for(int i=chrStart; i<chrEnd; i++){
					LDStream.addError(simDir+"sim.test.geno."+i, simDir+"sim.test.geno.error."+t+"."+i, seqError, rgen);	
				}
				*/
				
				
				

				
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginals(core, simDir+testName, indCols, t);
				
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");
				long startTime = System.nanoTime();		
				computePairwise(core, simDir+testName, indCols, relationships, t);	
				System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
				
				
			
			}
		
			
			
			//true path
			//sim.writeTruePathFile(simDir+"test12.ped", simDir+"10indiv.ped", ids);
			
			
			

		}
		
		


}
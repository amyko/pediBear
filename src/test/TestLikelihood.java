package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import likelihood.PairwiseLikelihoodCoreStream2;
import dataStructures.Relationship;
import dataStructures.Path;
import simulator.SimulatorStream;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class TestLikelihood {
		
		public static double r = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(526564L);


	
		
		public static void simulateCousins(String inPath, String outPath, int numGen, boolean full, int howMany, int nCluster, int start) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);

			
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
			SimulatorStream sim = new SimulatorStream(r);
			
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
		
		
		public static void computeMarginals(String dataDir, String simDir, String infoDir, String outPath, int[] indCols, int chrStart, int chrEnd) throws IOException{
			
			//compute
			PairwiseLikelihoodCoreStream2 pwStream = new PairwiseLikelihoodCoreStream2(seqError, r, back, indCols.length);
			double[] marginals = new double[indCols.length];
			for (int chr=chrStart; chr<chrEnd; chr++){
				
				String genoPath = simDir + chr;
				String infoPath = infoDir + chr;
				
				double[] temp = pwStream.computeMarginal(genoPath, infoPath, indCols);
				for (int i=0; i<marginals.length; i++) marginals[i] += temp[i];
				
			}

			//write header
			PrintWriter writer = DataParser.openWriter(outPath);		
			for (int i=0; i<indCols.length; i++){
				writer.write(i+"\t");
			}
			writer.write("\n");
			
			
			
			//write marginals
			for (double item : marginals){
				writer.write(String.format("%f\t",item));
			}
			writer.close();
		}
		
		
		public static void computePairwise(String dataDir, String simDir, String infoDir, String outPath, int[] indCols, List<Relationship> relationships, int chrStart, int chrEnd) throws IOException{
			
			int numIndiv = indCols.length;
			int numRel = relationships.size();
			
			//likelihood
			PairwiseLikelihoodCoreStream2 pwStream = new PairwiseLikelihoodCoreStream2(seqError, r, back, numIndiv);
			PrintWriter writer = DataParser.openWriter(outPath);
				
		
			double[][][] lkhd = new double[numRel][numIndiv][numIndiv];
			
			for (int chr=chrStart; chr<chrEnd; chr++){
			
				//compute likelihood
				String genoPath = simDir + chr;
				String infoPath =  infoDir +chr;
				double[][][] tempLkhd = pwStream.forwardAlgorithm(genoPath, infoPath, indCols, relationships);
				
				
				//add to total likelihood
				for(int i=0; i<numIndiv; i++){
					for(int j=i+1; j<numIndiv; j++){
						
						for(int k=0; k<numRel; k++){
							lkhd[k][i][j] += tempLkhd[k][i][j];	
						}
						
					}
				}
				
				
				//System.out.println();
				
			}
			
			
	
			//write to file
			for(int k=0; k<numRel; k++){
				
				for(Path path : relationships.get(k).getAllPaths()){
					//Path path = relationships.get(k).getPath();
					String header = String.format(">\t%d\t%d\t%d\n",path.getUp(), path.getDown(), path.getNumVisit());
					writer.write(header);
					DataParser.writeMatrixToFile(writer, lkhd[k]);
					writer.flush();					
				}

				/*
				//print lkhd
				System.out.print(header);
				for(int i=0; i<numIndiv; i++){
					for(int j=0; j<numIndiv; j++){
						System.out.format("%f\t",lkhd[k][i][j]);
					}
					System.out.println();
				}
				System.out.println();
				*/
			
			}
		
			//write to file
			writer.close();
			
		}
		
		
		

		
		public static void main(String[] args) throws IOException{
			
			////////////////////////////
			//dir
			String dataDir = dir + "unrelated/";
			String simDir = dir + "simulations/";

			
			//individuals
			int numIndiv = 6;
			int[] indCols = new int[numIndiv];
			for(int i=0; i<numIndiv; i++) indCols[i] = i;			
			
			
			int chrStart = 1;
			int chrEnd = 23;
			boolean full = true;
			int nSmallCluster = 2;
			int nBigCluster = 1;
			int nChildren = 3;
			int nGen = 4;
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);
			
			//relationships
			List<Relationship> relationships = new ArrayList<Relationship>();

			relationships.add(new Relationship(13d, new double[] {1-1d/512d, 1d/512d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
			//relationships.add(new Relationship(11d, new double[] {1d, 0, 0}, new Path[]{new Path(0,0,0)})); //unrelated
			
			// depth = 1 relationship
			relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
			relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
			relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
			
			//depth = 2 relationships
			//relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, 2,0,1)); //grand parents
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
			
			//depth = 3 relationships
			//relationships.add(new Relationship(3d, new double[] {.75, .25, 0d}, 3,0,1)); 
			//relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, 3,1,1)); 
			relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1), new Path(5,0,1)})); 
			relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1), new Path(5,1,1)}));
			//relationships.add(new Relationship(6d, new double[] {.75, .25, 0d}, 3,1,2));
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
			relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2), new Path(5,1,2)})); 
			
			//depth = 4 relationships
			//relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, 4,0,1)); 
			//relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, 4,1,1)); 
			//relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, 4,2,1)); 
			relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1), new Path(5,2,1)}));
			relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1), new Path(5,3,1)}));
			//relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, 4,1,2)); 
			//relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, 4,2,2)); 
			relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2), new Path(5,2,2)}));
			relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2), new Path(5,3,2)}));
			
			
			//depth = 5 relationships
			relationships.add(new Relationship(9d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,4,1)}));
			relationships.add(new Relationship(10d, new double[] {511d/512d, 1d/512d, 0d}, new Path[]{new Path(5,5,1)}));
			relationships.add(new Relationship(11d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(5,4,2)}));
			relationships.add(new Relationship(12d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,5,2)}));
			

			//depth = distantly related
			//relationships.add(new Relationship(11d, new double[] {1-1d/128d, 1d/128d, 0}, new Path[]{new Path(5,0,1), new Path(5,1,1), new Path(5,2,1), new Path(5,3,1), new Path(5,4,1), new Path(5,5,1), new Path(5,1,2), new Path(5,2,2), new Path(5,3,2), new Path(5,4,2), new Path(5,5,2)})); 
			
			//////////////////////////////
	
			
			//prune ld
			System.out.println("Pruning LD");
			for(int i=1; i<23; i++){
				LDStream.thin(dataDir+"msprime.geno."+i, dataDir+"msprime.geno.pruned."+i, .2, rgen);
			}
			
	
			
			//compute info
			System.out.println("Computing info");
			for (int i=chrStart; i<chrEnd; i++){
				String inPath = dataDir+"msprime.geno.pruned."+i;
				String outPath = dataDir+"msprime.info.pruned."+i;
				LDStream.writeLdOutfile(inPath, outPath, back);
			}
			
			
			

			
			
			for(int t=0; t<100; t++){
				
				System.out.println(t);
				
				
				// simulation
				System.out.println("Simulating pedigrees");
				for(int i=chrStart; i<chrEnd; i++){
					
					//files
					String unrelated = dataDir+"msprime.geno.pruned."+i;
					
					//simulateCousins(unrelated, simDir+"sim.test.geno."+i, nGen, full, nChildren, nSmallCluster, 8);	
					
					
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
				
				
				
				
				//add error
				System.out.println("Adding error");
				for(int i=chrStart; i<chrEnd; i++){
					LDStream.addError(simDir+"sim.test.geno."+i, simDir+"sim.test.geno.error."+i, seqError, rgen);	
				}
				
				
				

	
				
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginals(dataDir, simDir + "sim.test.geno.error.", dataDir+"msprime.info.pruned.", simDir+"pairwiseLikelihood/test7.marginal."+t, indCols, chrStart, chrEnd);
		
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");
				long startTime = System.nanoTime();		
				computePairwise(dataDir, simDir+"sim.test.geno.error.", dataDir+"msprime.info.pruned.", simDir+"pairwiseLikelihood/test7.pairwise."+t, indCols, relationships, chrStart, chrEnd);	
				System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
				
			
			}

		}
		
		


}
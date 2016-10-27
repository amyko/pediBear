package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStreamPed;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Relationship;
import dataStructures.Path;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class TestLikelihoodLD {
		
		public static double recombRate = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(526564L);


		
		public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, int t, String infoPath) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginal(fileName+"."+t+".tped", infoPath, ids);
			
			
			//open outfile
			PrintWriter writer = DataParser.openWriter(fileName+"."+t+".marginal");				
			
			//write marginals to file
			for (double item : marginals){
				writer.write(String.format("%f\n",item));
			}
			writer.close();
		}
		
		
		public static void computePairwise(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, List<Relationship> relationships, int t, String infoPath) throws IOException{
			
			int numIndiv = ids.length;
			int numRel = relationships.size();
			
			//likelihood
			PrintWriter writer = DataParser.openWriter(fileName+"."+t+".pairwise");
				
			double[][][] lkhd = core.forwardAlgorithm(fileName+"."+t+".tped", infoPath, ids, relationships);

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
			String simDir = dir + "simulations/cousins/";

			
			//individuals
			int numIndiv = 2;
			int[] indCols = new int[numIndiv];
			for(int i=0; i<numIndiv; i++) indCols[i] = i;			
			

			int[] myLengths = new int[]{40,30,20,10};
			String[] r_sqrs = new String[]{"0.100", "0.075", "0.050", "0.025", "0.013"};
			

			
			
			//pairwise core
			PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, recombRate, back, numIndiv);
	
			
			for(int myLength : myLengths){

				for(String r_sqr : r_sqrs){

					//compute info
					//System.out.println("Computing info");
					String infoPath = String.format(dir+"unrelated/msprime.%dmorgan.%s.info", myLength, r_sqr);
					//LDStreamPed.writeLdOutfile(String.format(dir+"unrelated/msprime.%dmorgan.%s", myLength, r_sqr), infoPath, back);	
					
					
					for(int gen=4; gen>=2; gen--){
					
					
						for(int t=0; t<50; t++){
							
							System.out.println(String.format("%d %d %s %d", gen, myLength, r_sqr, t));
							
							String testName = String.format("cousins%d.%dmorgan.%s", gen, myLength, r_sqr);
							
							//compute marginal
							//System.out.println("Computing marginals");
							computeMarginals(core, simDir+testName, indCols, t, infoPath);
							
							
							
							//compute pairwise
							//System.out.println("Computing pairwise likelihoods");
							//long startTime = System.nanoTime();		
							computePairwise(core, simDir+testName, indCols, relationships, t, infoPath);	
							//System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
							
							
							
						}
					}
					
					
				}
				
				
			}
			
			
			
			



			

		}	
		
		


}
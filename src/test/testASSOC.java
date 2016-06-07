package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import simulator.SimulatorStream;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class testASSOC {
		
		public static double r = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(526564L);


	
		
		public static void simulateFullCousins(String inPath, int numGen, int nCluster, int nChildren, int start, int t) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);

			//file names
			String tempOutPath = dir + "assoc/out";
			String genoOutPath = dir + "assoc/geno."+t;
			
			String[] famfiles = new String[nChildren];
			String[] clufiles = new String[nCluster];
			int[][] famcols = new int[nChildren][1];
			int[][] clucols = new int[nCluster][nChildren];
			for(int i=0; i<nChildren; i++){
				famfiles[i] = dir + "assoc/fam."+i;
				famcols[i] = new int[]{0};
			}
			for(int i=0; i<nCluster; i++){
				clufiles[i] = dir + "assoc/clu."+i;
				
				for(int j=0; j<nChildren; j++) 
					clucols[i][j] = j;
				
			}

			
			for(int clu=0; clu<nCluster; clu++){
					
				//first generation
				sim.makeChildren(inPath, inPath, inPath, tempOutPath, start, start+1, rgen, nChildren, "c1\n");
				start+=2;
				
				
				//later generations
				int currGen = numGen - 1;
				while(currGen > 0){
					
		
					for(int i=0; i<nChildren; i++){
						
						sim.makeChildren(inPath, tempOutPath, inPath, famfiles[i], i, start, rgen, 1, "c1\n");
						start++;
						
					}
	
					
					//concatenate into one output file
					if(currGen>1){
						DataParser.concatFiles(famfiles, tempOutPath, famcols);
					}
					else{
						DataParser.concatFiles(famfiles, clufiles[clu], famcols);
					}
					
					//increment
					currGen--;
					
					
				}
				
				
				
			}
			
			
			//concatenate all clusters
			DataParser.concatFiles(clufiles, genoOutPath, clucols);

 
		}
		
		
		public static void makeSibs(String inPath, int nCluster, int nChildren, int start, int t) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);

			//file names
			String genoOutPath = dir + "assoc/geno."+t;
			

			String[] clufiles = new String[nCluster];
			int[][] clucols = new int[nCluster][nChildren];
			for(int i=0; i<nCluster; i++){
				clufiles[i] = dir + "assoc/clu."+i;
				
				for(int j=0; j<nChildren; j++) 
					clucols[i][j] = j;
				
			}

			
			for(int clu=0; clu<nCluster; clu++){
					
				//first generation
				sim.makeChildren(inPath, inPath, inPath, clufiles[clu], start, start+1, rgen, nChildren, "c1\tc2\tc3\tc4\tc5\n");
				start+=2;

				
			}
			
			
			//concatenate all clusters
			DataParser.concatFiles(clufiles, genoOutPath, clucols);

 
		}

		


		public static void main(String[] args) throws IOException{

			//dir
			String dataDir = dir + "assoc/";
			String dataPath = dataDir+"msprime.geno.pruned.1";


			//simulation parameters
			int numGen = 3;
			int nFam = 20;
			int nChildren = 5;
			int start = 3;

	
			//init simulator
			SimulatorStream sim = new SimulatorStream(r);
			
			
			//prune ld
			System.out.println("Pruning LD");
			LDStream.thin(dataDir+"msprime.geno.1", dataPath, .2, rgen);


			//compute info
			System.out.println("Computing info");
			String inPath = dataDir+"msprime.geno.pruned.1";
			String outPath = dataDir+"msprime.info.pruned.1";
			LDStream.writeLdOutfile(inPath, outPath, back);
			
			
			
			
			for(int t=0; t<100; t++){
				System.out.println(t);
				//makeSibs(dataPath, nFam, nChildren, start, t);
				simulateFullCousins(dataPath, numGen, nFam, nChildren, start, t);

			}
			


		}
		
		


}
package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreMicrosatellites;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Relationship;
import dataStructures.Path;
import simulator.SimulatePedigreeGasbarra;
import simulator.SimulatePedigreeUnderPrior;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class TestLikelihood {
				
		//marginal for snp data
		public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, String infoPath, int[] ids, int t) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginalIndep(fileName+".tped", infoPath, ids);
			
			
			//open outfile
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.marginal", fileName, t));	

			//write marginals to file
			for (double item : marginals){
				writer.write(String.format("%f\n",item));
			}
			writer.close();
		}
		
		
		//pairwise for snp data
		public static void computePairwise(PairwiseLikelihoodCoreStreamPed core, String fileName, String infoPath, int[] ids, List<Relationship> relationships, int t) throws IOException{
			
			int numIndiv = ids.length;
			int numRel = relationships.size();
			
			//likelihood
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.pairwise", fileName, t));
	
			double[][][] lkhd = core.forwardAlgorithmIndep(fileName+".tped", infoPath, ids, relationships);

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
		
		
		//margina for microsats
		public static void computeMarginalsMicrosat(PairwiseLikelihoodCoreMicrosatellites core, String fileName, String freqPath, int[] ids, int t) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginalMicrosat(fileName+".tped", freqPath, ids);
			
			
			//open outfile
			//PrintWriter writer = DataParser.openWriter(String.format("%s.%d.marginal", fileName, t));	
			PrintWriter writer = DataParser.openWriter(String.format("%s.marginal", fileName));	//for grant only
			
			
			//write marginals to file
			for (double item : marginals){
				writer.write(String.format("%f\n",item));
			}
			writer.close();
		
		}
		
		
		//pairwise for snp data
		public static void computePairwiseMicrosat(PairwiseLikelihoodCoreMicrosatellites core, String fileName, String freqPath, int[] ids, List<Relationship> relationships, int t) throws IOException{
			
			int numIndiv = ids.length;
			int numRel = relationships.size();
			
			//likelihood
			//PrintWriter writer = DataParser.openWriter(String.format("%s.%d.pairwise", fileName, t));
			PrintWriter writer = DataParser.openWriter(String.format("%s.pairwise", fileName));
			
				
			double[][][] lkhd = core.computePairwiseMicrosat(fileName+".tped", freqPath, ids, relationships);

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
		
		
		
		public static void simulate(int N, int d, double alpha, double beta, Random rgen, String dir, boolean[][] sampled, int sampleSize, int sampleDepth, double recombRate, double epsilon1, double epsilon2, String freqPath, int t) throws IOException {
			
			
			//simulate sample pedigree
			//SimulatePedigreeGasbarra.simulatePedigreeForSamples(N, N/2, N/2, d, alpha, beta, rgen, dir, sampled);
			
			//simulate genotypes according to pedigree given in 200N.sample.fam; outputs ped file containing simulated individuals
			//d = 1 pedigree
			//SimulatePedigreeUnderPrior.simulatePopulationPed(dir+"microsat.20.founders.tped", dir+"20N.sample.fam", dir+"pop.ped", d, recombRate, rgen, N, dir+"temp");
			SimulatePedigreeUnderPrior.simulatePopulationPed(dir+"microsat.50.founders.tped", dir+"gasbarra.3d.sample.fam", dir+"pop.ped", d, recombRate, rgen, N, dir+"temp");
			
			//make tped/tfam for population pedigree; the second argument contains "map" file containing position info
			SimulatePedigreeUnderPrior.ped2tped(dir+"pop", dir+"microsat.50.founders.tped");
					
			//sample individuals. NOTE: the individuals must appear in order. 
			//List<String> samples = SimulatePedigreeUnderPrior.sampleFam(N, sampleSize, d, sampleDepth, rgen, dir);
			List<String> samples = new ArrayList<String>(Arrays.asList("0_0", "0_1", "0_2", "0_3", "0_4", "0_115", "0_116", "0_117", "0_118", "0_119"));

			
			
			//make true file
			//SimulatePedigreeUnderPrior.uninbred(String.format("%s20N.sample.fam", dir), String.format("%ssample.%d.true", dir, t), samples);
			SimulatePedigreeUnderPrior.uninbred(String.format("%sgasbarra.3d.sample.fam", dir), String.format("%ssample.%d.true", dir, t), samples);
			
			
			//make tped/tfam for sampled individuals. Output: temp (tped file), sample.t.tfam
			SimulatePedigreeUnderPrior.makeTpedTfam(N, sampleSize, d, sampleDepth, rgen, dir, t, samples);
		
			
			//add error
			SimulatePedigreeUnderPrior.addError(dir+"temp", freqPath, dir+"sample.tped", epsilon1, epsilon2, rgen);
			
			
			//make colony file
			makeColonyFile(dir, dir+"sample.tped", samples, t);
			
			
			
		}
		
		
		
		public static void makeColonyFile(String dir, String tpedFile, List<String> names, int t) throws IOException {
			
			int nSnp = DataParser.countLines(tpedFile);
			int n = names.size();
			
			String[][] snps = new String[n][2*nSnp];
			
			BufferedReader infile = DataParser.openReader(tpedFile);
			String line;
			
			int j = 0;
			while((line = infile.readLine() )!= null) {
				
				String[] fields = line.split("\\s");
				
				for(int i=0; i< n; i++) {
					
					snps[i][j] = fields[2*i + 4];
					snps[i][j+1] = fields[2*i + 5];
					
				}
				
				j+=2;
				
				
				
			}
			infile.close();
			
			
			//record ped colony file
			PrintWriter writer = DataParser.openWriter(dir + String.format("colony.%d.geno", t));
			
			for(int i=0; i<n; i++) {
				
				String name = names.get(i);
				writer.write(name + " ");
				
				for(j=0; j<snps[0].length; j++) {
					
					writer.write(snps[i][j] + " ");
					
				}
				
				writer.write("\n");
				
			}
			
			writer.close();
					
			
			
		}
		
		
		//for colony file
		public static void addError(String inPath, String outPath, double epsilon1, double epsilon2, Random rGen, int k) throws IOException {
			
			//open files
			BufferedReader reader = DataParser.openReader(inPath);
			PrintWriter writer = DataParser.openWriter(outPath);
			
			//error rates
			double e1 = epsilon1 / (1+epsilon1);
			
			
			String line;

			while((line=reader.readLine()) != null) {
				
				String[] fields = line.split("\\s");
				
				//number of alleles
				String[] alleles = new String[k];
				for(int i=0; i<k; i++) alleles[i] = (i+1) + "";
				
				
				//write header
				writer.write(fields[0] + " ");
				
				
				//for every SNP
				for(int i=0; i<fields.length/2; i++) {
					
					String a1 = fields[2*i + 1];
					String a2 = fields[2*i + 2];
					
					
					//dropout for heterozygote
					if(!a1.equals(a2)) {
						
						if(rGen.nextDouble() < 2*e1) { //drop out
							
							if(rGen.nextDouble() < .5)
								a1 = a2;
							else
								a2 = a1;
							
						}

					}

					
					//sequence error
					if(rGen.nextDouble() < epsilon2) {
						
						int oldIdx = a1.charAt(0) - '0';
						int toAdd = rGen.nextInt(k-1) + 1;
						int newIdx = (oldIdx + toAdd) % k;
						a1 = alleles[newIdx];
						
					}
					
					if(rGen.nextDouble() < epsilon2) {
						
						int oldIdx = a2.charAt(0) - '0';
						int toAdd = rGen.nextInt(k-1) + 1;
						int newIdx = (oldIdx + toAdd) % k;
						a2 = alleles[newIdx];
						
						
					}
					
					writer.write(String.format("%s %s ", a1, a2));	
					
					
				}
				
				writer.write("\n");
				
			}
			

			writer.close();

			
			
		}
		
		
		
		public static void main(String[] args) throws IOException{
			
			////////////////////////////
			
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
			
			/*
			//depth = 3 relationships
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(3,1,2)})); 
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
			*/
			
			//////////////////////////////
	
			
			//parameters
			double recombRate = 1.3e-8; 
			String dir = "/Users/amy/eclipse-workspace/mcmc/simulations/";
			Random rgen = new Random(5265646);
			int N = 200; //population size
			int n = 20; //number of samples in the current generation from which to generate a pedigree
			int d = 2; //maxDepth
			int sampleDepth = 1; //not relevant right now
			double beta = .01; //degree of monogamy
			//double alpha = N * beta / 2; //dominant father
			double alpha = .1;
			int sampleSize = 4;
			double epsilon1 = .001;
			double epsilon2 = .001;

			
			//individuals
			int[] indCols = new int[sampleSize];
			for(int i=0; i<sampleSize; i++) indCols[i] = i;
			
			//samples
			boolean[][] sampled = new boolean[N][N];
			for(int i=0; i<sampleSize; i++) {
				sampled[0][i] = true;
			}
			
			
			
			//SimulatePedigreeGasbarra.simulatePedigreeForSamples(N=20, N/2, N/2, d=2, alpha, beta, rgen, dir, sampled);
			
			
			String testName = "sample";
			
			
			//pairwise core
			//PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, sampleSize);
			PairwiseLikelihoodCoreMicrosatellites core = new PairwiseLikelihoodCoreMicrosatellites(sampleSize, epsilon1, epsilon2);
			
			//compute info
			//System.out.println("Computing info");
			//LDStreamPedMissing.writeLdOutfile(dir+"indepMarkers.founders", dir+"indepMarkers.founders.info", back, true);	
			
			
	
			
			for(int t=0; t<50; t++){
				
				System.out.println(t);
				
				//file names
				String fileName = String.format("%s%s",dir,testName);
				//String infoPath = dir+"indepMarkers.founders.info";
				String freqPath = dir + "microsat.50.founders.freq";
				
				
				System.out.println("Simulating");
				//makeColonyFile(N, d, alpha, beta, rgen, dir, sampled, sampleSize, sampleDepth, recombRate, epsilon1, epsilon2, freqPath, t);
				simulate(N, d, alpha, beta, rgen, dir, sampled, sampleSize, sampleDepth, recombRate, epsilon1, epsilon2, freqPath, t);
				
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginalsMicrosat(core, fileName, freqPath, indCols, t);
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");		
				computePairwiseMicrosat(core, fileName, freqPath, indCols, relationships, t);	

				
				/*
				//// this is for snp data///
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginals(core, fileName, infoPath, indCols, t);
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");
				//long startTime = System.nanoTime();		
				computePairwise(core, fileName, infoPath, indCols, relationships, t);	
				//System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
				*/
				
			
			}
		
			
			
			//true path
			//sim.writeTruePathFile(simDir+"test12.ped", simDir+"10indiv.ped", ids);
			
			
			
			/*
			 //this is for grant
			for(int t=0; t<50; t++){
				
				System.out.println(t);
				
				//file names
				String freqPath = "/Users/amy/eclipse-workspace/mcmc/grant/freq";
				String fileName = "/Users/amy/eclipse-workspace/mcmc/grant/mcmc.";

				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginalsMicrosat(core, fileName + t, freqPath, indCols, t);
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");		
				computePairwiseMicrosat(core, fileName + t, freqPath, indCols, relationships, t);	

			
			} 
			 */
			 
			
			

		}

		
		


}
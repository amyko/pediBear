package test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreMicrosatellites;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Relationship;
import dataStructures.SimplePair;
import dataStructures.Node;
import dataStructures.Path;
import simulator.SimulatePedigreeGasbarra;
import simulator.SimulatePedigreeUnderPrior;
import utility.ArrayUtility;
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
		public static void computeMarginalsMicrosat(PairwiseLikelihoodCoreMicrosatellites core, String fileName, String fPath, int[] ids, boolean withInbreeding, int t) throws IOException{
			
			//compute
			String freqPath = String.format("%s.%d.freq", fileName, t);
			double[] marginals = core.computeMarginalMicrosat(fileName+".tped", freqPath, fPath, ids, withInbreeding);
			
			
			//open outfile
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.marginal", fileName, t));	
			
			
			//write marginals to file
			for (double item : marginals){
				writer.write(String.format("%f\n",item));
			}
			writer.close();
		
		}
		
		
		//pairwise for snp data
		public static void computePairwiseMicrosat(PairwiseLikelihoodCoreMicrosatellites core, String fileName, String fPath, int[] ids, List<Relationship> relationships, boolean withInbreeding, int t) throws IOException{
			
			int numIndiv = ids.length;
			int numRel = relationships.size();
			String freqPath = String.format("%s.%d.freq", fileName, t);
			
			//likelihood
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.pairwise", fileName, t));
			//PrintWriter writer = DataParser.openWriter(String.format("%s.pairwise", fileName));
			
				
			double[][][] lkhd = core.computePairwiseMicrosat(fileName+".tped", freqPath, fPath, ids, relationships, withInbreeding);

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
		
		
		
		public static void simulate(int N, int d, double alpha, double beta, Random rgen, String dir, boolean[][] sampled, int sampleSize, int sampleDepth, double recombRate, double epsilon1, double epsilon2, String freqPath, boolean allowCycles, boolean estimateAF, int nSnp, int t) throws IOException {
			
			
			//simulate sample pedigree
			SimulatePedigreeGasbarra.simulatePedigreeForSamples(N, N/2, N/2, d, alpha, beta, rgen, dir, sampled, allowCycles);
			
			//simulate genotypes according to pedigree given in gasbarra.fam; outputs ped file containing simulated individuals
			SimulatePedigreeUnderPrior.simulatePopulationPed(String.format("%smicrosat.%dloci.%dN.founders.tped", dir, nSnp, N), dir+"gasbarra.fam", dir+"pop.ped", d, recombRate, rgen, N, dir+"temp");
			
			//make tped/tfam for population pedigree; the second argument contains "map" file containing position info
			SimulatePedigreeUnderPrior.ped2tped(dir+"pop", String.format("%smicrosat.%dloci.%dN.founders.tped", dir, nSnp, N));
					
			//sample individuals. NOTE: the individuals must appear in order.
			List<String> samples = new ArrayList<String>();
			for(int i=0; i<sampleSize; i++) {
				samples.add(String.format("0_%d", i));
			}

			//make true file
			//SimulatePedigreeUnderPrior.uninbred(String.format("%sgasbarra.fam", dir), String.format("%ssample.%d.true", dir, t), samples);
			makeTrueInbred(dir+"gasbarra.fam", samples, String.format("%ssample.%d.true", dir, t));
			
			
			//make tped/tfam for sampled individuals. Output: temp (tped file), sample.t.tfam
			SimulatePedigreeUnderPrior.makeTpedTfam(N, sampleSize, d, sampleDepth, rgen, dir, t, samples);
		
			
			//add error
			SimulatePedigreeUnderPrior.addError(dir+"temp", freqPath, dir+"sample.tped", epsilon1, epsilon2, rgen);
			
			
			//estimate allele frequency
			estimateAlleleFreq(dir+"sample.tped", freqPath, String.format("%ssample.%d.freq", dir, t), estimateAF);
			
			
			//make f file
			makeFfile(dir + "sample.tped", String.format("%ssample.%d.freq", dir, t), String.format("%ssample.%d.f", dir, t), sampleSize);
			
			
			//make colony file
			makeColonyFile(dir, dir+"sample.tped", samples, t);
			
			
			
		}
		

		public static void estimateAlleleFreq(String tpedPath, String popFreqPath, String outPath, boolean estimate) throws IOException {
			
			//open writer
			PrintWriter writer = DataParser.openWriter(outPath);
			BufferedReader popFreq = DataParser.openReader(popFreqPath);

			
			// use population frequencies
			if(!estimate) {
				String line;
				while((line = popFreq.readLine()) != null) {
					
					int k = line.split("\\s").length;
					String alleles = "";
					for(int i = 1; i<=k; i++) alleles += i+" ";
					
					writer.write(alleles+"\n");
					writer.write(line+"\n");
				}
			}
			
			// estimate from sample
			else {
				
				BufferedReader tped = DataParser.openReader(tpedPath);
				
				String line;
				
				while((line = tped.readLine()) != null) {
					
					Map<String, Integer> count = new HashMap<String, Integer>();
					

					double totalCount = 0;	
					
					
					String[] fields = line.split("\\s");
					for(int i=4; i<fields.length; i++) {
						
						String allele = fields[i];
						if(!count.containsKey(allele)) count.put(allele, 0);
						count.put(allele, count.get(allele)+1);
						totalCount++;
						
					}
					
					//write this locus
					String alleles = "";
					String freqs = "";
					for(String a : count.keySet()) {
						alleles += a + " ";
						freqs += String.format("%f ", count.get(a) / totalCount);
					}
					writer.write(alleles+"\n");
					writer.write(freqs+"\n");
					
					
				}
				
			}
			
			writer.close();
			
			
			
			
		}
		
		//get IBS
		

		
		//estimate inbreeding coefficient: f = (O(hom) - E(hom)) / (m - E(hom))
		public static void makeFfile(String tpedPath, String freqPath, String outPath, int n) throws IOException {
			
	
			// open files
			BufferedReader genoFile = DataParser.openReader(tpedPath);
			BufferedReader freqFile = DataParser.openReader(freqPath);
			
			int m = 0; // number of loci
			int[] obs = new int[n]; // observed number of homozygotes for each individual
			double[] exp = new double[n]; // expected number of homozygotes
			
			String genoLine;

			
			while((genoLine = genoFile.readLine()) != null) {
				
				freqFile.readLine();
				String[] freq = freqFile.readLine().split("\\s");
				String[] geno = genoLine.split("\\s");
				m++;
				
				for(int i=0; i<n; i++) {
					
					String g1 = geno[2*i + 4];
					String g2 = geno[2*i + 5];
					
					// count observed
					if(g1.equals(g2)) obs[i]++;
					
					// what's expected
					for(int j = 0; j<freq.length; j++)
						exp[i] += Math.pow(Double.parseDouble(freq[j]), 2);
					
					
				}	
				
			}
			
			//close files
			genoFile.close();
			freqFile.close();
			
			
			//write f
			PrintWriter writer = DataParser.openWriter(outPath);
			
			for(int i=0; i<n; i++) {
				
				double f = (obs[i] - exp[i]) / (m - exp[i]);
				f = Math.max(f, 0);
				//f = 0;
				writer.write(f+"\n");
				
			}
			
			writer.close();
			
			
			
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
		
		
		
		//true file for inbred pedigree; make the most recent relationship as true
		public static void makeTrueInbred(String famFile, List<String> samples, String outPath) throws IOException {
			
			int n = samples.size();
			Path[][] rel = new Path[n][n];
			Map<String, SimplePair<String, String>> id2parents = new HashMap<String, SimplePair<String, String>>(); // id -> (fid, mid)
			Map<String, ArrayList<String>> id2children = new HashMap<String, ArrayList<String>>();
			Map<String, Integer> id2idx = new HashMap<String, Integer>();
			List<String> parents = new ArrayList<String>(); //generation 1 individuals
			
			//read relationships; make nodes with children
			BufferedReader reader = DataParser.openReader(famFile);
			reader.readLine(); //skip header
			String line;
			int idx = 0;
			while((line = reader.readLine()) != null) {
				String[] fields = line.split("\\s");
				// record parents
				id2parents.put(fields[0], new SimplePair<String, String>(fields[1], fields[2]));
				if(fields[0].charAt(0)=='1') parents.add(fields[0]);
				// record children
				if(id2children.get(fields[1]) == null) id2children.put(fields[1], new ArrayList<String>());
				if(id2children.get(fields[2]) == null) id2children.put(fields[2], new ArrayList<String>());
				id2children.get(fields[1]).add(fields[0]);
				id2children.get(fields[2]).add(fields[0]);
				//record index
				id2idx.put(fields[0], idx);
				idx++;
				
			}
			reader.close();
			
			
			//record FS or HS in the current generation
			for(int i=0; i<n; i++) {
				
				SimplePair<String, String> p1 = id2parents.get(samples.get(i));
				
				for(int j=i+1; j<n; j++) {
					
					SimplePair<String, String> p2 = id2parents.get(samples.get(j));
					
					//FS
					if(p1.getFirst().equals(p2.getFirst()) && p1.getSecond().equals(p2.getSecond())) {
						rel[i][j] = new Path(1,1,2);		
						rel[j][i] = new Path(1,1,2);	
					}
					//HS
					else if(p1.getFirst().equals(p2.getFirst()) || p1.getSecond().equals(p2.getSecond())) {
						rel[i][j] = new Path(1,1,1);
						rel[j][i] = new Path(1,1,1);	
					}
					else {
						rel[i][j] = new Path(0,0,0);
						rel[j][i] = new Path(0,0,0);	
					}
					
				}
				
			}
			
			int numRel = 0;
			int numInbred = 0;
			
			
			//record FC, HC
			for(String id1 : parents) {
				
				SimplePair<String, String> p1 = id2parents.get(id1);
				
				for(String id2 : parents) {
					
					if(id1.equals(id2)) continue; 
					
					SimplePair<String, String> p2 = id2parents.get(id2);
					
					
					//FC
					if(p1.getFirst().equals(p2.getFirst()) && p1.getSecond().equals(p2.getSecond())) {
						
						boolean inbred = false;
						numRel++;

						//children of id1 and id2 are full cousins
						for(String c1 : id2children.get(id1)) {
							
							int i = id2idx.get(c1);
							
							for(String c2 : id2children.get(id2)) {
								
								int j = id2idx.get(c2);
								
								if (i>= j) continue;
								
								// unrelated or HC/FC
								if(rel[i][j].getNumVisit() != 0)
									inbred = true;
								if(rel[i][j].getNumVisit() == 0 || rel[i][j].getUp() == 2) {
									rel[i][j] = new Path(2,2,2);
									rel[j][i] = new Path(2,2,2);
								}
								
								
							}
						}
						
						if(inbred) numInbred++;

		
					}
					//HC
					else if(p1.getFirst().equals(p2.getFirst()) || p1.getSecond().equals(p2.getSecond())) {

						numRel++;
						boolean inbred = false;
						
						//children of id1 and id2 are full cousins
						for(String c1 : id2children.get(id1)) {
							
							int i = id2idx.get(c1);
							
							for(String c2 : id2children.get(id2)) {
								
								int j = id2idx.get(c2);
								
								if(i >= j) continue;
								
								// unrelated so far
								if(rel[i][j].getNumVisit() != 0)
									inbred = true;
								if(rel[i][j].getNumVisit() == 0) {
									rel[i][j] = new Path(2,2,1);
									rel[j][i] = new Path(2,2,1);
								}
	

								
							}
							
						}
						
						if(inbred) numInbred++;

					}
			
					
					
					
				}
				
			}
			
			
			//write true
			PrintWriter writer = DataParser.openWriter(outPath);

			
			//write true
			String pairwise = "";
			for(int i=0; i<rel.length; i++){

				for(int j=i+1; j<rel[0].length; j++){
					
					Path r = rel[i][j];
					
					pairwise += String.format("%d%d%d", r.getUp(),r.getDown(),r.getNumVisit());
					
				}
				
			}
			
			writer.write(pairwise);
			writer.close();
			
			System.out.println(String.format("%d %d %d %f", n + parents.size(), numRel, numInbred, numInbred / (double) numRel));
			

		}
		
		
		// testing prior distribution
		public static void test(int N, double alpha, double beta, Random rgen, String dir, boolean[][] sampled, int nIter) throws IOException {
			
			// 3 samples
			List<String> samples = new ArrayList<String>(Arrays.asList("0_0", "0_1", "0_2"));
			
			//distribution of pedigrees
			Map<String, Integer> ped2count = new HashMap<String, Integer>();
			
			for(int i=0; i<nIter; i++) {
			
				if(i % 1000 == 0)
					System.out.println(i);
				
				//simulate sample pedigree
				SimulatePedigreeGasbarra.simulatePedigreeForSamples(N, N/2, N/2, 1, alpha, beta, rgen, dir, sampled, false);
			
				//make true 
				makeTrueInbred(dir+"gasbarra.fam", samples, String.format("%ssample.true", dir));
				
				//read true
				String ped = DataParser.openReader(String.format("%ssample.true", dir)).readLine();
				
				//increment
				int count = ped2count.get(ped) == null? 0 : ped2count.get(ped);
				ped2count.put(ped, count+1);
				
			}
			
			//write results
			PrintWriter writer = DataParser.openWriter(dir+"ped.dist.simulated");
			for(String key : ped2count.keySet())
				writer.write(String.format("%s\t%d\n", key, ped2count.get(key)));
			
			writer.close();
			
			
		}
		
		public static void main(String[] args) throws IOException{
			
			////////////////////////////
			
			//relationships
			List<Relationship> relationships = new ArrayList<Relationship>();

			relationships.add(new Relationship(13d, new double[] {1d, 0d, 0d}, new Path[]{new Path(0,0,0)})); //unrelated
			
			// depth = 1 relationship
			relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
			relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
			relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1)})); //half-sib
			
			//depth = 2 relationships
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0}, new Path[]{new Path(2,1,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2)})); //first cousins
			
			/*
			relationships.add(new Relationship(13d, new double[] {1d, 0d, 0d}, new Path[]{new Path(0,0,0)})); //unrelated
	
			// depth = 1 relationship
			relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
			relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
			relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
			
			//depth = 2 relationships
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2)})); //first cousins
			
			
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
			String dir = "/Users/amy/eclipse-workspace/mcmc/simulations2/";
			Random rgen = new Random(5265646); 
			int N = 1000; //population size
			int nSnp = 20;
			int d = 2; //maxDepth
			int sampleDepth = 1; //not relevant right now
			double alpha = .2;
			double beta = .01; //degree of monogamy; but can't control polygamy of mothers 
			int sampleSize = 20;
			double epsilon1 = .05; //allele drop out rate (.05)
			double epsilon2 = .02; //other errors (.02)
			boolean allowCycles = true; // allow cycles in simulation
			boolean withInbreeding = true; // compute pairwise likelihood with inbreeding
			boolean estimateAF = false;

			
			//individuals
			int[] indCols = new int[sampleSize];
			for(int i=0; i<sampleSize; i++) indCols[i] = i;


			
			//SimulatePedigreeGasbarra.simulatePedigreeForSamples(N=200, N/2, N/2, d=2, alpha, beta, rgen, dir, sampled);
			
			
			String testName = "sample";
			
			
			//pairwise core
			//PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, back, sampleSize);
			PairwiseLikelihoodCoreMicrosatellites core = new PairwiseLikelihoodCoreMicrosatellites(sampleSize, epsilon1, epsilon2);
			
			//compute info
			//System.out.println("Computing info");
			//LDStreamPedMissing.writeLdOutfile(dir+"indepMarkers.founders", dir+"indepMarkers.founders.info", back, true);	

			
			for(int t=0; t<10; t++){
				
				
				System.out.println(t);
				
				//samples
				//TODO undo the comment below
				boolean[][] sampled = new boolean[N][N];
				for(int i=0; i<sampleSize; i++) {
					sampled[0][i] = true;
				}
				
				//file names
				String fileName = String.format("%s%s",dir,testName);
				//String infoPath = dir+"indepMarkers.founders.info";
				String popFreqPath = String.format("%smicrosat.%dloci.%dN.founders.freq", dir, nSnp, N);
				String fPath = String.format("%ssample.%d.f", dir, t);
				
				
				System.out.println("Simulating");
				//makeColonyFile(N, d, alpha, beta, rgen, dir, sampled, sampleSize, sampleDepth, recombRate, epsilon1, epsilon2, freqPath, t);
				simulate(N, d, alpha, beta, rgen, dir, sampled, sampleSize, sampleDepth, recombRate, epsilon1, epsilon2, popFreqPath, allowCycles, estimateAF, nSnp, t);
				
				
				//compute marginal probs
				System.out.println("Computing marginals");
				computeMarginalsMicrosat(core, fileName, fPath, indCols, withInbreeding, t);
				
				
				//compute pairwise
				System.out.println("Computing pairwise likelihoods");		
				computePairwiseMicrosat(core, fileName, fPath, indCols, relationships, withInbreeding, t);	

				
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


		}

		
		


}
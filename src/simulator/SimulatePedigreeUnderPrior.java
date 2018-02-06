package simulator;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.SimplePair;
import utility.DataParser;

public class SimulatePedigreeUnderPrior {
	
	
	//sample individuals from first d generations
	public static Set<String> sample(int N, int n, int d, Random rGen){
		
		Set<String> toReturn = new HashSet<String>();
		int numSampled = 0;
		double totalUnits = 0.0;
		for(int i=1; i<d+2; i++) totalUnits += i;
		
		while(numSampled < n){
			
			//choose depth
			double p = rGen.nextDouble();
			double cdf = 0;
			int depth = -1;
						
			for(int i=1; i<d+2; i++){
				
				cdf += i/totalUnits;
				
				if(p < cdf){
					depth = d - i + 1;
					break;
				}
				
			}
			
			
			//choose individual
			int indiv = rGen.nextInt(N);
			String id = String.format("%d_%d", depth, indiv);
			
			if(toReturn.contains(id)) continue;
			
			numSampled++;
			
			toReturn.add(id);
			
			//System.out.println(id);
			
		}
		

		return toReturn;
		
		
		
		
	}
	
	
	
	
	//sample individuals in the fam file uniformly at random
	public static List<String> sampleUniform(String famPath, Random rGen, int n) throws IOException{
		
		List<String> toReturn = new ArrayList<String>();
		List<String> candidates = new ArrayList<String>();
		
		BufferedReader infile = DataParser.openReader(famPath);
		infile.readLine(); //skip header
		
		
		//get candidate individuals
		String line;
		
		while((line=infile.readLine())!=null) {
			
			candidates.add(line.split("\\s")[0]);
			
		}
		infile.close();
		
		
		//sample n 
		int nSampled = 0;
		int nLeft = candidates.size();
		int nSampleNeeded = n;
		int i=0;
		while(nSampled < n) {
			
			int u = rGen.nextInt(nLeft);
			
			if(u < nSampleNeeded) {
				toReturn.add(candidates.get(i));
				nSampleNeeded--;
				nSampled++;
			}
				
			i++;
			nLeft--;
			
		}

		//TODO testing
		//List<String> testReturn = new ArrayList<String>(Arrays.asList("0_3","0_13","1_39","1_86","1_146"));
		//return testReturn;
		
		return toReturn;
		
		
		
		
		
	}
	
	//sample connected individuals from first d generations
	public static Set<String> sampleConnected(int N, int n, int d, Random rGen){
		
		Set<String> toReturn = new HashSet<String>();
		int numSampled = 0;
		double totalUnits = 0.0;
		for(int i=1; i<d+2; i++) totalUnits += i;
		
		while(numSampled < n){
			
			//choose depth
			double p = rGen.nextDouble();
			double cdf = 0;
			int depth = -1;
						
			for(int i=1; i<d+2; i++){
				
				cdf += i/totalUnits;
				
				if(p < cdf){
					depth = d - i + 1;
					break;
				}
				
			}
			
			
			//choose individual
			int indiv = rGen.nextInt(N);
			String id = String.format("%d_%d", depth, indiv);
			
			if(toReturn.contains(id)) continue;
			
			numSampled++;
			
			toReturn.add(id);
			
			//System.out.println(id);
			
		}
		
		
		return toReturn;
		
		
		
		
	}
	
	
	//construct pedigree under random mating
	public static void polygamous(int N, int d, Random rGen, String outPath) throws IOException{
		
		PrintWriter writer = DataParser.openWriter(String.format("%s%dN.%dd.pop.fam", outPath, N,d));
		writer.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		
		//for each depth
		for(int i=0; i<=d; i++){
			
			//for each individual
			for(int j=0; j<N; j++){
				
				String id = String.format("%d_%d", i,j);
				int sex = j%2==0? 1 : 7;
				
				String fid = "0";
				String mid = "0";
				
				if(i!=d){
					fid = String.format("%d_%d", i+1, rGen.nextInt(N/2));
					mid = String.format("%d_%d", i+1, rGen.nextInt(N/2) + N/2);
				}
				
				String sampled = "999999";
				
				writer.write(String.format("%s\t%s\t%s\t%d\t%s\n", id, fid, mid, sex, sampled));
				
				
				
			}
			
		}
		
		writer.close();
		
		
		
	}
	
	
	
	public static List<String> sampleFam(int N, int n, int d, int sampleDepth, Random rGen, String filePath) throws IOException{
		
		//sample from population
		List	<String> samples = sampleUniform(String.format("%s%dN.%dd.pop.fam", filePath, N, d), rGen, n);
		Set<String> ghosts = new HashSet<String>();
		
		//open pop fam
		BufferedReader reader = DataParser.openReader(String.format("%s%dN.%dd.pop.fam", filePath, N, d));
		PrintWriter writer = DataParser.openWriter(String.format("%s%dN.sample.fam", filePath, N));
		
		//header
		reader.readLine();
		writer.write(String.format("NAME\tFATHER\tMOTHER\tSEX\tSAMPLED\n"));
		
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//if connected to sample node
			if(samples.contains(fields[0])){
			
				fields[4] = "000000";
				
				for(String x : fields)
					writer.write(String.format("%s\t", x));
				
				writer.write("\n");

				
				ghosts.add(fields[1]);
				ghosts.add(fields[2]);
				
				
			}
			
			else if(ghosts.contains(fields[0])){
				
				
				for(String x : fields)
					writer.write(String.format("%s\t", x));
				
				writer.write("\n");

				
				ghosts.add(fields[1]);
				ghosts.add(fields[2]);
				
				
			}
			
			
		}
		
		reader.close();
		writer.close();
	
		return samples;
		
	}
	
	
	
	//make tped, and tfam files
	public static void makeTpedTfam(int N, int n, int d, int sampleDepth, Random rGen, String filePath, int t, List<String> samples) throws IOException{
				
		//make tfam files: individuals must appear in order
		BufferedReader reader = DataParser.openReader(String.format("%spop.ped", filePath));
		PrintWriter writer = DataParser.openWriter(String.format("%ssample.%d.tfam", filePath, t));
		String line;
		int lineNum = 0;
		
		//map name --> column number
		Map<String, SimplePair<String, Integer>> name2info = new HashMap<String, SimplePair<String, Integer>>();
		
		//iterate thorugh ped file. its header is same as the tfam file
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			name2info.put(fields[1], new SimplePair<String, Integer>(line, lineNum));
			
			lineNum++;
				
		}
		
		reader.close();
		
		//get relevant column numbers, write tfam in order
		int[] ids = new int[n]; //0 through ..
		int idx = 0;
		//write tfam and get relevant column numbers
		for(String id : samples) {
			
			//write tfam
			String[] fields = name2info.get(id).getFirst().split("\\s");
			for(int i=0; i<6; i++) 
				writer.write(fields[i] + " ");
			writer.write("\n");
			
			//column number
			ids[idx] = name2info.get(id).getSecond();
			idx++;
			
		}
		writer.close();
		
		
		
		reader = DataParser.openReader(String.format("%spop.tped", filePath));
		writer = DataParser.openWriter(String.format("%stemp", filePath, t));
		
		
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			//write header
			writer.write(String.format("%s %s %s %s ", fields[0],fields[1],fields[2],fields[3]));
			
			
			for(int i=0; i<ids.length; i++){
				
				String a1 = fields[2*ids[i]+4];
				String a2 = fields[2*ids[i]+5];
				
				writer.write(String.format("%s %s ", a1, a2));
				
			}
			
			writer.write("\n");
			
		}
		
		writer.close();
		reader.close();
		
		
		

		
	}
	
	
	
	//add error. NOTE: ONLY WORKS FOR NUMERIC ENCODING OF ALLELES
	public static void addError(String tpedPath, String freqPath, String outPath, double epsilon1, double epsilon2, Random rGen) throws IOException {
		
		//open files
		BufferedReader reader = DataParser.openReader(tpedPath);
		BufferedReader reader2 = DataParser.openReader(freqPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//error rates
		double e1 = epsilon1 / (1+epsilon1);
		
		
		String line;
		String line2;
		
		while((line=reader.readLine()) != null && (line2=reader2.readLine()) != null) {
			
			String[] fields = line.split("\\s");
			
			//number of alleles
			int k = line2.split("\\s").length;
			String[] alleles = new String[k];
			for(int i=0; i<k; i++) alleles[i] = (i+1) + "";
			
			
			//write header
			writer.write(String.format("%s %s %s %s ", fields[0], fields[1], fields[2], fields[3]));
			
			
			//for every ind in at ths locus
			for(int i=0; i<fields.length/2-2; i++) {
				
				String a1 = fields[2*i + 4];
				String a2 = fields[2*i + 5];
				
				
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
	
	
	//disconnect inbred connections
	public static Pedigree uninbred(String filePath, String truePath, List<String> samples) throws IOException{
		
		//make pedigree graph
		Pedigree ped = new Pedigree(filePath, truePath, samples);
		return ped;
		
		
	}
	
	
	
	public static void simulatePopulationPed(String founderGenotypePath, String pedPath, String outPath, int depth, double recomb, Random rGen, int N, String tempPath) throws IOException{
		
		//init simulator
		SimulatorStreamPed sim = new SimulatorStreamPed(recomb);
		
		//maps name to column
		Map<String, Integer> name2id = new HashMap<String, Integer>();
		
		//map for founder generation
		for(int i=0; i<N; i++) {
			String key = String.format("%d_%d", depth, i);
			name2id.put(key, i);
		}
		
		//init writer
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//first generation
		simulateOneGeneration(sim, writer, pedPath, founderGenotypePath, depth-1, rGen, N, tempPath, name2id);
		
		
		//other generations
		for(int i=depth-2; i>=0; i--){
			
			simulateOneGeneration(sim, writer, pedPath, tempPath, i, rGen, N, tempPath, name2id);
			
			
		}
		
		writer.close();
		
		
	}
	
	
	public static void simulateOneGeneration(SimulatorStreamPed sim, PrintWriter writer, String pedPath, String parentPath, int childDepth, Random rGen, int N, String tempPath, Map<String, Integer> name2id ) throws IOException{
		
		
		//open reader
		BufferedReader reader = DataParser.openReader(pedPath);
		reader.readLine();
	
			
		int nSnp = DataParser.countLines(parentPath);
		String[][] tped = new String[2*nSnp][N];
		String[] header = new String[nSnp];
		
		String line;
		int id=0;
		while((line=reader.readLine())!=null){
		
			String[] fields = line.split("\t");
			
			//check for current child generation
			int d = Integer.parseInt(fields[0].split("_")[0]);
			if(d!=childDepth) continue;
			
			//TODO fix this
			//get mom and dad IDs
			int dadID = name2id.get(fields[1]);
			int momID = name2id.get(fields[2]);
			
			//int dadID = Integer.parseInt(fields[1].split("_")[1]);
			//int momID = Integer.parseInt(fields[2].split("_")[1]);
			
			//simulate
			String pedLine = sim.makeChildrenReturnPedline(parentPath, parentPath, momID, dadID, rGen, 1);
			int sex = fields[3].equals("1") ? 1:2;
			
			
			//store tped; for every snp
			String[] geno = pedLine.split("\\s");
			for(int i=0; i<geno.length/2; i++){
				tped[2*i][id] = geno[2*i];
				tped[2*i+1][id] = geno[2*i+1];
			}
			
			
			//write
			writer.write(String.format("1 %s %s %s %d -9 %s\n", fields[0], fields[1], fields[2], sex, pedLine));
			writer.flush();
			
			//update
			name2id.put(fields[0], id); //child
			id++;
			//System.out.println(fields[0]);


			
		}
		
		reader.close();
		
			
		//get header file
		reader = DataParser.openReader(parentPath);
		for(int i=0; i<header.length; i++){
			String[] fields = reader.readLine().split("\\s");
			header[i] = String.format("%s %s %s %s ", fields[0], fields[1], fields[2], fields[3]);
		}
		reader.close();
		
	
		//write temp file
		PrintWriter tempFile = DataParser.openWriter(tempPath);
	
		
		//for every snp
		for(int i=0; i<tped.length/2; i++){
				
			//write header
			tempFile.write(header[i]);
			
			for(int j=0; j<tped[0].length; j++){
				
				tempFile.write(String.format("%s %s ",tped[2*i][j], tped[2*i+1][j]));
				
			}
			tempFile.write("\n");
			
		}
		
		reader.close();
		tempFile.close();
			
	}		

		

	public static void ped2tped(String fileName, String mapPath) throws IOException {
		
		BufferedReader infile = DataParser.openReader(fileName+".ped");
		int nSnp = DataParser.countLines(mapPath); 
		int nIndiv = DataParser.countLines(fileName+".ped");

		
		//store snps
		String[][] ped = new String[nIndiv][2*nSnp];
		int i = 0;
		
		String line;
		while((line=infile.readLine())!=null) {
			
			String[] fields = line.split("\\s");
			
			for(int j=0; j<2*nSnp; j++) {
				
				ped[i][j] = fields[j+6];
				
			}
			
			i++;
			
		}
		infile.close();
		
		
		
		//write tped
		PrintWriter outfile = DataParser.openWriter(fileName+".tped");
		infile = DataParser.openReader(mapPath);
		int j = 0;
		
		while((line=infile.readLine())!=null) {
			
			String[] fields = line.split("\\s");
			
			//writer header
			for(int k=0; k<4; k++) outfile.write(fields[k]+" ");
			
			
			for(i=0; i<ped.length; i++) {
				
				outfile.write(ped[i][2*j]+" "+ped[i][2*j+1]+" ");
				
			}
			outfile.write("\n");
			j++;
			
		}
		infile.close();
		outfile.close();
	}
	
	
	
	public static void main(String[] args) throws IOException{
		
		int N = 200;
		int n = 20;
		int d = 3;
		int sampleDepth = 2;
		double recomb = 1.3e-8;
		Random rGen = new Random(1239);
		double seqError = .0;
		String resultDir = "/Users/kokocakes/Google Drive/Research/pediBear/data/mcmc/";
		String dataDir = "/Users/kokocakes/Google Drive/Research/pediBear/data/mcmc/";
		
		//make population pedigree
		//polygamous(N,d,rGen, resultDir);
		
		/*
		String founderPath = String.format("%s%dfounders.%dN.10k.tped", dataDir, N, N);
		String famPath = String.format("%s%dN.%dd.pop.fam", resultDir, N, d);
		String outPath = String.format("%s%dN.pop.10k.ped", resultDir, N);
		String tempPath = String.format("%stemp.tped", resultDir);
		simulatePopulationPed(founderPath, famPath, outPath, d, recomb, rGen, N, tempPath);
		*/
		

		//sampleForReal(N, n, d, sampleDepth, rGen, resultDir, seqError, 0);
		
		
	
		
	
		
	}
	

}

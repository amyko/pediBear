package simulator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import dataStructures.Path;
import dataStructures.Pedigree;
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
			
			System.out.println(id);
			
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
	
	
	
	public static Set<String> sampleFam(int N, int n, int d, int sampleDepth, Random rGen, String filePath) throws IOException{
		
		//sample from population
		Set<String> samples = sample(N,n,sampleDepth,rGen);
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
	
	
	
	//sample, make true, tped, and tfam files
	public static void sampleForReal(int N, int n, int d, int sampleDepth, Random rGen, String filePath, double seqError) throws IOException{
		
		boolean badSample = true;
		Set<String> samples = null;
		Pedigree ped = null;
		
		while(badSample){
		
			//make fam file of the samples
			samples = sampleFam(N,n,d,sampleDepth,rGen,filePath);
		
			//un-inbred
			ped = uninbred(String.format("%s%dN.sample.fam", filePath, N), String.format("%ssample.true", filePath));
			
			badSample = ped.looped;
			
		}
		
		//relevant columns
		int[] ids = new int[n]; //0 through ...
		
		//make tfam files
		PrintWriter writer = DataParser.openWriter(String.format("%ssample.tfam", filePath, N));
		BufferedReader reader = DataParser.openReader(String.format("%s%dN.pop.tfam", filePath, N));
		String line;
		int lineNum = 0;
		int idx = 0;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			if(samples.contains(fields[1])){
				writer.write(line+"\n");
				ids[idx] = lineNum; //column number for tped
				idx++;
			}
			
			lineNum++;
				
		}
		
		reader.close();
		writer.close();
		
		for(int x : ids) System.out.println(x);
		
		
		reader = DataParser.openReader(String.format("%s%dN.pop.tped", filePath, N));
		writer = DataParser.openWriter(String.format("%ssample.tped", filePath));
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			//write header
			writer.write(String.format("%s %s %s %s ", fields[0],fields[1],fields[2],fields[3]));
			
			for(int i=0; i<ids.length; i++){
				
				String a1 = fields[2*i+4];
				String a2 = fields[2*i+5];
				
				//add error
				if(rGen.nextDouble() < seqError){
					a1 = a1.equals("A") ? "T" : "A";
				}
				
				if(rGen.nextDouble() < seqError){
					a2 = a2.equals("A") ? "T" : "A";
				}
				
				
				writer.write(String.format("%s %s ", a1, a2));
				
			}
			
			writer.write("\n");
			
		}
		
		writer.close();
		reader.close();
		
		
		

		
	}
	
	
	//disconnect inbred connections
	public static Pedigree uninbred(String filePath, String truePath) throws IOException{
		
		//make pedigree graph
		Pedigree ped = new Pedigree(filePath, truePath);
		return ped;
		
		
	}
	
	
	
	public static void simulatePopulationPed(String founderGenotypePath, String pedPath, String outPath, int depth, double recomb, Random rGen, int N, String tempPath) throws IOException{
		
		//init simulator
		SimulatorStreamPed sim = new SimulatorStreamPed(recomb);
		
		//init writer
		PrintWriter writer = DataParser.openWriter(outPath);
		
		//first generation
		simulateOneGeneration(sim, writer, pedPath, founderGenotypePath, depth-1, rGen, N, tempPath);
		
		
		//other generations
		for(int i=depth-2; i>=0; i--){
			
			simulateOneGeneration(sim, writer, pedPath, tempPath, i, rGen, N, tempPath);
			
			
		}
		
		writer.close();
		
		
	}
	
	
	public static void simulateOneGeneration(SimulatorStreamPed sim, PrintWriter writer, String pedPath, String parentPath, int childDepth, Random rGen, int N, String tempPath) throws IOException{
		
		
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
			
			//get mom and dad IDs
			int dadID = Integer.parseInt(fields[1].split("_")[1]);
			int momID = Integer.parseInt(fields[2].split("_")[1]);
			
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
			
			id++;
			System.out.println(id);


			
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

		
		
	
	
	
	
	public static void main(String[] args) throws IOException{
		
		int N = 200;
		int n = 20;
		int d = 3;
		int sampleDepth = 2;
		double recomb = 1.3e-8;
		Random rGen = new Random(1239);
		double seqError = .01;
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
		

		sampleForReal(N, n, d, sampleDepth, rGen, resultDir, seqError);
		
		
	
		
	
		
	}
	

}

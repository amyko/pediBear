package test;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import simulator.SimulatorStreamPed;
import utility.DataParser;

public class TestSim2 {
	
	static double recomb = 1.3e-8;
	static String dir = "/Users/kokocakes/Google Drive/Research/pediBear/data/unrelated/";
	static Random rGen = new Random(10479565L);
	
	
	public static int fam1(SimulatorStreamPed sim, String tpedPath, int startID) throws IOException{
		
		//first generation: make 2 children
		sim.makeChildren(tpedPath, tpedPath, dir+"gen1.out", startID, startID+1, rGen, 2);
		startID+=2;
		
		//for each branching family
		int id = 0;
		for(int i=0; i<2; i++){
			
			//make 2 children
			sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", i, startID, rGen, 2);
			startID++;
			
			//for each child
			for(int j=0; j<2; j++){
				String tempPath = String.format(dir+"c%d.out", id);
				id++;
				sim.makeChildren(dir+"gen2.out", tpedPath, tempPath, j, startID, rGen, 1);
				startID++;
			}			
		}
		
		
		//concatenate final children
		DataParser.concatFilesSpace(new String[]{dir+"c0.out", dir+"c1.out", dir+"c2.out",dir+"c3.out"}, dir+"fam1.out", new int[][]{{}, {4,5},{4,5},{4,5}});
		
		return startID;
		
	}
	
	
	
	public static int fam2(SimulatorStreamPed sim, String tpedPath, int startID) throws IOException{
		
		//first generation: make 2 children
		sim.makeChildren(tpedPath, tpedPath, dir+"gen1.out", startID, startID+1, rGen, 2);
		startID+=2;
		
		//first branching family
		//make 2 children
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", 0, startID, rGen, 2);
		startID++;
		
		//for each child
		for(int j=0; j<2; j++){
			sim.makeChildren(dir+"gen2.out", tpedPath, String.format(dir+"c%d.out", j), j, startID, rGen, 1);
			startID++;
		}
		
		//second branching family
		//make 1 child and save it
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"c2.out", 1, startID, rGen, 1);
		startID++;
		
		//make another child
		sim.makeChildren(dir+"c2.out", tpedPath, dir+"c3.out", 0, startID, rGen, 1);
		startID++;
		
		//concat files
		DataParser.concatFilesSpace(new String[]{dir+"c0.out", dir+"c1.out", dir+"c2.out",dir+"c3.out"}, dir+"fam2.out", new int[][]{{}, {4,5},{4,5},{4,5}});

		return startID;
		
		
		
	}
	
	
	
	public static int fam3(SimulatorStreamPed sim, String tpedPath, int startID) throws IOException{
		
		//first generation: make 2 children
		sim.makeChildren(tpedPath, tpedPath, dir+"gen1.out", startID, startID+1, rGen, 2);
		startID+=2;
		
		
		//first branching family
		//make child
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", 0, startID, rGen, 1);
		startID++;
		
		//make 2 half sibs
		for(int i=0; i<2; i++){
			sim.makeChildren(dir+"gen2.out", tpedPath, String.format(dir+"c%d.out", i), 0, startID, rGen, 1);
			startID++;
		}
		
		
		//second branching family
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"c2.out", 1, startID, rGen, 1);
		
		//another child
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", 1, startID, rGen, 1);
		startID++;
		
		sim.makeChildren(dir+"gen2.out", tpedPath, dir+"c3.out", 0, startID, rGen, 1);
		startID++;
		
		
		//concatenate
		DataParser.concatFilesSpace(new String[]{dir+"c0.out", dir+"c1.out", dir+"c2.out",dir+"c3.out"}, dir+"fam3.out", new int[][]{{}, {4,5},{4,5},{4,5}});
		
		return startID;
		
		
	}
	
	
	
	public static int fam4(SimulatorStreamPed sim, String tpedPath, int startID) throws IOException{
		
		//first generation: make 2 children
		sim.makeChildren(tpedPath, tpedPath, dir+"gen1.out", startID, startID+1, rGen, 2);
		startID+=2;
		
		
		//first branching family
		//make child
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"c0.out", 0, startID, rGen, 1);
		startID++;
		
		//another child
		sim.makeChildren(dir+"c0.out", tpedPath, dir+"gen2.out", 0, startID, rGen, 1);
		startID++;
		
		//c1
		sim.makeChildren(dir+"gen2.out", tpedPath, dir+"c1.out", 0, startID, rGen, 1);
		startID++;

		
		//second branching family
		sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", 1, startID, rGen, 1);
		startID++;
		
		sim.makeChildren(dir+"gen2.out", tpedPath, dir+"c2.out", 0, startID, rGen, 1);
		
		
		//another child
		sim.makeChildren(dir+"gen2.out", tpedPath, dir+"gen3.out", 0, startID, rGen, 1);
		startID++;
		
		sim.makeChildren(dir+"gen3.out", tpedPath, dir+"c3.out", 0, startID, rGen, 1);
		startID++;
		
		
		//concatenate
		DataParser.concatFilesSpace(new String[]{dir+"c0.out", dir+"c1.out", dir+"c2.out",dir+"c3.out"}, dir+"fam4.out", new int[][]{{}, {4,5},{4,5},{4,5}});
		
		return startID;
		
		
	}
	

	
	
	
	public static void sim2(SimulatorStreamPed sim, String tpedPath, int startID, int numUnrel) throws IOException{
		
		//simulate clusters
		int nextID = fam1(sim, tpedPath, startID);
		//System.out.println(nextID);
		nextID = fam2(sim, tpedPath, nextID);
		//System.out.println(nextID);
		nextID = fam3(sim, tpedPath, nextID);
		//System.out.println(nextID);
		nextID = fam4(sim, tpedPath, nextID);
		
		//concatenate clusters and unrel
		//System.out.println(nextID);
		
		String[] files = new String[]{dir+"fam1.out", dir+"fam2.out", dir+"fam3.out",dir+"fam4.out", tpedPath};
		int[][] cols = new int[5][];
		
		//fam1
		cols[0] = new int[]{};
		
		//fam2~fam4
		for(int i=1; i<4; i++){
			cols[i] = new int[]{4,5,6,7,8,9,10,11};		
		}
		
		//unrelated columns
		cols[4] = new int[2*numUnrel];
		for(int i=0; i<numUnrel; i++){
			cols[4][2*i] = 2*(nextID+2);
			cols[4][2*i+1] = 2*(nextID+2) + 1;
			nextID++;
		}
		
		DataParser.concatFilesSpace(files, dir+"sim2.out", cols);
		
		
		
	}
	
	
	
	public static void sim3(SimulatorStreamPed sim, String tpedPath, int startID, int numUnrel) throws IOException{
		
		//first generation: make 2 children
		sim.makeChildren(tpedPath, tpedPath, dir+"gen1.out", startID, startID+1, rGen, 2);
		startID+=2;
		
		
		//for each branching family
		for(int i=0; i<2; i++){
			
			//make a child
			sim.makeChildren(dir+"gen1.out", tpedPath, dir+"gen2.out", i, startID, rGen, 1);
			startID++;
			
			//make a child
			sim.makeChildren(dir+"gen2.out", tpedPath, dir+"gen3.out", 0, startID, rGen, 1);
			startID++;
			
			//make a child
			sim.makeChildren(dir+"gen3.out", tpedPath, dir+"gen4.out", 0, startID, rGen, 1);
			startID++;
			
			DataParser.concatFilesSpace(new String[]{dir+"gen2.out", dir+"gen3.out",dir+"gen4.out"}, String.format(dir+"fam%d.out" ,i), new int[][]{{},{4,5},{4,5}});
			
			
		}
		
		//unrelated columns
		int[] unrel = new int[2*numUnrel];
		for(int i=0; i<numUnrel; i++){
			unrel[2*i] = 2*(startID+2);
			unrel[2*i+1] = 2*(startID+2) + 1;
			startID++;
		}
		
		
		//concatenate final children
		DataParser.concatFilesSpace(new String[]{tpedPath, dir+"gen1.out", dir+"fam0.out", dir+"fam1.out", tpedPath}, dir+"sim3.out", new int[][]{{0,1,2,3,4,5}, {4,5,6,7},{4,5,6,7,8,9},{4,5,6,7,8,9}, unrel});
		
		
		
	}
	
	
	
	

	
	
	
	public static void addError(String tpedPath, String outPath, double errorRate) throws IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(tpedPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		
		String line;

		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			//write header
			writer.write(String.format("%s %s %s %s ", fields[0], fields[1], fields[2], fields[3]));
		
			
			//add error to alleles
			for(int i=4; i<fields.length; i++){

				
				String geno = fields[i];
				
				//if error
				if(Math.random() < errorRate){

					if(geno.equals("A")){
						geno = "T";
					}
					else{
						geno = "A";
					}
				}
				
				//write geno
				writer.write(String.format("%s ", geno));
				
			}
			
			writer.write("\n");
			
		}
		

		reader.close();
		writer.close();

		
	}
	
	
	public static void main(String[] args) throws IOException{
		
		//init simulator
		SimulatorStreamPed sim = new SimulatorStreamPed(recomb);
		double errorRate = .01;
		String outDir = "/Users/kokocakes/Google Drive/Research/pediBear/data/simulations/simPed3/";
		
		
		int[] lengths = new int[]{10,20,30,40};
		
		for(int len : lengths){
			
			//input ped file
			String tpedPath = String.format(dir + "msprime.%dmorgan.0.125.tped", len);
			
			
			//simulate replicate pedigrees
			for(int i=0; i<50; i++){
				
				System.out.println(i);
				
				//simulate
				sim3(sim, tpedPath, 0, 9);

				
				//add error and save
				addError(dir+"sim3.out", String.format(outDir+"sim3.%dmorgan.0.125.%d.tped", len, i), errorRate);
				
				
				
			}
			
			
		}
		


	
	}
	

}

package test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import likelihood.LDStream;
import likelihood.PairwiseLikelihoodCoreStream2;
import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.Relationship;
import utility.DataParser;

public class TestJays {
	
	static double r = 1.5e-8; //zebra finch recombination rate (Backstrom, 2010)
	static double seqError = 0.01;
	static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/";
	static int back = 30000;
	static Random rgen = new Random(526564L);
	
	//build true path
	public static void buildTrue() throws IOException{
		
		String inPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/pedigree.txt";
		String outPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/75jays.true";	
		String indexPath = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/jays/75jays.index2name";	
		int numIndiv = 75;
		
		//get map 
		Map<String, Integer> name2Index = new HashMap<String, Integer>();
		BufferedReader reader = DataParser.openReader(indexPath);
		reader.readLine();
		
		//exclude set
		Set<Integer> exclude = new HashSet<Integer>();
		int[] bad = new int[]{}; 
		for(int i : bad) exclude.add(i);
		
		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			name2Index.put(fields[1], Integer.parseInt(fields[0]));
			
		}
		
		new Pedigree(inPath, outPath, numIndiv, name2Index, exclude);
		
	}
	

	public static void computeMarginals(String simDir, String infoDir, String outPath, int[] indCols, int chrStart, int chrEnd) throws IOException{
		
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
	
	
	public static void computePairwise(String simDir, String infoDir, String outPath, int[] indCols, List<Relationship> relationships, int chrStart, int chrEnd) throws IOException{
		
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

		//relationships
		List<Relationship> relationships = new ArrayList<Relationship>();

		relationships.add(new Relationship(15d, new double[] {1-1d/2048d, 1d/2048d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
		// depth = 1 relationship
		relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
		relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
		relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
		
		//depth = 2 relationships
		relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
		relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
		relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
		relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
		
		//depth = 3 relationships
		relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1), new Path(5,0,1)})); 
		relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1), new Path(5,1,1), new Path(6,0,1)}));
		relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
		relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2), new Path(5,1,2)})); 
		
		//depth = 4 relationships 
		relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1), new Path(5,2,1), new Path(6,1,1)}));
		relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1), new Path(5,3,1), new Path(6,2,1)})); 
		relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2), new Path(5,2,2), new Path(6,1,2)}));
		relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2), new Path(5,3,2), new Path(6,2,2)}));
		
		
		//depth = 5 relationships
		relationships.add(new Relationship(9d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,4,1), new Path(6,3,1)}));
		relationships.add(new Relationship(10d, new double[] {511d/512d, 1d/512d, 0d}, new Path[]{new Path(5,5,1), new Path(6,4,1)}));
		relationships.add(new Relationship(11d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(5,4,2), new Path(6,3,2)}));
		relationships.add(new Relationship(12d, new double[] {255d/256d, 1d/256d, 0d}, new Path[]{new Path(5,5,2), new Path(6,4,2)}));
		
		//depth = 6 relationships
		relationships.add(new Relationship(11d, new double[] {1023d/1024d, 1d/1024d, 0d}, new Path[]{new Path(6,5,1)}));
		relationships.add(new Relationship(12d, new double[] {2047d/2048d, 1d/2048d, 0d}, new Path[]{new Path(6,6,1)}));
		relationships.add(new Relationship(13d, new double[] {511d/512d, 1d/512d, 0d}, new Path[]{new Path(6,5,2)}));
		relationships.add(new Relationship(14d, new double[] {1023d/1024d, 1d/1024d, 0d}, new Path[]{new Path(6,6,2)}));
		
		//parameters
		int chrStart = 1;
		int chrEnd = 35;
		
		//individuals
		int numIndiv = 75;
		int[] indCols = new int[numIndiv];
		for(int i=0; i<numIndiv; i++) indCols[i] = i+3;	
		
		
		//compute info
		System.out.println("Computing info");
		for (int i=chrStart; i<chrEnd; i++){
			String inPath = dir+"75jays.geno."+i;
			String outPath = dir+"75jays.info."+i;
			LDStream.writeLdOutfile(inPath, outPath, back);
		}
		
		
		//compute marginal
		System.out.println("Computing marginals");
		computeMarginals(dir + "75jays.geno.", dir+"75jays.info.", dir+"75jays.marginal", indCols, chrStart, chrEnd);
		
		
		//compute pairwise
		System.out.println("Computing pairwise likelihoods");
		long startTime = System.nanoTime();		
		computePairwise(dir+"75jays.geno.", dir+"75jays.info.", dir+"75jays.pairwise", indCols, relationships, chrStart, chrEnd);	
		System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");

	
		
	
	}

	
	
}

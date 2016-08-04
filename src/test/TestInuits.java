package test;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import likelihood.LDStreamPed;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import dataStructures.Path;
import dataStructures.Relationship;
import utility.DataParser;

public class TestInuits {
	
	static double recombRate = 1.3e-8; //zebra finch recombination rate (Backstrom, 2010)
	static double seqError = 0.01;
	static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/inuits/";
	static int back = 500000;
	static Random rgen = new Random(526564L);


	public static void computeMarginals(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids) throws IOException{
		
		//compute
		double[] marginals = core.computeMarginal(fileName+".tped", fileName+".info", ids);

		//open outfile
		PrintWriter writer = DataParser.openWriter(fileName+".marginal");				
		
		//write marginals to file
		for (double item : marginals){
			writer.write(String.format("%f\n",item));
		}
		writer.close();
	}
	
	
	public static void computePairwise(PairwiseLikelihoodCoreStreamPed core, String fileName, int[] ids, List<Relationship> relationships) throws IOException{
		
		int numIndiv = ids.length;
		int numRel = relationships.size();
		
		//likelihood
		PrintWriter writer = DataParser.openWriter(fileName+".pairwise");
			
		double[][][] lkhd = core.forwardAlgorithm(fileName+".tped", fileName+".info", ids, relationships);
		

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

		relationships.add(new Relationship(15d, new double[] {1d, 0d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
		
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
		 
		
		
		//individuals
		int numIndiv = 100;
		int[] indCols = new int[numIndiv];
		for(int i=0; i<numIndiv; i++) indCols[i] = i;
		
		//pairwise core
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(seqError, recombRate, back, numIndiv);
		
		String fileName = dir + "100tasiilaq.admixed0.5.aims0.05.prune0.1";

		
		//compute info
		System.out.println("Computing info");
		LDStreamPed.writeLdOutfile(fileName, fileName + ".info", back);
		
		
		//compute marginal
		System.out.println("Computing marginals");
		computeMarginals(core, fileName, indCols);
		core.setMarginals(fileName+".marginal");
		 
		
		//compute pairwise
		System.out.println("Computing pairwise likelihoods");
		long startTime = System.nanoTime();		
		computePairwise(core, fileName, indCols, relationships);	
		System.out.println("Total time was " + (System.nanoTime() - startTime) / 1e9 / 60d/ 60d + " hours");
		
		
		
	
	}

	
	
}
package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.io.PrintWriter;


import likelihood.PairwiseLikelihoodCoreMicrosatellites;
import dataStructures.Relationship;
import dataStructures.Path;

import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class TestFrogs{
				
		
		//margina for microsats
		public static void computeMarginalsMicrosat(PairwiseLikelihoodCoreMicrosatellites core, String fileName, String freqPath, int[] ids, int t) throws IOException{
			
			//compute
			double[] marginals = core.computeMarginalMicrosat(fileName+".tped", freqPath, ids);
			
			
			//open outfile
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.marginal", fileName, t));				
			
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
			PrintWriter writer = DataParser.openWriter(String.format("%s.%d.pairwise", fileName, t));
				
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
			*/
			
			//////////////////////////////
	
			
			//parameters 
			String dir = "/Users/amy/eclipse-workspace/mcmc/frogs/";
			int sampleSize = 90;
			double epsilon1 = .005;
			double epsilon2 = .01;

			
			//individuals
			int[] indCols = new int[sampleSize];
			for(int i=0; i<sampleSize; i++) indCols[i] = i;
			

			
			String testName = "frogs.juv";
			
			
			//pairwise core
			PairwiseLikelihoodCoreMicrosatellites core = new PairwiseLikelihoodCoreMicrosatellites(sampleSize, epsilon1, epsilon2);
	
			
			//compute marginal probs
			System.out.println("Computing marginals");
			computeMarginalsMicrosat(core, dir + testName, dir + "frogs.all.freq", indCols, 0);
			
			
			//compute pairwise
			System.out.println("Computing pairwise likelihoods");		
			computePairwiseMicrosat(core, dir + testName, dir + "frogs.all.freq", indCols, relationships, 0);	

		}


	
}
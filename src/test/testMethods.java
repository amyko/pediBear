package test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import statistic.Accuracy;
import Unused.PairwiseLikelihoodCoreStream;
import likelihood.PairwiseLikelihoodCoreStream2;
import dataStructures.Node;
import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.Relationship;


public class testMethods {
	
	
	
	public static void main(String[] args) throws IOException{

		
		//parameters
		int depth = 4;
		int numIndiv = 6;
		double seqError = 0.01;
		double r = 1.3e-8;
		double genTime = 16;
		int back = 100;
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/pairwiseLikelihood/test5.";
		Random rGen = new Random(1489864090);
		int numMaxNodes = 2;
		PairwiseLikelihoodCoreStream2 core2 = new PairwiseLikelihoodCoreStream2(seqError, r, back, numIndiv);
		String marginalPath = dir + ".marginal."+0;
		String lkhdPath = dir + ".pairwise."+0;

		
		//core2.setMarginals(marginalPath);
		//core2.setLikelihoods(lkhdPath);
		
		
		//relationships
		List<Relationship> relationships = new ArrayList<Relationship>();

		
		relationships.add(new Relationship(11d, new double[] {1-1d/128d, 1d/128d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
		// depth = 1 relationship
		relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
		relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
		relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib
		
		//depth = 2 relationships
		//relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, 2,0,1)); //grand parents
		relationships.add(new Relationship(3d, new double[] {.75, .25, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
		relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
		relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
		relationships.add(new Relationship(6d, new double[] {.75, .25, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
		
		//depth = 3 relationships
		//relationships.add(new Relationship(3d, new double[] {.75, .25, 0d}, 3,0,1)); 
		//relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, 3,1,1)); 
		relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1)})); 
		relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1)}));
		//relationships.add(new Relationship(6d, new double[] {.75, .25, 0d}, 3,1,2));
		relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
		relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2)})); 
		
		//depth = 4 relationships
		//relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, 4,0,1)); 
		//relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, 4,1,1)); 
		//relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, 4,2,1)); 
		relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1)}));
		relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1)}));
		//relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, 4,1,2)); 
		//relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, 4,2,2)); 
		relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2)}));
		relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2)}));
		
		
		int c=0;
		
		for(int t=7; t<8; t++){
			
			lkhdPath = dir +"pairwise."+t;
			core2.setLikelihoods(lkhdPath);
			

			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Relationship bestRel = null;
					Relationship secondBest = bestRel;
					double bestLkhd = Double.NEGATIVE_INFINITY;
					double secondLkhd = bestLkhd;
					
					for(Relationship rel : relationships){
						
						double currLkhd = core2.getLikelihood(new Node(true, i), new Node(true, j), rel.getPath());
						
						if(currLkhd > bestLkhd){
							secondBest = bestRel;
							secondLkhd = bestLkhd;
							bestLkhd = currLkhd;
							bestRel = rel;
						}
						else if(currLkhd > secondLkhd){
							secondBest = rel;
							secondLkhd = currLkhd;
						}
						
					}
					
					//most likely relationship
					//System.out.println(String.format("(%d, %d), (%d, %d, %d)", i,j, bestRel.getPath().getUp(), bestRel.getPath().getDown(), bestRel.getPath().getNumVisit()));
					
					
					if(!bestRel.getPath().equals(new Path(0,0,0))){
						System.out.println(String.format("(%d, %d), (%d, %d, %d)", i,j, bestRel.getPath().getUp(), bestRel.getPath().getDown(), bestRel.getPath().getNumVisit()));
						System.out.println(String.format("(%d, %d), (%d, %d, %d)", i,j, secondBest.getPath().getUp(), secondBest.getPath().getDown(), secondBest.getPath().getNumVisit()));
						System.out.println();
						c++;
					}

					
					
				}
			
			}
			
		}
		
		System.out.println(c);

		





		

		
	}
	
	
}

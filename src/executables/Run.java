package executables;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

import dataStructures.Path;
import dataStructures.Relationship;
import likelihood.LDStreamPedMissing;
import likelihood.PairwiseLikelihoodCoreStreamPed;
import likelihood.PreProcess;
import utility.DataParser;

public class Run {
	
	//private static String[] validOptions = new String[]{ "iterpertemp","tempfact","maxiter","conv"};

	//SA parameters
	public static String fileName = "";
	public static String refPopFileName = "";
	public static String ageFileName = "";
	public static double maf = 0.01;
	public static double errorRate = 0.01;
	public static int maxDepth = 4;
	public static int sampleDepth = maxDepth;
	public static double back = 0.04;
	public static double startTemp = 100;
	public static double tempFact = 1.01;
	public static int iterPerTemp = 20000;
	public static int maxIter = 4000000;
	public static double conv = 1;
	public static int numIndiv = 0;
	public static boolean conditional = true;
	public static int numRun = 3;
	
	//misc
	public static int maxNumNodes = 200;
	
	
	//relationships for likelihood computation
	private static List<Relationship> relationships = new ArrayList<Relationship>();


	public static void computeLikelihoods() throws IOException{
		
		//init relationships
		initRelationships();
		
		//compute info
		LDStreamPedMissing.writeLdOutfile(refPopFileName, refPopFileName + ".info", back, conditional);
		
		//init core
		PairwiseLikelihoodCoreStreamPed core = new PairwiseLikelihoodCoreStreamPed(errorRate, back, numIndiv);
		
		
		
		//compute marginal likelihood
		System.out.println("Computing marginal likelihoods");
		int[] indCols = new int[numIndiv];
		for(int i=0; i<numIndiv; i++) indCols[i] = i;
		computeMarginals(core, fileName, indCols);
		
		//compute pairwise likelihood
		System.out.println("Computing pairwise likelihoods");
		computePairwise(core, fileName, indCols, relationships);
		
		
	}
	
	
	
	public static void initRelationships(){
		
		relationships.add(new Relationship(15d, new double[] {1d, 0d, 0}, new Path[]{new Path(0,0,0)})); //unrelated
		
		
		// depth = 1 relationship
		relationships.add(new Relationship(1d, new double[] {0d, 1d, 0d}, new Path[]{new Path(1,0,1)})); //parent
		relationships.add(new Relationship(4d, new double[] {.25, .5, .25}, new Path[]{new Path(1,1,2)})); //full-sib	
		relationships.add(new Relationship(2d, new double[] {.5, .5, 0d}, new Path[]{new Path(1,1,1), new Path(2,0,1)})); //half-sib
		
		//depth = 2 relationships
		if(maxDepth >= 2){
			relationships.add(new Relationship(3d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,1,1), new Path(3,0,1)})); //half uncle
			relationships.add(new Relationship(4d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(2,2,1), new Path(3,1,1), new Path(4,0,1)})); //half cousins
			relationships.add(new Relationship(5d, new double[] {.5, .5, 0d}, new Path[]{new Path(2,1,2)})); //uncle
			relationships.add(new Relationship(6d, new double[] {3d/4d, 1d/4d, 0d}, new Path[]{new Path(2,2,2), new Path(3,1,2)})); //first cousins
		}
		
		//depth = 3 relationships
		if(maxDepth >= 3){
			relationships.add(new Relationship(5d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,2,1), new Path(4,1,1)})); 
			relationships.add(new Relationship(6d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(3,3,1), new Path(4,2,1)}));
			relationships.add(new Relationship(7d, new double[] {7d/8d, 1d/8d, 0d}, new Path[]{new Path(3,2,2), new Path(4,1,2)})); 
			relationships.add(new Relationship(8d, new double[] {15d/16d, 1d/16d, 0d}, new Path[]{new Path(3,3,2), new Path(4,2,2)})); 
		}
		
		//depth = 4 relationships 
		if(maxDepth >= 4){
			relationships.add(new Relationship(7d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,3,1)}));
			relationships.add(new Relationship(8d, new double[] {127d/128d, 1d/128d, 0d}, new Path[]{new Path(4,4,1)})); 
			relationships.add(new Relationship(9d, new double[] {31d/32d, 1d/32d, 0d}, new Path[]{new Path(4,3,2)}));
			relationships.add(new Relationship(10d, new double[] {63d/64d, 1d/64d, 0d}, new Path[]{new Path(4,4,2)}));
		
		}
		
	}
	

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
	
		System.out.println("Preprocessing...");
		PreProcess.processOptionfile(args);
		PreProcess.checkInputFiles(fileName, refPopFileName);
		
		//TODO use plink to filter maf & geno and prune
		
		
		//computeLikelihoods();
		
		//annealing schedule
		double[] heat = new double[maxIter / iterPerTemp];
		heat[0] = 1d / startTemp;
		for(int i=1; i<heat.length; i++) heat[i] = heat[i-1]*tempFact;
		
		int[] coolingSchedule = new int[heat.length-1];
		for(int i=0; i<heat.length-1; i++){
			coolingSchedule[i] = iterPerTemp;
		}
		
		
		
		
		
		
		
	}
	

}

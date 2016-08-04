package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Writer;
import java.util.Map;

import statistic.Accuracy;
import utility.DataParser;
import dataStructures.Path;


//infers the relationship from pairwise likelihood
public class RunPairwiseTest {

	//get path with the highest likelihood
	public static Path getHighestLkhdPath(String lkhdPath, int i, int j) throws NumberFormatException, IOException{

		
		//open reader
		BufferedReader reader = DataParser.openReader(lkhdPath);
		
		Path bestPath = null;
		double bestLkhd = Double.NEGATIVE_INFINITY;
		Path currPath = null;
		double currLkhd = Double.NEGATIVE_INFINITY;
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//compare relationship
			if(fields[0].equals(">")){
				currPath = new Path(Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), Integer.parseInt(fields[3]));
				//System.out.print(String.format("%d\t%d\t%d\n", currPath.getUp(), currPath.getDown(), currPath.getNumVisit()));
				continue;
			}
			
			
			int currI = Integer.parseInt(fields[0]);
			int currJ = Integer.parseInt(fields[1]);
			
			//find i and j
			if(currI==i && currJ==j){

				currLkhd = Double.parseDouble(fields[2]);
				
				//System.out.println(currLkhd);

				
				if(currLkhd > bestLkhd){
					bestLkhd = currLkhd;
					bestPath = currPath;
				}
				
			}


			
		}

		reader.close();
		
		return bestPath;
					
		
	}
	
	
	
	//get path with the highest likelihood
	public static Path[][] getHighestLkhdPath(String lkhdPath, int numIndiv) throws NumberFormatException, IOException{

		double[][] bestLkhd = new double[numIndiv][numIndiv];
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				bestLkhd[i][j] = Double.NEGATIVE_INFINITY;
			}
		}
		Path[][] bestPath = new Path[numIndiv][numIndiv];
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				bestPath[i][j] = null;
			}
		}
		
		//open reader
		BufferedReader reader = DataParser.openReader(lkhdPath);

		Path currPath = null;
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//compare relationship
			if(fields[0].equals(">")){
				currPath = new Path(Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), Integer.parseInt(fields[3]));
				continue;
			}
			
			
			int i = Integer.parseInt(fields[0]);
			int j = Integer.parseInt(fields[1]);
			double currLkhd = Double.parseDouble(fields[2]);
			

			
			if(currLkhd > bestLkhd[i][j]){
				
				bestLkhd[i][j] = currLkhd;
				bestPath[i][j] = currPath;
				
			}
			
		}

		reader.close();
		
		return bestPath;
					
		
	}
	
	
	
	
	public static void main(String[] args) throws IOException{
		
		//param
		int numIndiv = 100;
		
		//files
		/* for simulations
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String testName = "test13";
		String inPath = dir + "pairwiseLikelihood/"+testName+".pairwise.";
		String outPath = dir + "results/"+testName+".out";
		String mapAccPath = dir + "results/"+testName+".pairwise.map.acc";
		String pathToOmegaPath = dir + "pathToOmega.txt";
		//String truePath = dir + testName + ".true";
		String truePath = dir + "results/test12.true";
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmegaPath);
		*/ 
		
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/inuits/";
		String testName = dir + "100tasiilaq.admixed0.5.aims0.05.prune0.1";
		String inPath = testName+".pairwise";
		String outPath = testName+".out";


		//open outfiles
		//Writer mapWriter = DataParser.openWriter(mapAccPath);
		
		//run pairwise test
		for(int t=0; t<1; t++){
			
			//open writer
			Writer writer2 = DataParser.openWriter(outPath);
			writer2.write(String.format(">\t%d\n", t));
			
			//write to acc files
			//mapWriter.write(String.format(">\t%d\n", t));
			
			//String lkhdPath = inPath + t;
			String lkhdPath = inPath;

			
			
			//pairwise result
			Path[][] bestPaths = getHighestLkhdPath(lkhdPath, numIndiv);
			
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					//Path bestPath = getHighestLkhdPath(lkhdPath, i, j);
					Path bestPath = bestPaths[i][j];
					
					writer2.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));	
				}
			}
			
			writer2.close();
			
			/*
			//accuracy based on mcmc
			double[][] mapAcc = Accuracy.kinshipAccuracy(outPath, truePath, numIndiv, pathToKinship);

			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					mapWriter.write(String.format("%d\t%d\t%.3f\n", i, j, mapAcc[i][j]));
				}
			}
			*/
			
			
		}
		
		//mapWriter.close();
		
		
		
		
	}
	
	
}

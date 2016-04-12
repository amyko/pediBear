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

	
	public static Path getHighestLkhdPath(String lkhdPath, int i, int j) throws NumberFormatException, IOException{

		
		//open reader
		BufferedReader reader = DataParser.openReader(lkhdPath);
		
		Path bestPath = null;
		double bestLkhd = Double.NEGATIVE_INFINITY;
		Path currPath = null;
		double currLkhd = Double.NEGATIVE_INFINITY;
		
		String line;
		int currLineNum = 0;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			//compare relationship
			if(fields[0].equals(">")){
				currPath = new Path(Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), Integer.parseInt(fields[3]));
				currLineNum = 0;
				continue;
			}
			
			
			//find i and j
			if(currLineNum==i){

				currLkhd = Double.parseDouble(fields[j]);
				
				if(currLkhd > bestLkhd){
					bestLkhd = currLkhd;
					bestPath = currPath;
				}
				
			}

			currLineNum++;

			
		}

		reader.close();
		
		return bestPath;
					
		
	}
	
	
	
	
	
	public static void main(String[] args) throws IOException{
		
		//param
		int numIndiv = 6;
		
		//files
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String testName = "test5";
		String inPath = dir + "pairwiseLikelihood/"+testName+".pairwise.";
		String outPath = dir + "results/"+testName+".out";
		String accPath = dir + "results/"+testName+".kinship.pairwise.acc";
		String pathToOmegaPath = dir + "pathToOmega.txt";
		String truePath = dir + "results/" +testName + ".true";
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmegaPath);
		
		//open outfiles
		Writer writer = DataParser.openWriter(accPath);
		
		//run pairwise test
		for(int t=0; t<100; t++){
			
			Writer writer2 = DataParser.openWriter(outPath);
			writer2.write(String.format(">\t%d\n", t));
			writer.write(String.format(">\t%d\n", t));
			
			String lkhdPath = inPath + t;
			
			
			//pairwise result
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					
					Path bestPath = getHighestLkhdPath(lkhdPath, i, j);
					
					writer2.write(String.format("%d\t%d\t%d\t%d\t%d\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
					
					
				}
			}
			
			writer2.close();
			
			//accuracy
			double[][] kinshipAcc = Accuracy.kinshipAccuracy(outPath, truePath, numIndiv, pathToKinship);
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					writer.write(String.format("%d\t%d\t%.3f\n", i, j, kinshipAcc[i][j]));
				}
			}
			
			
		}
		
		
		
		writer.close();

		
		
		
	}
	
	
}

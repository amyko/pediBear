package statistic;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import utility.DataParser;
import dataStructures.Path;


public class Accuracy {

	//output file format
	final static int IND1 = 0;
	final static int IND2 = 1;
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	
	
	//ordered = finds the number of samples with the exact path, unordered = exact path + complement
	public static int numSamplesWithGivenRel(String outPath, int i, int j, Path path, boolean ordered) throws IOException{
		
		BufferedReader reader = DataParser.openReader(outPath);
		
		int n = 0;
		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")) continue;
			
			if(Integer.parseInt(fields[IND1])!=i || Integer.parseInt(fields[IND2])!=j) continue;
			

			if(ordered){
				if(Integer.parseInt(fields[NUMVISIT])==path.getNumVisit()){
					if(Integer.parseInt(fields[UP])==path.getUp() && Integer.parseInt(fields[DOWN])==path.getDown()){
						n++;				
					}
				}
			}
			
			else{
				if(Integer.parseInt(fields[NUMVISIT])==path.getNumVisit()){
					if(Integer.parseInt(fields[UP])==path.getUp() && Integer.parseInt(fields[DOWN])==path.getDown()){
						n++;				
					}
					
					if(Integer.parseInt(fields[UP])==path.getDown() && Integer.parseInt(fields[DOWN])==path.getUp()){
						n++;	
					}
				}	
			}
			

			
		}
		
		reader.close();
		
		return n;
		
	}
	
	
	//returns path with the highest number of hits
	public static void mostLikelyPath(String outPath, int i, int j) throws NumberFormatException, IOException{
		
		Map<Path, Integer> count = new HashMap<Path, Integer>();
		
		BufferedReader reader = DataParser.openReader(outPath);
		

		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")) continue;
			
			if(Integer.parseInt(fields[0])!=i || Integer.parseInt(fields[1])!=j) continue;
			
			
			Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
			
			//increment count
			if(count.containsKey(key)){
				count.put(key, count.get(key) + 1);
			}
			else{
				count.put(key, 1);
			}
			
			
		}
		
		reader.close();
		
		
		//determine most likely path
		Path bestPath = null;
		int bestCount = 0;
		int currCount = 0;
		
		for(Path key : count.keySet()){
			
			currCount = count.get(key);
			if(currCount > bestCount){
				bestCount = currCount;
				bestPath = key;
			}
			
		}
		
		//print
		System.out.println(String.format("(%d, %d) : (%d, %d, %d) %d\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit(), bestCount));
		
	}
	
	
	public static double[][] kinshipAccuracy(String outPath, String truePath, int numIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		double[][] kinshipDist = new double[numIndiv][numIndiv];
		
		//read true relationship
		double[][][] trueOmega = getTrueOmega(truePath, numIndiv, pathToOmega);
			
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		String line;
		int numSample = 0;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				numSample++;
				continue;
			}
			
			int i = Integer.parseInt(fields[IND1]);
			int j = Integer.parseInt(fields[IND2]);
			Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
			
			int c = 0;
			for(int k=0; k<3; k++){
				if(Math.abs(trueOmega[i][j][k] - pathToOmega.get(key)[k]) < 1e-13){
					c++;
				}
			}
			
			int add = c==3 ? 1 : 0;
			
			kinshipDist[i][j] += add;
			
		}
		
		reader.close();
		
		
		//normalize by sample size
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				kinshipDist[i][j] /= numSample;
			}
		}
		
		return kinshipDist;
		
		
	}
	
	public static double[][] relAccuracy(String outPath, String truePath, int numIndiv) throws NumberFormatException, IOException{
		
		double[][] accuracy = new double[numIndiv][numIndiv];
		
		//read true paths
		Path[][] trueRel = getTruePath(truePath, numIndiv);
		
		
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		String line;
		int numSample = 0;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				numSample++;
				continue;
			}
			
			int i = Integer.parseInt(fields[IND1]);
			int j = Integer.parseInt(fields[IND2]);
			Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
			
			if(trueRel[i][j].equals(key)){
				accuracy[i][j] += 1;
			}

		}
		
		reader.close();
		
		
		//normalize by sample size
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				accuracy[i][j] /= numSample;
			}
		}
		
		
		return accuracy;
		
	}
	
	
	
	public static Map<Path, double[]> getPathToOmega(String pathToOmegaPath) throws NumberFormatException, IOException{
		
		Map<Path, double[]> pathToKinship = new HashMap<Path, double[]>();
		
		
		BufferedReader reader = DataParser.openReader(pathToOmegaPath);
		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			Path key = new Path(Integer.parseInt(fields[0]) , Integer.parseInt(fields[1]) , Integer.parseInt(fields[2]));
			
			//omega 0
			String[] frac = fields[4].split("/");
			double omega0 = Double.parseDouble(frac[0]) / Double.parseDouble(frac[1]);
			
			//omega 1
			frac = fields[5].split("/");
			double omega1 = Double.parseDouble(frac[0]) / Double.parseDouble(frac[1]);
			
			//omega 2
			frac = fields[6].split("/");
			double omega2 = Double.parseDouble(frac[0]) / Double.parseDouble(frac[1]);

			//kinship coeff
			double[] val = new double[]{omega0, omega1, omega2};
			
			pathToKinship.put(key, val);
		
		}
		
		reader.close();
		

		return pathToKinship;
	}
	
	

	public static double[][][] getTrueOmega(String inPath, int numIndiv, Map<Path, double[]> pathToKinship) throws NumberFormatException, IOException{
		
		//read true relationship
		double[][][] trueKinship = new double[numIndiv][numIndiv][3];
		
		BufferedReader reader = DataParser.openReader(inPath);
		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			int i = Integer.parseInt(fields[IND1]);
			int j = Integer.parseInt(fields[IND2]);
			Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
			
			for(int k=0; k<3; k++)
				trueKinship[i][j][k] = pathToKinship.get(key)[k];
			
		}
		
		reader.close();
		
		return trueKinship;
		
	}
	
	
	public static Path[][] getTruePath(String inPath, int numIndiv) throws NumberFormatException, IOException{
		
		//read true relationship
		Path[][] trueRel = new Path[numIndiv][numIndiv];
		
		BufferedReader reader = DataParser.openReader(inPath);
		String line;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			int i = Integer.parseInt(fields[IND1]);
			int j = Integer.parseInt(fields[IND2]);
			Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
			
			trueRel[i][j] = key;
			
		}
		
		reader.close();
		
		return trueRel;
		
	}
	
	
}

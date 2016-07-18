package statistic;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
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
	
	
	//best paths based on pairwise test
	public static Path[][] pairwiseBestPed(String outPath, int totalIndiv, int numIndiv) throws NumberFormatException, IOException{
		
		//to return
		Path[][] bestPath = new Path[numIndiv][numIndiv];
		double[][] bestLkhd = new double[numIndiv][numIndiv];
		
		for(int i=0; i<numIndiv; i++){
			for(int j=0; j<numIndiv; j++){
				bestLkhd[i][j] = Double.NEGATIVE_INFINITY;
			}
		}
		
		
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		String line;
		Path rel = null;
		int i=0;
		
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				rel = new Path(Integer.parseInt(fields[1]) , Integer.parseInt(fields[2]) , Integer.parseInt(fields[3]));
				i=0;
				continue;
			}
			
			if(i>=numIndiv) continue;
			
			for(int j=0; j<numIndiv; j++){
				double currLkhd = Double.parseDouble(fields[j]);
				if(currLkhd > bestLkhd[i][j]){
					bestLkhd[i][j] = currLkhd;
					bestPath[i][j] = rel;
				}
			}
			
			i++;
			
			
		}
		
		reader.close();
		
		
		
		return bestPath;
		
		
		
	}
	
	
	
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
					
					else if(Integer.parseInt(fields[UP])==path.getDown() && Integer.parseInt(fields[DOWN])==path.getUp()){
						n++;	
					}
				}	
			}
			

			
		}
		
		reader.close();
		
		return n;
		
	}
	
	
	//ordered = finds the number of samples with the exact path, unordered = exact path + complement
	public static int numSamplesWithGivenRel(String outPath, Path path) throws IOException{
		
		BufferedReader reader = DataParser.openReader(outPath);
		
		int n = 0;
		String line;
		int c = 0;
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				if(c==3) n++;
				c = 0;
				continue;
			}
			
			
			if(Integer.parseInt(fields[NUMVISIT])==path.getNumVisit()){
				if(Integer.parseInt(fields[UP])==path.getUp() && Integer.parseInt(fields[DOWN])==path.getDown()){
					c++;				
				}
				
				else if(Integer.parseInt(fields[UP])==path.getDown() && Integer.parseInt(fields[DOWN])==path.getUp()){
					c++;	
				}
			}
			
			

			
		}
		
		reader.close();
		
		return n;
		
	}
	
	
	//returns path with the highest number of hits
	public static Path mostLikelyPath(String outPath, int i, int j) throws NumberFormatException, IOException{
		
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
		
		return bestPath;
		
		
	}
	
	
	
	//accuracy based on MAP estimate
	public static double[][] mapAccuracy(String outPath, String truePath, int totalIndiv, int numIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		double[][] accuracy = new double[numIndiv][numIndiv];
		Path[][] mapPed = new Path[numIndiv][numIndiv];
		
		//read true relationship
		double[][][] trueOmega = getTrueOmega(truePath, totalIndiv, pathToOmega);
			
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		String line;
		double bestLkhd = Double.NEGATIVE_INFINITY;
		boolean process = false;
		
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				
				double currLkhd = Double.parseDouble(fields[1]);
				
				if(currLkhd > bestLkhd){
					bestLkhd = currLkhd;
					process = true;
				}
				else
					process = false;
				
				continue;
			}
			
			
			if(process){
				
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
				
				accuracy[i][j] = add;
				
				//best path
				mapPed[i][j] = key;
				
			}

			

			
		}
		
		reader.close();
		
		//print map path
		System.out.println();
		System.out.println(String.format("MCMC MAP likelihood: %.2f", bestLkhd));
				
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				Path bestPath = mapPed[i][j];
				//System.out.println(String.format("(%d, %d) : (%d, %d, %d)\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
			}
		}
		
		
		return accuracy;
		
		
	}
	
	
	public static void printMAP(String outPath, int numIndiv) throws NumberFormatException, IOException{
		

		Path[][] mapPed = new Path[numIndiv][numIndiv];

		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		String line;
		double bestLkhd = Double.NEGATIVE_INFINITY;
		boolean process = false;
		
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				
				double currLkhd = Double.parseDouble(fields[1]);
				
				if(currLkhd > bestLkhd){
					bestLkhd = currLkhd;
					process = true;
				}
				else
					process = false;
				
				continue;
			}
			
			
			if(process){
				
				int i = Integer.parseInt(fields[IND1]);
				int j = Integer.parseInt(fields[IND2]);
				Path key = new Path(Integer.parseInt(fields[UP]) , Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[NUMVISIT]));
				
				//best path
				mapPed[i][j] = key;
				
			}

			

			
		}
		
		reader.close();
		
		//print map path
		System.out.println();
		System.out.println(String.format("MCMC MAP likelihood: %.2f", bestLkhd));
				
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				Path bestPath = mapPed[i][j];
				System.out.println(String.format("(%d, %d) : (%d, %d, %d)\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
			}
		}

	}
	
	
	
	//acc=1 if inferred ibd coefficeints are the same as truth
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
	
	
	//average r = .25*k1 + .5k2 across posterior samples
	public static double[][] ibdAccuracy(String outPath, int numIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		double[][] kinshipDist = new double[numIndiv][numIndiv];
		
			
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
			
			kinshipDist[i][j] += .25*pathToOmega.get(key)[1] + .5*pathToOmega.get(key)[2];
			
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
			
			Path key2 = new Path(Integer.parseInt(fields[DOWN]) , Integer.parseInt(fields[UP]) , Integer.parseInt(fields[NUMVISIT]));
			trueRel[j][i] = key2;
			
		}
		
		reader.close();
		
		return trueRel;
		
	}
	
	
	public static boolean hasSameIBD(Map<Path, double[]> pathToOmega, Path i, Path j){
	
		double[] omega1 = pathToOmega.get(i);
		double[] omega2 = pathToOmega.get(j);
		
		if(omega1[0]==omega2[0] && omega1[1]==omega2[1] && omega1[2]==omega2[2])
			return true;
		
		return false;
		
	}
	
	
}

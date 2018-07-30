package statistic;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;

import utility.DataParser;
import dataStructures.Path;


public class Accuracy {

	//output file format
	final static int IND1 = 0;
	final static int IND2 = 1;
	final static int UP = 2;
	final static int DOWN = 3;
	final static int NUMVISIT = 4;
	

	//returns (1): true if MAP is true; (2) #MAP_samples / #total
	public static String[] inMAP(String countPath, String truePed) throws IOException{
			
		//get MAP
		int bestCount = 0;
		String bestPed = "";
		
		
		BufferedReader reader = DataParser.openReader(countPath);
		
		String line;
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			int currCount = Integer.parseInt(fields[1]);
			
			if(currCount > bestCount){
				bestPed = fields[0];
				bestCount = currCount;
			}
			
			
		}
		reader.close();
	
		
		//TODO PRINT
		//print best map
		//System.out.println(bestPed);
		
		return null;
	}
	
	
	//returns (1): true if true pedigree is in the credible interval; (2) #pedigrees in credible interval
	public static String[] inCI(String countPath, String truePed, double CI, double denom) throws NumberFormatException, IOException{
	
		//build map: count -> ped
		Map<Integer, String> count2ped = new HashMap<Integer, String>();
		
		BufferedReader reader = DataParser.openReader(countPath);
		String line;
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			int key = -Integer.parseInt(fields[1]);
			String val = fields[0];
			count2ped.put(key, val);
			
		}
		reader.close();
		
		
		//sort map
		Map<Integer, String> sortedMap = new TreeMap<Integer, String>(count2ped); 
		
		//credible interval
		int totalCount = 0;
		int pedCount = 0;
		int inCI = 0;
		int inCIExceptFalsePositives = 0;
		
		for(int key : sortedMap.keySet()){
			
			//get count
			int count = Math.abs(key);
			totalCount += count;
			pedCount++;
			
			//check if this pedigree is the true one
			String testPed = sortedMap.get(key);
			if(inCI==0 && truePed.equals(testPed)){
				inCI=1;
			}
			
			//check if this pedigree is the true one ignoring the false positives
			if(inCIExceptFalsePositives==0) {
			
				int same = 1;
				
				for(int i=0; i<testPed.length()/3; i++) {
					
					String testPair = "" + testPed.charAt(3*i) + testPed.charAt(3*i+1) + testPed.charAt(3*i + 2);
					String truePair = "" + truePed.charAt(3*i) + truePed.charAt(3*i+1) + truePed.charAt(3*i + 2);
					
					
					if(! testPair.equals(truePair) && ! truePair.equals("000")) {
						
						same *= 0;
						break;
						
					}
					
				}
				
				if(same==1)
					inCIExceptFalsePositives = 1;
			
			}
			
			if(totalCount / denom > CI) break;
			
		}
		

		//record results
		String[] toReturn = new String[3];
		toReturn[0] = String.format("%d", inCI);
		toReturn[1] = String.format("%d", inCIExceptFalsePositives);
		toReturn[2] = String.format("%d", pedCount);

		
		
		return toReturn;	
		
		
	}
	
	//returns (1): true if true pedigree is in the credible interval; (2) #pedigrees in credible interval
	public static String[] inCIPair(String countPath, String truePed, double CI, double denom) throws NumberFormatException, IOException{
	
		//build map: count -> ped
		Map<Integer, String> count2ped = new HashMap<Integer, String>();
		
		BufferedReader reader = DataParser.openReader(countPath);
		String line;
		
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\\s");
			
			int key = -Integer.parseInt(fields[1]);
			String val = fields[0];
			count2ped.put(key, val);
			
		}
		reader.close();
		
		
		//sort map
		Map<Integer, String> sortedMap = new TreeMap<Integer, String>(count2ped); 
		
		//credible interval
		int totalCount = 0;
		int pedCount = 0;
		int inCI = 0;
		int numPairs = truePed.length() / 3;
		int[] pairCorrect = new int[numPairs];

		
		for(int key : sortedMap.keySet()){
			
			//get count
			int count = Math.abs(key);
			totalCount += count;
			pedCount++;
			
			//check if this pedigree is the true one
			String testPed = sortedMap.get(key);
			if(inCI==0 && truePed.equals(testPed)){
				inCI=1;
			}
			
			//check if each pair is correct
			for(int i=0; i<testPed.length()/3; i++) {
				
				if(pairCorrect[i]==0) {
					
					if(testPed.charAt(3*i)==truePed.charAt(3*i) && testPed.charAt(3*i+1)==truePed.charAt(3*i+1) && testPed.charAt(3*i+2)==truePed.charAt(3*i+2))
						pairCorrect[i] = 1;
				}
				
			}
			
			if(totalCount / denom > CI) break;
			
		}
		
		String pairAccuracy = "";
		for(int i : pairCorrect) pairAccuracy += String.format("%d", i);
		
		
		//record results
		String[] toReturn = new String[3];
		toReturn[0] = String.format("%d", inCI);
		toReturn[1] = pairAccuracy;
		toReturn[2] = String.format("%d", pedCount);

		
		
		return toReturn;	
		
		
	}
	
	
	
	//accuracy based on MAP estimate; acc=1 if the kinship coefficients are equal
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
		
		/*
		//print map path
		System.out.println();
		System.out.println(String.format("MCMC MAP likelihood: %.2f", bestLkhd));
				
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				Path bestPath = mapPed[i][j];
				System.out.println(String.format("(%d, %d) : (%d, %d, %d)\n", i, j, bestPath.getUp(), bestPath.getDown(), bestPath.getNumVisit()));
			}
		}
		*/
		
		
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
	
	
	
	//returns squared distance from true kinship coefficient
	public static double[][] kinshipDist(String outPath, String truePath, int totalIndiv, int numIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		
		
		double[][] accuracy = new double[numIndiv][numIndiv];
		
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
				
				double trueKinship = .25*trueOmega[i][j][1] + .5*trueOmega[i][j][2];
				double inferredKinship = .25*pathToOmega.get(key)[1] + .5*pathToOmega.get(key)[2];

				
				accuracy[i][j] = Math.abs(trueKinship-inferredKinship);

			}

			

			
		}
		
		reader.close();
		

		
		
		return accuracy;

		
		
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
			
			if(i>=numIndiv || j>=numIndiv) continue;
			
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
		
		if(i.getNumVisit()==-1 || j.getNumVisit()==-1)
			return false;
		
		else if(omega1[0]==omega2[0] && omega1[1]==omega2[1] && omega1[2]==omega2[2])
			return true;
		
		return false;
		
	}
	
	
}

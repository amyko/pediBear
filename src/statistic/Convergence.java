package statistic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Map;

import dataStructures.Path;
import utility.DataParser;

public class Convergence {

	
	
	public static double[] getTwoHighestLikelihoods(String inPath) throws NumberFormatException, IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		
		//find two highest likelihoods
		double bestLkhd = Double.NEGATIVE_INFINITY;
		double secondLkhd = Double.NEGATIVE_INFINITY;
		
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")) 
				continue;
			
			
			double currLkhd = Double.parseDouble(fields[1]);
			
			
			if(currLkhd > bestLkhd){
				secondLkhd = bestLkhd;
				bestLkhd = currLkhd;
			}
			else if(currLkhd!=bestLkhd && currLkhd > secondLkhd){
				secondLkhd = currLkhd;
			}
			
			
			
		}
		
		reader.close();
		
		
		return new double[]{bestLkhd, secondLkhd};
		
		
	}
	
	
	// get #Ha/#Hb
	public static void getSampleProportion(String inPath, String outPath, double h1, double h2, int sampleRate) throws IOException{

		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		

		//posterior ratio (#h1/#h2)
		double n1 = 0;
		int n2 = 0;
		int n = 0;
		String line;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")) 
				continue;
			
			double currLkhd = Double.parseDouble(fields[1]);

			//count
			if(currLkhd==h1) 
				n1++;
			else if(currLkhd==h2) 
				n2++;
			n++;
			
			//if(n%100000==0){
				if(n2!=0)
					writer.write(String.format("%d\t%.2f\n", n*sampleRate, n1/n2));
			
				//n1 = 0;
				//n2 = 0;
				
			
			//}
			
			
			
			
		}
		
		
		reader.close();
		writer.close();

		
		
		
	}
	
	
	//write likelihood trajectory of the posterior samples
	public static void likelihoodConvergence(String inPath, String outPath, int sampleRate) throws NumberFormatException, IOException{
		
		//open files
		BufferedReader reader = DataParser.openReader(inPath);
		PrintWriter writer = DataParser.openWriter(outPath);
		

		//read likelihood
		String line;
		int n=0;
		while((line=reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			
			if(!fields[0].equals(">")){ 
				n++;
				continue;
			}
			
			double currLkhd = Double.parseDouble(fields[1]);
			
			writer.write(String.format("%d\t%.2f\n", n*sampleRate, currLkhd));
			
			
			
		}
		
		
		reader.close();
		writer.close();
		
		
	}
	
	

		
	//accuracy based on MAP estimate
	public static void distanceFromTruth(String inPath, String outPath, String truePath, int numIndiv, int totalIndiv, Map<Path, double[]> pathToOmega) throws NumberFormatException, IOException{
		
		
		int nPairs = numIndiv*(numIndiv-1)/2;
		
		
		//read true relationship
		double[][][] trueOmega = Accuracy.getTrueOmega(truePath, totalIndiv, pathToOmega);
			
		//read output
		BufferedReader reader = DataParser.openReader(outPath);
		PrintWriter writer = DataParser.openWriter(inPath);
		String line;
		double error = 0;
		
		
		while((line = reader.readLine())!=null){
			
			String[] fields = line.split("\t");
			if(fields[0].equals(">")){
				writer.write(String.format("%.5f", error/nPairs));
				continue;
			}
			

				
			int i = Integer.parseInt(fields[0]);
			int j = Integer.parseInt(fields[1]);
			Path key = new Path(Integer.parseInt(fields[2]) , Integer.parseInt(fields[3]) , Integer.parseInt(fields[4]));
			
			int c = 0;
			for(int k=0; k<3; k++){
				if(Math.abs(trueOmega[i][j][k] - pathToOmega.get(key)[k]) < 1e-13){
					c++;
				}
			}
			
			int add = c==3 ? 1 : 0;
			
			error += add;


			

			
		}
		
		reader.close();
		writer.close();
			
		
		
	}
	
	
	
	
	
}

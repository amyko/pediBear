package statistic;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;


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
			
			if(n2!=0)
				writer.write(String.format("%d\t%.2f\n", n*sampleRate, n1/n2));
			
			
			
		}
		
		
		reader.close();
		writer.close();

		
		
		
	}
	
	
	
	
	
}

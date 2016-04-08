package simulator;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import org.apache.commons.math3.distribution.PoissonDistribution;

import utility.DataParser;

	public class SimulatorStream {
		
		//input columns
		private static int POS = 0;
		private static int A1 = 1;
		private static int A2 = 2;
		
		double recombRate; 
		
		public SimulatorStream(double r) throws IOException{

			recombRate = r;

		}
		
		
		
		//outfile: c1 c2 ... cn
		public void makeChildren(String infoPath, String momPath, String dadPath, String outPath, int momCol, int dadCol, Random rgen, int numChildren, String header) throws IOException{
			
			int recomb_count1=0;
			int recomb_count2=0;
			
			//open files
			BufferedReader momReader = DataParser.openReader(momPath);
			BufferedReader dadReader = DataParser.openReader(dadPath);
			BufferedReader infoReader = DataParser.openReader(infoPath);
			PrintWriter writer = DataParser.openWriter(outPath);
			
			//skip header
			momReader.readLine();
			dadReader.readLine();
			infoReader.readLine();
			
			//write header
			writer.write(header);
			
			//initalize
			int[] curr1hap = new int[numChildren];
			int[] curr2hap = new int[numChildren];
			for (int i=0; i< numChildren; i++){
				curr1hap[i] = rgen.nextBoolean() ? 1 : 0;
				curr2hap[i] = rgen.nextBoolean() ? 1 : 0;
			}
			String[] momFields = momReader.readLine().split("\\s");
			String[] dadFields = dadReader.readLine().split("\\s");
			String[] infoFields = infoReader.readLine().split("\\s");
			int prevPos = Integer.parseInt(infoFields[POS]);
			String momGeno = momFields[momCol];
			String dadGeno = dadFields[dadCol];
			//writer.write(String.format("%s\t%s\t",momGeno,dadGeno));
			for (int i=0; i< numChildren; i++){
				String childGeno = "" + momGeno.charAt(curr1hap[i]) + dadGeno.charAt(curr2hap[i]);
				writer.write(childGeno+"\t");
			}
			writer.write("\n");
	
			
			//for all snps
			String momLine;
			String dadLine;
			String infoLine;
			while((momLine = momReader.readLine())!=null && (dadLine = dadReader.readLine())!=null &&(infoLine = infoReader.readLine())!=null){
				
				//read data
				momFields = momLine.split("\\s");
				dadFields = dadLine.split("\\s");
				infoFields = infoLine.split("\\s");
				momGeno = momFields[momCol];
				dadGeno = dadFields[dadCol];
				
				int currPos = Integer.parseInt(infoFields[POS]);
				double dist = (currPos - prevPos) * recombRate;
				
				if(currPos==prevPos) continue;
				
				//determine recombination
				for (int i=0; i<numChildren; i++){
					PoissonDistribution p = new PoissonDistribution(dist);
					if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
						recomb_count1++;
						curr1hap[i] = 1- curr1hap[i];
					}
					if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
						recomb_count2++;
						curr2hap[i] = 1- curr2hap[i];
					}	
				}
					
				// write result to temp file
				//writer.write(String.format("%s\t%s\t",momGeno,dadGeno));
				for (int i=0; i< numChildren; i++){
					String childGeno = "" + momGeno.charAt(curr1hap[i]) + dadGeno.charAt(curr2hap[i]);
					writer.write(childGeno+"\t");
				}
				writer.write("\n");
				
				//update prev info
				prevPos = currPos;

			}
		
			
			//close files
			writer.close();
			momReader.close();
			dadReader.close();
			
			//System.out.format("%f\t%f\n", recomb_count1/(double)numChildren, recomb_count2/(double)numChildren);
	
			
		}	
		
		
		
		public void makeIndividuals(String outPath, int nInd, int nLocus, int dist, double maf, Random rgen, String header) throws IOException{
			
			PrintWriter writer = DataParser.openWriter(outPath);
			writer.write(header);
			
			for(int locus=0; locus < nLocus; locus++){
				
				//write pos, a1, a2
				writer.write(String.format("%d\t%s\t%s\t", dist*locus , "A", "T"));
				
				for(int i=0; i<nInd; i++){
					
					String geno = "";
					for(int j=0; j<2; j++){
						if(rgen.nextDouble() < maf){
							geno += "A";
						}
						else{
							geno += "T";
						}
					}
					
					writer.write(geno+"\t");
					
				}
				
				writer.write("\n");
				
			}
			
			writer.close();
			
			
		}
		
		
		
		public static String correctMissingData(String geno, String[] fields){
			
			char a1 = fields[A1].charAt(0);
			char a2 = fields[A2].charAt(2);
			
			String toReturn = geno;
			if ((geno.charAt(0)!=a1 && geno.charAt(0)!=a2) || (geno.charAt(1)!=a1 && geno.charAt(1)!=a2)){
				toReturn = ""+a1+a2;
			}
			
			return toReturn;
			
		}
		
}
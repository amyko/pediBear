package simulator;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Random;

import org.apache.commons.math3.distribution.ExponentialDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;

import utility.DataParser;

	
public class SimulatorStreamPed {
		
		//input columns
		private static int CHROM = 0;
		private static int POS = 3;
		
		
		double recombRate; 
		
		public SimulatorStreamPed(double r) throws IOException{

			recombRate = r;

		}
		
		
		
		//outfile in tped format: info_fields, c1 c2 ... cn
		//first individual column = 0
		public void makeChildren(String momPath, String dadPath, String outPath, int momID, int dadID, Random rgen, int numChildren) throws IOException{
	
			//open files
			BufferedReader momReader = DataParser.openReader(momPath);
			BufferedReader dadReader = DataParser.openReader(dadPath);
			PrintWriter writer = DataParser.openWriter(outPath);

			
			//adjusted mom and dad columns
			int dadCol = 2*(dadID+2);
			int momCol = 2*(momID+2);
			
			//String[] mPath = momPath.split("/");
			//String[] dPath = dadPath.split("/");
			
			//System.out.println(String.format("%s %s", mPath[mPath.length-1], dPath[dPath.length-1]));
			//System.out.println(String.format("%d %d", momID, dadID));
			
			
			//initialize
			int[] curr1hap = new int[numChildren];
			int[] curr2hap = new int[numChildren];

			//for all snps
			String momLine;
			String dadLine;
			int prevChrom = -1;
			int prevPos = 0;
			//int n = 0;


			while((momLine = momReader.readLine())!=null && (dadLine = dadReader.readLine())!=null){
				
				//read data
				String[] momFields = momLine.split("\\s");
				String[] dadFields = dadLine.split("\\s");
				int currChrom = Integer.parseInt(momFields[CHROM]);
				int currPos = Integer.parseInt(momFields[POS]);

				
				if(currPos==prevPos) continue;

				//if first snp, randomly choose start haplotype	
				if(currChrom != prevChrom){
					
					//randomly choose start haplotype
					for (int i=0; i< numChildren; i++){
						curr1hap[i] = rgen.nextBoolean() ? 1 : 0;
						curr2hap[i] = rgen.nextBoolean() ? 1 : 0;
					}
					
					

				}
				
				
				

				//if not, determine recombination
				else{
	
					
					double dist = (currPos - prevPos) * recombRate;
					
					for (int i=0; i<numChildren; i++){
						PoissonDistribution p = new PoissonDistribution(dist);
						if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
							curr1hap[i] = 1- curr1hap[i];
							//n++;
						}
						if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
							curr2hap[i] = 1- curr2hap[i];
						}	
					}
					
						
					
				}
				
				
					
				// write to outfile
				writer.write(String.format("%s %s %s %s ", momFields[0], momFields[1], momFields[2], momFields[3])); //header
				for (int i=0; i< numChildren; i++){				
					writer.write(String.format("%s %s ", momFields[momCol + curr1hap[i]], dadFields[dadCol + curr2hap[i]]));
				}
				writer.write("\n");
				
				
				//update pos & chrom
				prevChrom = currChrom;
				prevPos = currPos;

			}
		
			
			//close files
			writer.close();
			momReader.close();
			dadReader.close();
			//System.out.println(n);
			
			
			
		}
	

		//outfile in tped format: info_fields, c1 c2 ... cn
		//first individual column = 0
		public void makeChildrenExp(String momPath, String dadPath, String outPath, int momID, int dadID, Random rgen) throws IOException{
	
			//open files
			BufferedReader momReader = DataParser.openReader(momPath);
			BufferedReader dadReader = DataParser.openReader(dadPath);
			PrintWriter writer = DataParser.openWriter(outPath);

			
			//adjusted mom and dad columns
			int dadCol = 2*(dadID+2);
			int momCol = 2*(momID+2);

			
			//initialize
			int curr1hap = 0;
			int curr2hap = 0;
			ExponentialDistribution exp = new ExponentialDistribution(recombRate);

			//for all snps
			String momLine;
			String dadLine;
			int prevChrom = -1;
			int prevPos = 0;
			int nextBreakPoint1 = 0;
			int nextBreakPoint2 = 0;
			//int n = 0;


			while((momLine = momReader.readLine())!=null && (dadLine = dadReader.readLine())!=null){
				
				//read data
				String[] momFields = momLine.split("\\s");
				String[] dadFields = dadLine.split("\\s");
				int currChrom = Integer.parseInt(momFields[CHROM]);
				int currPos = Integer.parseInt(momFields[POS]);

				
				if(currPos==prevPos) continue;

				//if first snp, randomly choose start haplotype	
				if(currChrom != prevChrom){
					
					System.out.println(currChrom);
					
					//randomly choose start haplotype
					curr1hap = rgen.nextBoolean() ? 1 : 0;
					curr2hap = rgen.nextBoolean() ? 1 : 0;
				
					//choose next break points
					nextBreakPoint1 = currPos + (int) (1/exp.sample());
					nextBreakPoint2 = currPos + (int) (1/exp.sample());
				
				}
				
				
				
				//if not, determine recombination
				else{
					//recombine
					if(currPos >= nextBreakPoint1){
						curr1hap = 1 - curr1hap;
						nextBreakPoint1 = currPos + (int) (1/exp.sample());
					}
					if(currPos >= nextBreakPoint2){
						curr2hap = 1 - curr2hap;
						nextBreakPoint2 = currPos + (int) (1/exp.sample());
					}						
					
				}
				
				
					
				// write to outfile
				writer.write(String.format("%s %s %s %s ", momFields[0], momFields[1], momFields[2], momFields[3])); //header			
				writer.write(String.format("%s %s ", momFields[momCol + curr1hap], dadFields[dadCol + curr2hap]));	
				writer.write("\n");

				
				//update pos & chrom
				prevChrom = currChrom;
				prevPos = currPos;
				
				
				

			}
		
			
			//close files
			writer.close();
			momReader.close();
			dadReader.close();
			//System.out.println(n);
			
			
			
		}
	

		
		
		public void addError(String inPath, String outPath, double epsilon, Random rGen) throws IOException{
			
			//open files
			BufferedReader reader = DataParser.openReader(inPath);
			PrintWriter writer = DataParser.openWriter(outPath);
			

			String line;
			while((line = reader.readLine())!=null){
				
				//read line
				String[] fields = line.split("\\s");
				
				//info fields
				for(int i=0; i<4; i++){
					writer.write(fields[i]+" ");
				}
				
				//for each genotype
				for(int i=4; i<fields.length; i++){

						
					if(rGen.nextDouble() < epsilon){
						if(fields[i].charAt(0)=='A')
							writer.write("T");
						else
							writer.write("A");
					}
					else{
						writer.write(fields[i].charAt(0));
					}
						
					
					
					writer.write(" ");
					
					
				}
				
				writer.write("\n");
					
				
			}
			
			
			//close files
			reader.close();
			writer.close();
			
		}
		
		
		//returns ped file line of the offspring
		public String makeChildrenReturnPedline(String momPath, String dadPath, int momID, int dadID, Random rgen, int numChildren) throws IOException{
			
			//open files
			BufferedReader momReader = DataParser.openReader(momPath);
			BufferedReader dadReader = DataParser.openReader(dadPath);

			//adjusted mom and dad columns
			int dadCol = 2*(dadID+2);
			int momCol = 2*(momID+2);
			
			
			//initialize
			int[] curr1hap = new int[numChildren];
			int[] curr2hap = new int[numChildren];

			//for all snps
			String momLine;
			String dadLine;
			int prevChrom = -1;
			int prevPos = 0;
			String toReturn = "";


			while((momLine = momReader.readLine())!=null && (dadLine = dadReader.readLine())!=null){
				
				//read data
				String[] momFields = momLine.split("\\s");
				String[] dadFields = dadLine.split("\\s");
				int currChrom = Integer.parseInt(momFields[CHROM]);
				int currPos = Integer.parseInt(momFields[POS]);

			
				
				//if first snp, randomly choose start haplotype	
				if(currChrom != prevChrom){
					
					//randomly choose start haplotype
					for (int i=0; i< numChildren; i++){
						curr1hap[i] = rgen.nextBoolean() ? 1 : 0;
						curr2hap[i] = rgen.nextBoolean() ? 1 : 0;
					}
					
					

				}
				
				
				

				//if not, determine recombination
				else{
	
					
					double dist = (currPos - prevPos) * recombRate;
					
					for (int i=0; i<numChildren; i++){
						PoissonDistribution p = new PoissonDistribution(dist);
						if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
							curr1hap[i] = 1- curr1hap[i];
							//n++;
						}
						if(p.sample() % 2 == 1){	//odd number of recombinations happened, switch
							curr2hap[i] = 1- curr2hap[i];
						}	
					}
					
						
					
				}
				
				
					

				for (int i=0; i< numChildren; i++){	
					toReturn += String.format("%s %s ", momFields[momCol + curr1hap[i]], dadFields[dadCol + curr2hap[i]]);
				}
				
				//update pos & chrom
				prevChrom = currChrom;
				prevPos = currPos;

			}
		
			
			//close files
			momReader.close();
			dadReader.close();

			return toReturn;
			
			
		}

}




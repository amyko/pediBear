package simulator;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.distribution.PoissonDistribution;

import dataStructures.Pedigree;
import utility.DataParser;

	public class SimulatorStream {
		
		//input columns
		private static int POS = 0;
		
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
		

		
		public void simulatePopulation(String inPath, String outPath, String pedPath, int oldestGen, Random rGen, int chrStart, int chrEnd, int[] cols, int howManyUnrelated) throws IOException{
			
			
			//variables
			int marriageNum = 0;
			int currIdx = 2;
			int unrelatedIdx = 8;
			int currCol = 2; //current column being written
			List<Integer> breeders = new ArrayList<Integer>();
			List<Integer> breedersRealCol = new ArrayList<Integer>();
			
			
			//open files
			PrintWriter pedfile = DataParser.openWriter(pedPath);
			
			//outfiles
			String[][] outPaths = new String[2][chrEnd-chrStart];
			for(int chr=chrStart; chr < chrEnd; chr++){
				for(int i=0; i<2; i++){
					outPaths[i][chr-1] = String.format("%stemp.%d.%d", outPath, chr, i);
				}
			}
			
			int curr = 0;
			int prev = 1;
			
			
			for(int chr=chrStart; chr < chrEnd; chr++){
				
				String genoPath = inPath + chr;
				
				//first generation
				makeChildren(genoPath, genoPath, genoPath, outPaths[curr][chr-1], unrelatedIdx , unrelatedIdx+1, rGen, 3, String.format("%d\t%d\t%d\n", currIdx, currIdx+1, currIdx+2));
				//concatenate parents with children
				DataParser.concatFiles(new String[]{genoPath, outPaths[curr][chr-1]}, outPaths[prev][chr-1], new int[][]{new int[]{unrelatedIdx, unrelatedIdx+1}, new int[]{0,1,2}});
			}

			//write relationship
			writePedFile(pedfile, currIdx-2, currIdx-1, new int[]{currIdx, currIdx+1, currIdx+2});
			
			//add moms to breeder list
			for(int i : new int[]{currCol, currCol+1, currCol+2})
				breeders.add(i);
			for(int i : new int[]{currIdx, currIdx+1, currIdx+2})
				breedersRealCol.add(i);
			
			
			//increment
			unrelatedIdx += 2; //for mom and dad
			currIdx += 3; //3 children
			currCol += 3; //3 breeders
			prev = curr;
			curr = (curr+1) % 2;
			marriageNum++;
			
			
			//later generations
			for(int gen=oldestGen-1; gen>0; gen--){
				
				List<Integer> nextBreeders = new ArrayList<Integer>();
				List<Integer> nextBreedersRealCol = new ArrayList<Integer>();
				
				for(int idx=0; idx<breeders.size(); idx++){
					
					//file to store children
					String childrenPath = outPath + "children";
					
					//increment current index for making a dad
					int dadIdx = currIdx;
					currIdx++;
					int momIdx = breedersRealCol.get(idx);
					int momCol = breeders.get(idx);
					
					//first marriage
					for(int chr=chrStart; chr < chrEnd; chr++){
						
						String genoPath = inPath + chr;
						
						makeChildren(genoPath, outPaths[curr][chr-1], genoPath, childrenPath, momCol, unrelatedIdx, rGen, 3, String.format("%d\t%d\t%d\n", currIdx, currIdx+1, currIdx+2));
						DataParser.concatFiles(new String[]{outPaths[curr][chr-1], childrenPath}, outPaths[prev][chr-1], new int[][]{new int[]{}, new int[]{0,1,2}});
					}

					//add children to breeder list
					for(int i : new int[]{currCol, currCol+1, currCol+2})
						nextBreeders.add(i);
					for(int i : new int[]{currIdx, currIdx+1, currIdx+2})
						nextBreedersRealCol.add(i);
					
					//write relationship
					writePedFile(pedfile, momIdx, dadIdx, new int[]{currIdx, currIdx+1, currIdx+2});
					
					//increment
					unrelatedIdx++; //for dad
					currIdx += 3; //for 3 children
					currCol += 3; //for 3 breeders
					prev = curr;
					curr = (curr+1) % 2;
					marriageNum++;

					
					//another marriage
					if(marriageNum % 8 == 0){
						
						//increment current index for making a dad
						dadIdx = currIdx;
						currIdx++;
						
						//first marriage
						for(int chr=chrStart; chr < chrEnd; chr++){
							
							String genoPath = inPath + chr;
							makeChildren(genoPath, outPaths[curr][chr-1], genoPath, childrenPath, momCol, unrelatedIdx, rGen, 2, String.format("%d\t%d\n", currIdx, currIdx+1));
							DataParser.concatFiles(new String[]{outPaths[curr][chr-1], childrenPath}, outPaths[prev][chr-1], new int[][]{new int[]{}, new int[]{0,1}});
						}
						
						
						//add children to breeder list
						for(int i : new int[]{currCol, currCol+1})
							nextBreeders.add(i);
						for(int i : new int[]{currIdx, currIdx+1})
							nextBreedersRealCol.add(i);
						
						//write relationship
						writePedFile(pedfile, momIdx, dadIdx, new int[]{currIdx, currIdx+1});
						
						//increment
						unrelatedIdx++; //for dad
						currIdx += 2; //for 2 children
						currCol += 2; //for 2 breeders
						prev = curr;
						curr = (curr+1) % 2;
						marriageNum++;

					}

				}				
				
				//increment
				breeders = new ArrayList<Integer>(nextBreeders);
				breedersRealCol = new ArrayList<Integer>(nextBreedersRealCol);
				
				
			}
			
			
			//close files
			pedfile.close();
			
			System.out.println(unrelatedIdx);
			
			//sample
			sampleIndiv(outPath, curr, rGen, chrStart, chrEnd, cols, howManyUnrelated, unrelatedIdx, inPath);
		
			
			
		}
		
		
		public void writePedFile(PrintWriter pedfile, int mom, int dad, int[] children){
			
			for(int i : children){
				pedfile.write(String.format("%d\t%d\t%d\n", i, mom, dad));
			}
			
			pedfile.flush();
			
			
		}
		
		
		
		public void sampleIndiv(String filePath, int curr, Random rGen, int chrStart, int chrEnd, int[] cols, int howManyUnrelated, int unrelatedIdx, String unrelatedGenoPath) throws IOException{
			
			
			
			for(int chr=chrStart; chr<chrEnd; chr++){
			
				//open files
				BufferedReader infile = DataParser.openReader(String.format("%stemp.%d.%d", filePath, chr, curr));
				BufferedReader unrelFile = DataParser.openReader(unrelatedGenoPath + chr);
				PrintWriter outfile = DataParser.openWriter(filePath + chr);
			
				//skip header
				String[] headerFields = infile.readLine().split("\t");
				String[] headerFields2 = unrelFile.readLine().split("\t");
				
				String header = "";
				
				for(int i=0; i<cols.length; i++){
					header += headerFields[cols[i]] + "\t";
				}
				for(int i=0; i<howManyUnrelated; i++){
					header += headerFields2[unrelatedIdx + i] + "\t";
				}
				header += "\n";
				
				
				//write header
				outfile.write(header);
				System.out.println(header);
				
			
				//copy sampled individuals to outfile
				String line1;
				String line2;
				while((line1 = infile.readLine()) != null && (line2 = unrelFile.readLine()) != null){
					
					//simulated individuals
					String[] fields = line1.split("\t");
					for(int i : cols){
						outfile.write(String.format("%s\t", fields[i]));
					}
					
					//append unrelated
					fields = line2.split("\\s");
					for(int i=0; i<howManyUnrelated; i++){
						outfile.write(String.format("%s\t", fields[unrelatedIdx+i]));
					}
					
					
					
					outfile.write("\n");
					
					
				}
				
				
				//close file
				infile.close();
				outfile.close();
				
				
			}

			
		}
		
		

		
		public void writeTruePathFile(String inPath, String outPath, int[] cols) throws IOException{
			
			new Pedigree(inPath, outPath, 10, cols);

			
		}
		
		
		
		
}
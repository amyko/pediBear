package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.Random;
import java.io.PrintWriter;

import likelihood.LDStream;
import simulator.SimulatorStream;
import utility.DataParser;

//makes a test pedigree and then computes its likelihood under the pairwise core and using lander-green
public class testIBD {
		
		public static double r = 1.3e-8; 
		public static double seqError = 0.01;
		public static String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/";
		public static int back = 30000;
		public static Random rgen = new Random(5265704L);


	
		
		public static void simulateTwoHalfCousins(String inPath, String commonAncPath, String unrelatedAncPath, int numGen, int start, int t, String testName) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);

			int commonAnc = start;
			start++;

			//file names
			String genoOutPath = dir + "ibd/"+ testName + ".geno."+t;
			String ancOutPath = dir + "ibd/" + testName + ".anc";
			String tempGeno1 = dir + "ibd/temp1.geno";
			String tempGeno2 = dir + "ibd/temp2.geno";
			String tempAnc1 = dir + "ibd/temp1.anc";
			String tempAnc2 = dir + "ibd/temp2.anc";
			
					
			//first generation
			sim.makeChildrenIBD(inPath, inPath, inPath, commonAncPath, unrelatedAncPath, tempGeno1, tempAnc1, commonAnc, start, 0,0, rgen, 1, "c1\n");
			start++;
			sim.makeChildrenIBD(inPath, inPath, inPath, commonAncPath, unrelatedAncPath, tempGeno2, tempAnc2, commonAnc, start, 0,0, rgen, 1, "c1\n");
			start++;
			
			//concatenate into one output file
			DataParser.concatFiles(new String[]{tempGeno1, tempGeno2}, genoOutPath, new int[][]{new int[]{0}, new int[]{0}});
			DataParser.concatFiles(new String[]{tempAnc1, tempAnc2}, ancOutPath, new int[][]{new int[]{0}, new int[]{0}});


			
			//later generations
			int currGen = numGen - 1;
			while(currGen > 0){
				
	
				sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno1, tempAnc1, 0, start, 0, 0, rgen, 1, "c1\n");
				start++;
				sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno2, tempAnc2, 1, start, 1, 0, rgen, 1, "c1\n");
				start++;
				
				//concatenate into one output file
				DataParser.concatFiles(new String[]{tempGeno1, tempGeno2}, genoOutPath, new int[][]{new int[]{0}, new int[]{0}});
				DataParser.concatFiles(new String[]{tempAnc1, tempAnc2}, ancOutPath, new int[][]{new int[]{0}, new int[]{0}});
				
				//increment
				currGen--;
				
				
			}
			
			
			
			//record ibd status
			PrintWriter ibdWriter = DataParser.openWriter(ancOutPath+"."+t);
			BufferedReader ibdReader = DataParser.openReader(ancOutPath);
			
			ibdWriter.write("ibd\n");
			ibdReader.readLine();
			String line;
			while((line=ibdReader.readLine())!=null){
				String[] fields = line.split("\t");
				if((fields[0].charAt(0)=='1' || fields[0].charAt(1)=='1') && (fields[1].charAt(0)=='1' || fields[1].charAt(1)=='1')){
					ibdWriter.write(String.format("%d\n", 1));
				}
				else if((fields[0].charAt(0)=='2' || fields[0].charAt(1)=='2') && (fields[1].charAt(0)=='2' || fields[1].charAt(1)=='2')){
					ibdWriter.write(String.format("%d\n", 1));
				}
				else{
					ibdWriter.write(String.format("%d\n", 0));
				}
			}
			
			ibdWriter.close();
			ibdReader.close();

			
			
		}
		
		
		
		public static void simulateTwoHalfCousinsPedigree(String inPath, String commonAncPath, String unrelatedAncPath, int start, int t, String testName) throws IOException{
			
			//initialize simulator
			SimulatorStream sim = new SimulatorStream(r);

			int commonAnc = start;
			start++;

			//file names
			String genoOutPath = dir + "ibd/"+ testName + ".geno."+t;
			String ancOutPath = dir + "ibd/" + testName + ".anc";
			String tempGeno1 = dir + "ibd/temp1.geno";
			String tempGeno2 = dir + "ibd/temp2.geno";
			String tempAnc1 = dir + "ibd/temp1.anc";
			String tempAnc2 = dir + "ibd/temp2.anc";
			String tempGeno3 = dir + "ibd/temp3.geno";
			String tempGeno4 = dir + "ibd/temp4.geno";
			String temp = dir + "ibd/temp";
			
					
			//first generation
			sim.makeChildrenIBD(inPath, inPath, inPath, commonAncPath, unrelatedAncPath, tempGeno1, tempAnc1, commonAnc, start, 0,0, rgen, 1, "c1\n");
			start++;
			sim.makeChildrenIBD(inPath, inPath, inPath, commonAncPath, unrelatedAncPath, tempGeno2, tempAnc2, commonAnc, start, 0,0, rgen, 1, "c1\n");
			start++;
			
			//concatenate into one output file
			DataParser.concatFiles(new String[]{tempGeno1, tempGeno2}, genoOutPath, new int[][]{new int[]{0}, new int[]{0}});
			DataParser.concatFiles(new String[]{tempAnc1, tempAnc2}, ancOutPath, new int[][]{new int[]{0}, new int[]{0}});


			//second generation
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno1, tempAnc1, 0, start, 0,0, rgen, 2, "c1\tc2\n");
			start++;
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno2, tempAnc2, 1, start, 1,0, rgen, 2, "c1\tc2\n");
			start++;
			
			//concatenate into one output file
			DataParser.concatFiles(new String[]{tempGeno1, tempGeno2}, genoOutPath, new int[][]{new int[]{0,1}, new int[]{0,1}});
			DataParser.concatFiles(new String[]{tempAnc1, tempAnc2}, ancOutPath, new int[][]{new int[]{0,1}, new int[]{0,1}});
			
			
			//TODO start here
			//third generations
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno1, temp, 0, start, 0, 0, rgen, 1, "c1\n");
			start++;
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno2, tempAnc1, 1, start, 1, 0, rgen, 1, "c1\n"); 
			start++;
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno3, tempAnc2, 2, start, 2, 0, rgen, 1, "c1\n");
			start++;
			sim.makeChildrenIBD(inPath, genoOutPath, inPath, ancOutPath, unrelatedAncPath, tempGeno4, temp, 3, start, 3, 0, rgen, 1, "c1\n");
			start++;
			
			//concatenate into one output file
			DataParser.concatFiles(new String[]{tempGeno2, tempGeno3, tempGeno1, tempGeno4}, genoOutPath, new int[][]{new int[]{0}, new int[]{0}, new int[]{0}, new int[]{0}});
			DataParser.concatFiles(new String[]{tempAnc1, tempAnc2}, ancOutPath, new int[][]{new int[]{0}, new int[]{0}});

			
			
			//record ibd status
			PrintWriter ibdWriter = DataParser.openWriter(ancOutPath+"."+t);
			BufferedReader ibdReader = DataParser.openReader(ancOutPath);
			
			ibdWriter.write("ibd\n");
			ibdReader.readLine();
			String line;
			while((line=ibdReader.readLine())!=null){
				String[] fields = line.split("\t");
				if((fields[0].charAt(0)=='1' || fields[0].charAt(1)=='1') && (fields[1].charAt(0)=='1' || fields[1].charAt(1)=='1')){
					ibdWriter.write(String.format("%d\n", 1));
				}
				else{
					ibdWriter.write(String.format("%d\n", 0));
				}
			}
			
			ibdWriter.close();
			ibdReader.close();

			
			
		}
		

		
		public static void makePosFile(String inPath, String outPath) throws IOException{
			
			
			BufferedReader reader = DataParser.openReader(inPath);
			PrintWriter writer = DataParser.openWriter(outPath);
			
			//skip header
			reader.readLine();
			
			String line;
			while((line=reader.readLine())!=null){
				
				String[] fields = line.split("\t");
				
				double pos = Double.parseDouble(fields[0]) / 1e6 * 1.3;
				
				writer.write(String.format("%f\n", pos));
				
			}
			
			reader.close();
			writer.close();
			
			
		}
		

		//makes anc file for an individual with a given anc
		public static void makeAncFile(String outPath, int numSNP, int anc1, int anc2) throws IOException{
			
			
			PrintWriter writer = DataParser.openWriter(outPath);
			
			// write header
			writer.write("anc\n");
			
			//write anc
			for(int i=0; i<numSNP; i++)
				writer.write(String.format("%d%d\n", anc1, anc2));
			
			writer.close();
			
			
			
		}

		
		public static void main(String[] args) throws IOException{

			//dir
			String dataDir = dir + "unrelated/";
			String simDir = dir + "ibd/";
			String dataPath = dataDir+"msprime.geno.pruned.1";


			//simulation parameters
			int numGen = 3;
			int start = 8;
			String testName = "test";
			
			/*
			int numSnps = DataParser.countLines(dataPath) - 1;
			makeAncFile(simDir+"common.anc", numSnps, 1, 2);
			makeAncFile(simDir+"unrelated.anc", numSnps, 0, 0);
			*/
	
			/*
			//init simulator
			SimulatorStream sim = new SimulatorStream(r);
			
			
			//prune ld
			System.out.println("Pruning LD");
			LDStream.thin(dataDir+"msprime.geno.1", dataPath, .2, rgen);



			 
			//make anc files: 1 for common ancestor, 0 for everyone else
			int numSnps = DataParser.countLines(dataPath) - 1;
			sim.makeAncFile(simDir+"common.anc", numSnps, 1);
			sim.makeAncFile(simDir+"unrelated.anc", numSnps, 0);
			
			
			
			//compute info
			System.out.println("Computing info");
			String inPath = dataDir+"msprime.geno.pruned.1";
			String outPath = dataDir+"msprime.info.pruned.1";
			LDStream.writeLdOutfile(inPath, outPath, back);
			*/
			 
			
			//make position file
			makePosFile(dataPath, simDir + testName+".pos");			
			
			for(int t=0; t<300; t++){
				
				System.out.println(t);
		
				simulateTwoHalfCousins(dataDir+"msprime.geno.pruned.1", simDir+"common.anc", simDir+"unrelated.anc", numGen, start, t, testName);

			}
			


		}
		
		


}
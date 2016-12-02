package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import statistic.Accuracy;
import utility.DataParser;
import dataStructures.Path;


//infers the relationship from pairwise likelihood
public class PairwiseLRT {

	
	public static Map<Path, double[][]> path2table(String filePath, int numIndiv) throws NumberFormatException, IOException{
					
		Map<Path, double[][]> toReturn = new HashMap<Path, double[][]>();
		
		//open files
		BufferedReader reader = DataParser.openReader(filePath);
		
		String line;
		Path path = null;
		double toAdd = Math.log(numIndiv)/numIndiv;
		toAdd = 1;
		
		
		while((line=reader.readLine())!=null){

			
			String[] fields = line.split("\t");
			
			//update relationship
			if(fields[0].equals(">")){
				
				int up = Integer.parseInt(fields[1]);
				int down = Integer.parseInt(fields[2]);
				int nVisit = Integer.parseInt(fields[3]);
				
				double[][] lkhdTable = new double[numIndiv][numIndiv];
				path = new Path(up, down, nVisit);

				toReturn.put(path, lkhdTable);

				continue;
				
			}

			//read likelihood
			int i = Integer.parseInt(fields[0]);
			int j = Integer.parseInt(fields[1]);
			double lkhd = Double.parseDouble(fields[2]);

			if(i<numIndiv && j<numIndiv){
				
				//penalty
				double coeff = 1;
				if(path.getNumVisit()==0){
					 coeff = 1;
				}
				
				
				toReturn.get(path)[i][j] = lkhd + coeff*toAdd;
				toReturn.get(path)[j][i] = lkhd + coeff*toAdd;
			}
	
			

		}
		
		//close file
		reader.close();
			
		
		
		return toReturn;
		
		
	}
	
	
	//returns relationships that are not rejected 
	public static List<Path> pairwiseLRT(Map<Path, double[][]> path2lkhd, int i, int j, double c){
		
		
		List<Path> toReturn = new ArrayList<Path>();
		

		
		for(Path nullRel : path2lkhd.keySet()){
			
			//TODO
			if(nullRel.getDown() + nullRel.getUp() > 8) continue;
			
			double nullLkhd = path2lkhd.get(nullRel)[i][j];
			
			//get the supremum of L(H1)
			double maxL1 = Double.NEGATIVE_INFINITY;
			
			for(Path altRel : path2lkhd.keySet()){
				
				//TODO
				if(altRel.getDown() + altRel.getUp() > 8) continue;
				
				
				double altLkhd = path2lkhd.get(altRel)[i][j];
			
				//skip null lkhd
				if(altLkhd==nullLkhd) continue;
				
				if(altLkhd > maxL1) maxL1 = altLkhd;
			}
			
			
			//cannot reject null
			if(nullLkhd - maxL1 >= c){
				toReturn.add(nullRel);
			}
			
			
		}
		
		return toReturn;
		
		
		
	}
	
	
	//write accuracy and distance
	public static void LRT(String lkhdPath, double[][][] trueOmega, Map<Path, double[]> pathToKinship, PrintWriter pathWriter, PrintWriter mapWriter, PrintWriter distWriter, int numIndiv, double c) throws NumberFormatException, IOException{
		
		//get path2lkhd table
		Map<Path, double[][]> path2lkhd = path2table(lkhdPath, numIndiv);


		
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				
				pathWriter.write(String.format("%d\t%d\t", i , j));
				mapWriter.write(String.format("%d\t%d\t", i , j));
				distWriter.write(String.format("%d\t%d\t", i , j));
				
				List<Path> inferredRel = pairwiseLRT(path2lkhd, i, j, c);

				if(inferredRel.size()==0){
					
					pathWriter.write("x\n");
					mapWriter.write("0\n");
					distWriter.write("0\n");
					
				}
				
				else{
					
					double acc = 0;
					double dist = 0;
					
					for(Path rel : inferredRel){
					
						pathWriter.write(String.format("%d\t%d\t%d\t", rel.getUp(), rel.getDown(), rel.getNumVisit()));

						//accuracy
						int count = 0;
						for(int k=0; k<3; k++){
							if(Math.abs(trueOmega[i][j][k] - pathToKinship.get(rel)[k]) < 1e-13){
								count++;
							}
						}
						
						int add = count==3 ? 1 : 0;
						acc += add;
						
						//dist
						double trueKinship = .25*trueOmega[i][j][1] + .5*trueOmega[i][j][2];
						double inferredKinship = .25*pathToKinship.get(rel)[1] + .5*pathToKinship.get(rel)[2];
						dist += Math.abs(trueKinship - inferredKinship);
						
						
						
					}
					
					acc = acc / inferredRel.size();
					dist = dist / inferredRel.size();
					mapWriter.write(String.format("%f\n", acc));
					distWriter.write(String.format("%f\n", dist));
					

					
					
				}
				

				pathWriter.write("\n");
				
				
			}
		}
		
		
		
		
		
	}
	
	
	
	public static void main(String[] args) throws IOException{
		
		//param
		int numIndiv = 18;
		double c = 0;
		
		//files
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String pathToOmegaPath = dir + "pathToOmega.txt";
		String truePath = dir + "results/sim4.true";

		
		//get true omega
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmegaPath);
		double[][][] trueOmega = Accuracy.getTrueOmega(truePath, numIndiv, pathToKinship);
		
		//test name
		String testName = "sim4";
		
		
		//open outfiles
		String outPath = dir + "results/" +testName + ".12.pairwise";
		PrintWriter mapWriter = DataParser.openWriter(outPath+".mapAcc");
		PrintWriter outWriter = DataParser.openWriter(outPath+".out");
		PrintWriter distWriter = DataParser.openWriter(outPath+".kinshipDist");
		
		//run pairwise test
		for(int t=0; t<100; t++){
			
			String lkhdPath = String.format(dir + "simPed4/sim4.%d.pairwise", t);
			
			//write header
			outWriter.write(String.format(">\t%d\n", t));
			mapWriter.write(String.format(">\t%d\n", t));
			distWriter.write(String.format(">\t%d\n", t));
			
			LRT(lkhdPath, trueOmega, pathToKinship, outWriter, mapWriter, distWriter, numIndiv, c);
	
		}
		
		mapWriter.close();
		outWriter.close();
		distWriter.close();

		/*
		int[] myLengths = new int[]{10,20,30,40};
		String[] r_sqrs = new String[]{"0.100", "0.075", "0.050", "0.025"};
		
		for(int cl : myLengths){
			
			for(String rsqr : r_sqrs){
				
				System.out.println(String.format("%d %s", cl, rsqr));
				
				
				
				for(int gen=2; gen<=4; gen++){
					
					//test name
					String testName = String.format("cousins%d.%dmorgan.%s", gen,cl, rsqr);
					
					//open outfiles
					String outPath = dir + "results/" +testName + ".pairwise";
					PrintWriter mapWriter = DataParser.openWriter(outPath+".mapAcc");
					PrintWriter outWriter = DataParser.openWriter(outPath+".out");
					PrintWriter distWriter = DataParser.openWriter(outPath+".kinshipDist");
				
					//true
					truePath = String.format(dir + "results/cousins%d.true", gen);
					trueOmega = Accuracy.getTrueOmega(truePath, numIndiv, pathToKinship);
					
					//lkhd test					
					String lkhdPath = String.format(dir + "cousins/cousins%d.%dmorgan.%s.pairwise", gen, cl, rsqr);
					
					//write header
					outWriter.write(String.format(">\t%d\n", 0));
					mapWriter.write(String.format(">\t%d\n", 0));
					distWriter.write(String.format(">\t%d\n", 0));
					
					LRT(lkhdPath, trueOmega, pathToKinship, outWriter, mapWriter, distWriter, numIndiv, c);
				
					
					
					mapWriter.close();
					outWriter.close();
					distWriter.close();
					
				}
				
				
	
				
			}
			
		}
		*/
		
		
		
		
			
			
	
		
		

		

		
		
		
		
	}
	
	
}
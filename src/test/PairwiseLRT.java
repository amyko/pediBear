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
				toReturn.get(path)[i][j] = lkhd;
				toReturn.get(path)[j][i] = lkhd;
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
			
			double nullLkhd = path2lkhd.get(nullRel)[i][j];
			
			//get the supremum of L(H1)
			double maxL1 = Double.NEGATIVE_INFINITY;
			
			for(Path altRel : path2lkhd.keySet()){
				
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
		String testName = "sim3";
		String outPath = dir + "results/" +testName + ".pairwise";
		String pathToOmegaPath = dir + "pathToOmega.txt";
		String truePath = dir + "results/sim3.true";

		
		//get true omega
		Map<Path, double[]> pathToKinship = Accuracy.getPathToOmega(pathToOmegaPath);
		double[][][] trueOmega = Accuracy.getTrueOmega(truePath, numIndiv, pathToKinship);

		//open outfiles
		PrintWriter mapWriter = DataParser.openWriter(outPath+".mapAcc");
		PrintWriter outWriter = DataParser.openWriter(outPath+".out");
		PrintWriter distWriter = DataParser.openWriter(outPath+".kinshipDist");
		
		//run pairwise test
		for(int t=0; t<1; t++){
			
			String lkhdPath = dir + "simPed3/sim3.10morgan.0.050."+t+".pairwise";
			
			//write header
			outWriter.write(String.format(">\t%d\n", t));
			mapWriter.write(String.format(">\t%d\n", t));
			distWriter.write(String.format(">\t%d\n", t));
			
			LRT(lkhdPath, trueOmega, pathToKinship, outWriter, mapWriter, distWriter, numIndiv, c);
	
		}
		
		mapWriter.close();
		outWriter.close();
		distWriter.close();
		

		
		
		
		
	}
	
	
}
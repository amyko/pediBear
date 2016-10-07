package test;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import dataStructures.Path;
import dataStructures.Pedigree;
import statistic.Accuracy;
import utility.DataParser;

public class primusAcc {
	
	
	public static void main(String[] args) throws IOException{
		
		//directories
		String testName = "test12.pruned.10k";
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String outPath = dir + "results/test12.pruned.10k.plink";
		String truePath = dir + "results/test12.true";
		
		//open files
		PrintWriter accWriter = DataParser.openWriter(outPath+".mapAcc");
		PrintWriter distWriter = DataParser.openWriter(outPath+".kinshipDist");
		Path[][] trueRel = Accuracy.getTruePath(truePath, 20);
		Map<Path, double[]> pathToOmega = Accuracy.getPathToOmega(dir + "pathToOmega.txt");
		
		for(int t=0; t<100; t++){
			
			System.out.println(t);
			
			String resultPath = dir + String.format("test12.10k.primus.plink.results/result.%d/", t);
			
			//accuracy matrix (i,j,acc)
			double[][] accMat = new double[20][20];
			double[][] distMat = new double[20][20];
			
			//////////////////////////////////////////////////////////////
			//write accuracy for unrelated individuals
			BufferedReader reader = DataParser.openReader(resultPath + String.format("%s.%d.genome_unrelated_samples.txt", testName, t));
			reader.readLine();
			
			String line;
			List<Integer> unrel = new ArrayList<Integer>();
			while((line = reader.readLine())!=null){
				String[] fields = line.split("\t");
				unrel.add(Integer.parseInt(fields[1]) - 1);
			}
			reader.close();
			
			for(int ind1 : unrel){
				for(int ind2 : unrel){
		
					if(ind1>=ind2) continue;
					
					//accuracy
					int acc = 0;
					if(trueRel[ind1][ind2].getNumVisit()==0) acc = 1;
					accMat[ind1][ind2] = acc;
					
					//distance
					Path trueKey = trueRel[ind1][ind2];
					double trueKinship = .25*pathToOmega.get(trueKey)[1] + .5*pathToOmega.get(trueKey)[2];

					
					distMat[ind1][ind2] = trueKinship;
	
								
				}
			}
			
			
			/////////////////////////////////////////////////////////////
			// analyze networks
			
			//get number of networks
			reader = DataParser.openReader(resultPath + String.format("Summary_%s.%d.genome.txt", testName, t));
			reader.readLine();
			int nNetworks = Integer.parseInt(reader.readLine().split(": ")[1]);
			reader.close();
			
			//for each network
			for(int net=1; net<=nNetworks; net++){
				
				//network directory
				String netDir = resultPath + String.format("%s.%d.genome_network%d/", testName, t, net);
				
				//get num max score pedigrees and indiv names
				reader = DataParser.openReader(netDir + String.format("Summary_%s.%d.genome_network%d.txt", testName, t, net));
				int nPed = 1;
				List<String> names = new ArrayList<String>();
				
				while((line = reader.readLine())!=null){
					String[] fields = line.split("\\s");
					
					if(fields.length < 2) break;
					
					
					if(fields[1].equals("Num_at_max_score:")) 
						nPed = Integer.parseInt(fields[2]);
					
					else if(fields[1].equals("Sample_IIDs:")){
						String[] inds = fields[2].split(",");
						for(String name : inds){
							names.add(name);
						}
					}	
				}
				reader.close();
				
				
				//score unrelated between this network and others
				for(String name1 : names){
					
					int rid1 = Integer.parseInt(name1.split("__")[1]) - 1;
					
					for(int rid2=0; rid2<20; rid2++){
						
						if(names.contains(String.format("%d__%d", 1, rid2+1)) || rid1==rid2) continue;

						//acc
						int acc = 0;
						if(trueRel[rid1][rid2].getNumVisit()==0) acc = 1;		
						accMat[rid1][rid2] = acc;
						
						//distance
						Path trueKey = trueRel[rid1][rid2];
						double trueKinship = .25*pathToOmega.get(trueKey)[1] + .5*pathToOmega.get(trueKey)[2];
						distMat[rid1][rid2] = trueKinship;
						
						
					}
					
					
					
				}
				
				
				//get number of good pedigrees
				int numGoodPed = nPed;
				for(int p=1; p<=nPed; p++){
					
					String famDir = netDir + String.format("%s.%d.genome_network%d_%d.fam",testName,t,net,p);
					
					//make name2Index map
					Map<String, Integer> name2Index = new HashMap<String, Integer>();
					reader = DataParser.openReader(famDir);
					
					int idx = 0;
					while((line = reader.readLine())!=null){
						String name = line.split("\t")[1];
						name2Index.put(name, idx++);
					}
					reader.close();
					
					
					//build pedigree
					Pedigree ped = new Pedigree(famDir, name2Index);
					
					//skip if pedigree is looped
					if(ped.looped){
						numGoodPed--;
					}
				}
				
				if(numGoodPed==0) continue;
				
				
				//for each max score pedigree
				for(int p=1; p<=nPed; p++){
					
					String famDir = netDir + String.format("%s.%d.genome_network%d_%d.fam",testName,t,net,p);
					
					//make name2Index map
					Map<String, Integer> name2Index = new HashMap<String, Integer>();
					reader = DataParser.openReader(famDir);
					
					int idx = 0;
					while((line = reader.readLine())!=null){
						String name = line.split("\t")[1];
						name2Index.put(name, idx++);
					}
					reader.close();
					
					
					//build pedigree
					Pedigree ped = new Pedigree(famDir, name2Index);
					
					//skip if pedigree is looped
					if(ped.looped){
						continue;
					}
					
					Path[][] rel = ped.getRelationships();
					
					//score pairwise accuracy within network
					for(String name1 : names){
						for(String name2 : names){
							
							int id1 = name2Index.get(name1);
							int id2 = name2Index.get(name2);
							
							if(id1 >= id2) continue;
							
							Path inferred = rel[id1][id2];
							
							//real ID
							int rid1 = Integer.parseInt(name1.split("__")[1]) - 1;
							int rid2 = Integer.parseInt(name2.split("__")[1]) - 1;
							Path truth = trueRel[rid1][rid2];

							//acc
							if(Accuracy.hasSameIBD(pathToOmega, inferred, truth))
								accMat[rid1][rid2] += 1.0/numGoodPed;							

							
							//dist
							double trueKinship = .25*pathToOmega.get(truth)[1] + .5*pathToOmega.get(truth)[2];
							double inferredKinship = .25*pathToOmega.get(inferred)[1] + .5*pathToOmega.get(inferred)[2];
							distMat[rid1][rid2] += Math.abs(trueKinship - inferredKinship)/numGoodPed;
							
							
							
						}
						
					}
					
					
					
				}
				
											

				
			}
			
			//////////////////////////////////////
			//Kinship accuracy
			//write header for map file
			accWriter.write(String.format(">\t%d\n", t));
			
			for(int i=0; i<accMat.length; i++){
				for(int j=i+1; j<accMat.length; j++){
					accWriter.write(String.format("%d\t%d\t%f\n", i, j, accMat[i][j]));
				}
			}
			
			//kinship distance
			distWriter.write(String.format(">\t%d\n", t));
			
			for(int i=0; i<accMat.length; i++){
				for(int j=i+1; j<accMat.length; j++){
					distWriter.write(String.format("%d\t%d\t%.8f\n", i, j, distMat[i][j]));
				}
			}
			
			

		}
		
		accWriter.close();
		distWriter.close();
		
		
	}
	

}

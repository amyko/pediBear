package test;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
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
		
		int numIndiv = 18;
		
		//directories
		String testName = "sim4";
		String dir = System.getProperty("user.home") + "/Google Drive/Research/pediBear/data/simulations/";
		String outPath = dir + "results/sim4.primus.relate";
		String truePath = dir + "results/sim4.true";
		
		//open files
		PrintWriter accWriter = DataParser.openWriter(outPath+".mapAcc");
		PrintWriter distWriter = DataParser.openWriter(outPath+".kinshipDist");
		PrintWriter nWriter = DataParser.openWriter(outPath+".nCorrect");
		Path[][] trueRel = Accuracy.getTruePath(truePath, numIndiv);
		Map<Path, double[]> pathToOmega = Accuracy.getPathToOmega(dir + "pathToOmega.txt");
		
		for(int t=0; t<100; t++){
			
			//System.out.println(t);
			
			String resultPath = dir + String.format("sim4_primus_relate/result.%d/", t);
			
			//if didn't finish running, skip
			BufferedReader mySummary;
			try{
				mySummary = DataParser.openReader(String.format(resultPath + "Summary_sim4.%d.relate.txt",t));
			}
			catch(FileNotFoundException e){
				System.out.println(String.format("Did not complete: %d", t));
				continue;
			}
			
			boolean hasError = false;
			String myLine;
			//if there's an error message or no pedigree was constructed, skip
			while((myLine=mySummary.readLine())!=null){
				String[] fields = myLine.split("\t");
				
				if(fields.length>=7 && !fields[0].equals("Network_num")){
					if(!fields[7].equals("none") || Integer.parseInt(fields[2])==0){
						hasError = true;
					}
				}

			}
			
			if(hasError){
				System.out.println(String.format("Error: %d", t));
				continue;
			}
			
			
			//accuracy matrix (i,j,acc)
			double[][] accMat = new double[numIndiv][numIndiv];
			double[][] distMat = new double[numIndiv][numIndiv];
			
			//////////////////////////////////////////////////////////////
			//write accuracy for unrelated individuals
			BufferedReader reader = DataParser.openReader(resultPath + String.format("%s.%d.relate_unrelated_samples.txt", testName, t));
			reader.readLine();
			
			String line;
			List<Integer> unrel = new ArrayList<Integer>();
			while((line = reader.readLine())!=null){
				String[] fields = line.split("\t");
				unrel.add(Integer.parseInt(fields[0]) - 1);
			}
			reader.close();
			
			//System.out.println("Unrelated: ");
			
			for(int ind1 : unrel){
				for(int ind2 : unrel){
		
					if(ind1>=ind2) continue;
					
					//System.out.println(String.format("%d %d", ind1, ind2 ));
					
					//accuracy
					int acc = 0;
					if(trueRel[ind1][ind2].getNumVisit()==0) acc = 1;
					accMat[ind1][ind2] = acc;
					accMat[ind2][ind1] = acc;
					
					//distance
					Path trueKey = trueRel[ind1][ind2];
					double trueKinship = .25*pathToOmega.get(trueKey)[1] + .5*pathToOmega.get(trueKey)[2];

					
					distMat[ind1][ind2] = trueKinship;
					distMat[ind2][ind1] = trueKinship;
	
								
				}
			}
			
			
			/////////////////////////////////////////////////////////////
			// analyze networks
			
			//get number of networks
			reader = DataParser.openReader(resultPath + String.format("Summary_%s.%d.relate.txt", testName, t));
			reader.readLine();
			int nNetworks = Integer.parseInt(reader.readLine().split(": ")[1]);
			reader.close();
			int nCorrect = 0;
			int nReported = 0;
			
			//for each network
			for(int net=1; net<=nNetworks; net++){
				
				
				//network directory
				String netDir = resultPath + String.format("%s.%d.relate_network%d/", testName, t, net);
				
				
				//make ped2rank map
				Map<Integer, Integer> pednum2rank = new HashMap<Integer,Integer>();
				
								
				//get ped2rank and indiv names
				int nTotalPed = 0;
				reader = DataParser.openReader(netDir + String.format("Summary_%s.%d.relate_network%d.txt", testName, t, net));
				List<String> names = new ArrayList<String>();
				boolean passedHeader = false;
				
				while((line = reader.readLine())!=null){
					String[] fields = line.split("\\s");
					
					if(fields.length < 2){
						passedHeader=true;
						reader.readLine();
						continue;
					}
					
					if(fields[1].equals("Num_pedigrees:")){
						String[] fields2 = fields[2].split(",");
						nTotalPed = Integer.parseInt(fields2[0]);
					}	

					else if(fields[1].equals("Sample_IIDs:")){
						String[] inds = fields[2].split(",");
						for(String name : inds){
							names.add(name);
						}
					}
					
					
					
					//ped line
					if(passedHeader){
						pednum2rank.put(Integer.parseInt(fields[0]), Integer.parseInt(fields[3]));
					}
					
				}
				
				reader.close();
				
				
				//score unrelated between this network and others
				//System.out.println("Between networks: ");
				
				for(String name1 : names){
					
					//int rid1 = Integer.parseInt(name1.split("__")[1]) - 1;
					int rid1 = Integer.parseInt(name1) - 1;
					
					for(int rid2=0; rid2<numIndiv; rid2++){
						
						//if(names.contains(String.format("%d__%d", 1, rid2+1)) || rid1==rid2) continue;
						if(names.contains(""+(rid2+1)) || rid1==rid2) continue;
						
						//System.out.println(String.format("%d %d", rid1, rid2));

						//acc
						double acc = 0;
						if(trueRel[rid1][rid2].getNumVisit()==0) acc = 1;		
						accMat[rid1][rid2] = acc;
						accMat[rid2][rid1] = acc;
						
						//distance
						Path trueKey = trueRel[rid1][rid2];
						double trueKinship = .25*pathToOmega.get(trueKey)[1] + .5*pathToOmega.get(trueKey)[2];
						distMat[rid1][rid2] = trueKinship;
						distMat[rid2][rid1] = trueKinship;
						
					}
					
					
					
				}
				
				
				//get which rank to look at, and how many there are
				int currRank = 1;
				int currNum = 1;
				int maxRank = pednum2rank.get(nTotalPed);
				int numGoodPed = 0;
				
				
				while(currRank <= maxRank){
								
					for(int p=currNum; p<=nTotalPed; p++){
						
						//break if the rank exceeds currRank
						if(pednum2rank.get(p) > currRank){
							currNum = p;
							break;
						}
						
						String famDir = netDir + String.format("%s.%d.relate_network%d_%d.fam",testName,t,net,p);
						
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
						if(!ped.looped){
							numGoodPed++;
						}
						
						
					}
					
					//if >0 outbred pedigrees exist, we know which rank and number to look at
					if(numGoodPed>0){
						break;
					}
					else currRank++;
		
					
				}
				
			
				
				if(numGoodPed==0){
					System.out.println(String.format("No outbred pedigrees: %d", t));
					nCorrect = 0;
					nReported = 0;
					break;
				}
				
				//System.out.println("Within networks: ");
				
				//for each max score pedigree (within network)
				for(int p=1; p<=nTotalPed; p++){
					
					
					int nCorrectCurr = 1;
					
					//break if the rank exceeds currRank
					if(pednum2rank.get(p) < currRank){
						continue;
					}
					if(pednum2rank.get(p) > currRank){
						break;
					}
					
					String famDir = netDir + String.format("%s.%d.relate_network%d_%d.fam",testName,t,net,p);
					
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
					
					
					nReported++;
							
					Path[][] rel = ped.getRelationships();
					
					//score pairwise accuracy within network
					for(String name1 : names){
						for(String name2 : names){
							
							int id1 = name2Index.get(name1);
							int id2 = name2Index.get(name2);
							
							if(id1 >= id2) continue;
							
							if(p==1){
								//System.out.println(String.format("%s %s", name1, name2));
							}
							
							
							Path inferred = rel[id1][id2];
							
							//System.out.println(String.format("%d %d %d %d %d", id1, id2, inferred.getUp(), inferred.getDown(), inferred.getNumVisit()));
							
							//real ID
							//int rid1 = Integer.parseInt(name1.split("__")[1]) - 1;
							//int rid2 = Integer.parseInt(name2.split("__")[1]) - 1;
							int rid1 = Integer.parseInt(name1) - 1;
							int rid2 = Integer.parseInt(name2) - 1;
							Path truth = trueRel[rid1][rid2];

							
							//correct?
							if(inferred.getUp()!=truth.getUp() || inferred.getDown()!=truth.getDown() || inferred.getNumVisit()!=truth.getNumVisit()){
								nCorrectCurr *= 0;
							}
							else{
								nCorrectCurr *= 1;
							}
							

							
							if(Accuracy.hasSameIBD(pathToOmega, inferred, truth)){
								accMat[rid1][rid2] += 1.0/numGoodPed;			
								accMat[rid2][rid1] += 1.0/numGoodPed;
							}

							
							//dist
							double trueKinship = .25*pathToOmega.get(truth)[1] + .5*pathToOmega.get(truth)[2];
							
							if(inferred.getNumVisit()==-1){
								distMat[rid1][rid2] = -1;
								distMat[rid2][rid1] = -1;
								System.out.println("Shouldn't be here");
							}
							
							else{
								double inferredKinship = .25*pathToOmega.get(inferred)[1] + .5*pathToOmega.get(inferred)[2];
								distMat[rid1][rid2] += Math.abs(trueKinship - inferredKinship)/numGoodPed;
								distMat[rid2][rid1] += Math.abs(trueKinship - inferredKinship)/numGoodPed;
							}
							
							
							
						}
						
					}
					
					
					//System.out.println(nCorrectCurr);
				
					if(nCorrectCurr==1){
						nCorrect = 1;
					}
				
				}
												

				
			}
			
			if(nReported==0) continue;
			
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
			
			nWriter.write(String.format("%d\t%d\t%d\n", nReported, nCorrect, t));
			
			

		}
		
		accWriter.close();
		distWriter.close();
		nWriter.close();
		
		
	}
	

}

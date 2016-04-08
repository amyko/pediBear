package Unused;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import utility.ArrayUtility;
import utility.DataParser;
import dataStructures.EmissionKey;
import dataStructures.Genotype;
import dataStructures.Node;
import dataStructures.Relationship;
import dataStructures.SimplePair;

public class UnusedMethods {
	 	//pairwise likelihood for dependent sites
	public double[][][] forwardAlgorithm(String genoPath, String infoPath, int[] indCols, List<Relationship> rel) throws IOException{
		
		int numIndiv = indCols.length;
		int numRel = rel.size();
		double[][][] toReturn = new double[numRel][numIndiv][numIndiv];
		int numSnp = DataParser.countLines(infoPath)-1;

		
		//open file
		BufferedReader genoFile = DataParser.openReader(genoPath);
		BufferedReader infoFile = DataParser.openReader(infoPath);
		genoFile.readLine(); //skip header
		infoFile.readLine(); //skip header
		
		// data structures
		Map<Integer,String[]> posToGeno = new HashMap<Integer,String[]>(back); //holds genotypes for the last 50 snps
		Map<Integer, String[]> posToInfo = new HashMap<Integer, String[]>(back); //holds info line for the last 50 snps

		//previous forward probabilities (only 3 distinct ibd states)
		double[][][] alpha = new double[numIndiv][numIndiv][3];
		
		//initialization:
		//read file
		String[] geno = genoFile.readLine().split("\t");
		String[] info = infoFile.readLine().split("\t");
	
		Map<EmissionKey, Double> emissionMap = computePossibleOneLocusEmissionWithError(info);
		
		//initial alpha
		for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
			for(int i=0; i<numIndiv; i++){
				int i1 = indCols[i];
				for(int j=i+1; j<numIndiv; j++){
					int i2 = indCols[j];
					EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(genotypeKey.get(geno[i1]), genotypeKey.get(geno[i2])), null, ibd);
					alpha[i][j][ibd] = rel.getMarginalProbs()[ibd] * emissionMap.get(key);
				}
			}
		}

		// scale alphas
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				if(ArrayUtility.sum(alpha[i][j]) == 0d){
					toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
				}
				
				double scalingFactor = 1d / ArrayUtility.sum(alpha[i][j]);

				for(int ibd = 0; ibd < alpha[0][0].length; ibd++) alpha[i][j][ibd] = alpha[i][j][ibd] * scalingFactor;//scale alpha
				
				//add to lkhd
				toReturn[i][j] -= Math.log(scalingFactor);

			}
		}
		
	
		//update previous data
		int prevPos = Integer.parseInt(info[POS]);
		posToGeno.put(prevPos,geno);
		posToInfo.put(prevPos, info);
		

		//recursion:
		for(int locus = 1; locus < numSnp; locus++){
			//if (locus%10000==0) System.out.println(locus);
			
			double[][][] tempAlpha = new double[numIndiv][numIndiv][3];
			
			//read info
			info = infoFile.readLine().split("\t");
			int ldPos = Integer.parseInt(info[LDPOS]);
			String[] ldInfo = posToInfo.get(ldPos);
			
			//read geno
			geno = genoFile.readLine().split("\t");
			String[] ldGeno = posToGeno.get(ldPos);
			
			//read relevant data
			int currPos = Integer.parseInt(info[POS]);
			double dist = (currPos - prevPos) * recombRate;

			
			//compute all possible emission probs
			Map<EmissionKey, Double> conditionalEmissionMap = computePossibleConditionalEmissionWithError(ldInfo, info);
			
			
			for (int ibdNew = 0; ibdNew < alpha[0][0].length; ibdNew++){
				//transition
				for(int ibdPrev = 0; ibdPrev < alpha[0][0].length; ibdPrev++){	
					
					for(int i=0; i<numIndiv; i++){
						for(int j=i+1; j<numIndiv; j++){
							double transitionDensity = transitionDensity(ibdNew, ibdPrev, dist, rel);
							tempAlpha[i][j][ibdNew] += alpha[i][j][ibdPrev] * transitionDensity;
						}
					}

				}
				
				
				//emission	
				for(int i=0; i<numIndiv; i++){
					int i1 = indCols[i];
					for(int j=i+1; j<numIndiv; j++){
						int i2 = indCols[j];
						EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(new Genotype(ldGeno[i1]), new Genotype(ldGeno[i2])), new SimplePair<Genotype,Genotype>(new Genotype(geno[i1]), new Genotype(geno[i2])), ibdNew);
						double emissionDensity = conditionalEmissionMap.get(key); //should never be zero
						//double emissionDensity = emissionDensityWithError(info, geno[i1], geno[i2], ibdNew); //should never be zero		
						tempAlpha[i][j][ibdNew] *= emissionDensity;
						
						
					}
				}
				

			}
			
			
			//scale alphas
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					if(ArrayUtility.sum(tempAlpha[i][j]) == 0d){
						toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
					}
					
					double scalingFactor = 1d / ArrayUtility.sum(tempAlpha[i][j]);

					for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
						alpha[i][j][ibd] = tempAlpha[i][j][ibd] * scalingFactor;//scale alpha
					}
					
					//add to lkhd
					toReturn[i][j] -= Math.log(scalingFactor);
					
				}
			}

			
			//add previous info
			prevPos = currPos;
			posToGeno.put(prevPos,geno);
			posToInfo.put(prevPos, info);
			

			//remove useless
			if(posToGeno.size()>back){
				int minKey = Collections.min(posToGeno.keySet());
				posToGeno.remove(minKey);
				posToInfo.remove(minKey);
			}
			
		}
		
		return toReturn;
		
	}	
	 
	
	//pairwise likelihood for independent sites
	public double[][] forwardAlgorithmIndep(String genoPath, String infoPath, int[] indCols, Relationship rel) throws IOException{
		
		int numIndiv = indCols.length;
		double[][] toReturn = new double[numIndiv][numIndiv];
		int numSnp = DataParser.countLines(infoPath)-1;

		
		//open file
		BufferedReader genoFile = DataParser.openReader(genoPath);
		BufferedReader infoFile = DataParser.openReader(infoPath);
		genoFile.readLine(); //skip header
		infoFile.readLine(); //skip header
		

		//previous forward probabilities (only 3 distinct ibd states)
		double[][][] alpha = new double[numIndiv][numIndiv][3];
		
		//initialization:
		//read file
		String[] geno = genoFile.readLine().split("\t");
		String[] info = infoFile.readLine().split("\t");
	
		Map<EmissionKey, Double> emissionMap = computePossibleOneLocusEmissionWithError(info);
		
		//initial alpha
		for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
			for(int i=0; i<numIndiv; i++){
				int i1 = indCols[i];
				for(int j=i+1; j<numIndiv; j++){
					int i2 = indCols[j];
					EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(genotypeKey.get(geno[i1]), genotypeKey.get(geno[i2])), null, ibd);
					alpha[i][j][ibd] = rel.getMarginalProbs()[ibd] * emissionMap.get(key);
				}
			}
		}

		// scale alphas
		for(int i=0; i<numIndiv; i++){
			for(int j=i+1; j<numIndiv; j++){
				if(ArrayUtility.sum(alpha[i][j]) == 0d){
					toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
				}
				
				double scalingFactor = 1d / ArrayUtility.sum(alpha[i][j]);

				for(int ibd = 0; ibd < alpha[0][0].length; ibd++) alpha[i][j][ibd] = alpha[i][j][ibd] * scalingFactor;//scale alpha
				
				//add to lkhd
				toReturn[i][j] -= Math.log(scalingFactor);

			}
		}
		
		int prevPos = Integer.parseInt(info[POS]);
		

		//recursion:
		for(int locus = 1; locus < numSnp; locus++){

			
			double[][][] tempAlpha = new double[numIndiv][numIndiv][3];
			
			//read info
			info = infoFile.readLine().split("\t");
			geno = genoFile.readLine().split("\t");

			
			//read relevant data
			int currPos = Integer.parseInt(info[POS]);
			double dist = (currPos - prevPos) * recombRate;

			
			//compute all possible emission probs
			emissionMap = computePossibleOneLocusEmissionWithError(info);
			
			
			for (int ibdNew = 0; ibdNew < alpha[0][0].length; ibdNew++){
				//transition
				for(int ibdPrev = 0; ibdPrev < alpha[0][0].length; ibdPrev++){	
					
					for(int i=0; i<numIndiv; i++){
						for(int j=i+1; j<numIndiv; j++){
							double transitionDensity = transitionDensity(ibdNew, ibdPrev, dist, rel);
							tempAlpha[i][j][ibdNew] += alpha[i][j][ibdPrev] * transitionDensity;
						}
					}

				}
				
				
				//emission
				for(int i=0; i<numIndiv; i++){
					int i1 = indCols[i];
					for(int j=i+1; j<numIndiv; j++){
						int i2 = indCols[j];
						EmissionKey key = new EmissionKey(new SimplePair<Genotype,Genotype>(genotypeKey.get(geno[i1]), genotypeKey.get(geno[i2])), null, ibdNew);
						tempAlpha[i][j][ibdNew] *= emissionMap.get(key);

					}
				}
				

				

			}
			
			
			//scale alphas
			for(int i=0; i<numIndiv; i++){
				for(int j=i+1; j<numIndiv; j++){
					if(ArrayUtility.sum(tempAlpha[i][j]) == 0d){
						toReturn[i][j] = Double.NEGATIVE_INFINITY;	//It is impossible for these genotypes to have this relationship (should never happen with sequencing errors)
					}
					
					double scalingFactor = 1d / ArrayUtility.sum(tempAlpha[i][j]);

					for(int ibd = 0; ibd < alpha[0][0].length; ibd++){
						alpha[i][j][ibd] = tempAlpha[i][j][ibd] * scalingFactor;//scale alpha
					}
					
					//add to lkhd
					toReturn[i][j] -= Math.log(scalingFactor);
					
				}
			}
			
			prevPos = currPos;


		}
		
		return toReturn;
	}
	
	
	
	public Node getLastExistingNode(Node node){
		
		clearVisit();
		nDeleted = 0;
		return lastExistingNode(node);
		
	}
	
	//returns the last node visited that was NOT deleted, and the number of deleted nodes
	//clear visit before calling this method
	private Node lastExistingNode(Node node){//works

		//mark visit
		node.numVisit++;
		
		if(node.sampled || node.getNumUnvisitedNeighbors() > 1){
			return node;
		}

		else{

			nDeleted++;
			
			//recurse
			for(Node i : node.getParents()){
				if(i.numVisit > 0) 
					continue;
				else 
					return lastExistingNode(i);
			}
			
			for(Node i : node.getChildren()){
				if(i.numVisit > 0) 
					continue;
				else 
					return lastExistingNode(i);
			}
			
		}
		
		throw new RuntimeException("last existing node not found!");

	}
	
	
	
	
	
	public void stretch(Node parent){
		
		//get cluster
		clearVisit();
		List<Node> ped = parent.getConnectedNodes(new ArrayList<Node>());

		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);	

		
		//stretch
		List<Node> ghostParents = new ArrayList<Node>(parent.getChildren().size());
	
		for(Node child : parent.getChildren()){
			Node ghostParent = makeNewNode(parent.getDepth(), parent.getSex());
			ghostParents.add(ghostParent);
			
			child.removeParent(parent);
			child.addParent(ghostParent);
			ghostParent.addChild(child);
			ghostParent.addParent(parent);
		}
		
		parent.setDepth(parent.getDepth()+1);
		parent.setChildren(ghostParents);
		
		
		//update adj matrix 
		for(Node ind : inds){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
	
		
	}
	
	
	public void compress(Node grandParent){
		
		//get cluster
		clearVisit();
		List<Node> ped = grandParent.getConnectedNodes(new ArrayList<Node>());

		//subtract old terms
		this.logLikelihood[curr] -= likelihoodLocalPedigree(ped);	
		
		
		//compress
		//get grand children
		List<Node> grandChildren = new ArrayList<Node>();
		List<Node> ghostParents = new ArrayList<Node>(grandParent.getChildren().size());
		for(Node ghostParent : grandParent.getChildren()){
			
			ghostParents.add(ghostParent);
			
			//update every grand child
			for(Node grandChild : ghostParent.getChildren()){
				
				grandChildren.add(grandChild);
			
			}
			
		}
		
		//delete ghost parents
		for(Node gc : ghostParents){
			deleteNode(gc);
		}

		
		//update edges
		for(Node gc : grandChildren){
			gc.addParent(grandParent);
			grandParent.addChild(gc);
		}
		grandParent.setDepth(grandParent.getDepth()-1);
		

		
		//update adj matrix 
		for(Node ind : inds){
			updateAdjMat(ind);
		}


		//add new terms
		this.logLikelihood[curr] += likelihoodLocalPedigree(ped);
		
		
	}

	
	
}

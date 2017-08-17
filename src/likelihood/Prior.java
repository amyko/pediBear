package likelihood;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import dataStructures.Node;
import utility.ArrayUtility;

//gasbarra prior
public class Prior {
	
	
	//gasbarra prior
	private double alpha = 0;
	private double beta = 0;
	private int[] F;
	private int[] M;
	private int[] Fs;
	private int[] Ms;
	private Map<Integer, Integer> Cf;
	private Map<Integer, HashMap<Integer, Integer>> Cfm;
	private int[] nFSampled;
	private int[] nMSampled;
	private int[] kSoFar;
	Set<Integer> oldF = new HashSet<Integer>();
	Set<Integer> oldM = new HashSet<Integer>();
	
	private Random rGen;
	private int minN;
	private int maxN;
	private int stepSize;
	private int maxDepth;
	
	
	public Prior(Random rGen, int maxDepth, int minN, int maxN, int stepSize) {
				
		//prior 
		F = new int[maxDepth+1];
		M = new int[maxDepth+1];
		Fs = new int[maxDepth+1];
		Ms = new int[maxDepth+1];
		Cf = new HashMap<Integer,Integer>();
		Cfm = new HashMap<Integer, HashMap<Integer, Integer>>();
		nFSampled = new int[maxDepth+1];
		nMSampled = new int[maxDepth+1];
		kSoFar = new int[maxDepth+1];
		
		this.rGen = rGen;
		this.minN = minN;
		this.maxN = maxN;
		this.stepSize = stepSize;
		this.maxDepth = maxDepth;
		
		
	}
	
	
	
	
	public double computePrior(List<Node> nodes) {
		
		//hyper parameters
		int N = rGen.nextInt((maxN - minN)/stepSize)*stepSize + minN;
		int Nf = N/2;
		int Nm = Nf;
		beta = rGen.nextDouble(); //TODO limit
		alpha = rGen.nextDouble();
		
		//clear
		ArrayUtility.clear(F);
		ArrayUtility.clear(M);
		ArrayUtility.clear(Fs);
		ArrayUtility.clear(Ms);
		Cf.clear();
		Cfm.clear();
		ArrayUtility.clear(nFSampled);
		ArrayUtility.clear(nMSampled);
		countSampledNodes(nodes, nFSampled, nMSampled);
		ArrayUtility.clear(kSoFar);
		
		
		Set<Integer> oldF = new HashSet<Integer>();
		Set<Integer> oldM = new HashSet<Integer>();
		
		double totalProb = 0;
		
		//for every node
		for(Node child : nodes) {
			
			
			double childProb = getProbForChild(child, N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta);
			
			totalProb += childProb;
			
		}
		
		
		
		return totalProb;
		
	}
	
	
	
	private double getProbForChild(Node child, int N, int Nf, int Nm,  int[] F, int[] M, int[] Fs, int[] Ms, int[] kSoFar, int[] nFSampled, int[] nMSampled, Map<Integer, Integer> Cf, Map<Integer, HashMap<Integer, Integer>> Cfm, Set<Integer> oldF, Set<Integer> oldM, double alpha, double beta) {
		
		if(child.getDepth()==maxDepth) return 0;
		
		double toReturn = 0;
	
		Node father = child.getParentWithSex(1);
		Node mother = child.getParentWithSex(0);
		int parentDepth = child.getDepth()+1;
		int k = kSoFar[child.getDepth()];
		int f = father==null? -1 : father.getIndex();
		int m = mother==null? -1 : mother.getIndex();
		
		
		//father assignment
		if(oldF.contains(f))  //choose old father
			toReturn += Math.log(alpha + Cf.get(f)) - Math.log(Nf*alpha + k);
		
		else {//new father

			//unsampled (including ghost father)
			if(f==-1 || !father.sampled)
				toReturn += 	Math.log(alpha) + Math.log(Nf - F[parentDepth] - nFSampled[parentDepth] + Fs[parentDepth]) - Math.log(Nf*alpha + k);
			
			//sampled
			else{
				toReturn += Math.log(alpha) - Math.log(Nf*alpha + k);
				Fs[parentDepth]++;
			}
			
			//update father count
			F[parentDepth]++;


		}
		

	
		
		if(oldM.contains(m)) {//old mother
			
			//new father
			if(!oldF.contains(f))
				toReturn += -Math.log(Nm);
			else
				toReturn += Math.log(beta + Cfm.get(f).get(m)) - Math.log(Nm*beta + Cf.get(f));
			
		}
		
		else {//new mother
		
			if(oldF.contains(f)) { //old father
				
				//unsampled
				if(m==-1 || !mother.sampled)
					toReturn += Math.log(beta) + Math.log(Nm - M[parentDepth]- nMSampled[parentDepth] + Ms[parentDepth]) - Math.log(Nm*beta + Cf.get(f));
				
				
				//sampled 
				else {
					toReturn += Math.log(beta) - Math.log(Nm*beta + Cf.get(f));
					Ms[parentDepth]++;
				}
				
			}
			
			else {//new father
				
				//unsampled
				if(m==-1 || !mother.sampled)
					toReturn += Math.log(Nm - M[parentDepth]- nMSampled[parentDepth] + Ms[parentDepth]) - Math.log(Nm);
				
				//sampled 
				else {
					toReturn += -Math.log(Nm);
					Ms[parentDepth]++;
				}
				
			}
			
			//update mother count
			M[parentDepth]++;

			
		}
	
		
		
		//update counts
		kSoFar[child.getDepth()]++;
		if(f!=-1) {
			if(!Cf.containsKey(f)) Cf.put(f, 1);
			else Cf.put(f, Cf.get(f)+1);
			oldF.add(f);
		}
		
		if(m!=-1)
			oldM.add(m);
		
		if(f!=-1 && m!=-1) {
			
			if(!Cfm.containsKey(f))
				Cfm.put(f, new HashMap<Integer, Integer>());	
			if(!Cfm.get(f).containsKey(m))
				Cfm.get(f).put(m, 1);
			else
				Cfm.get(f).put(m, Cfm.get(f).get(m)+1);
		}
		
		
		//recurse on ghost parents
		if(f==-1)
			toReturn += getProbForChild(new Node("ghost", "ghost", 0, false, 0, child.getDepth()+1, -1), N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta);
		if(m==-1)
			toReturn += getProbForChild(new Node("ghost", "ghost", 0, false, 0, child.getDepth()+1, -1), N, Nf, Nm, F, M, Fs, Ms, kSoFar, nFSampled, nMSampled, Cf, Cfm, oldF, oldM, alpha, beta);
		
		
		return toReturn;
		
		
	}
	
	
	private static void countSampledNodes(List<Node> nodes, int[] nFSampled, int[] nMSampled) {
		

		for(Node x : nodes) {
			
			//count sampled nodes
			if(x.sampled) {
				
				if(x.getSex()==0)
					nMSampled[x.getDepth()]++;
				
				else
					nFSampled[x.getDepth()]++;
				
			}
				

			
			
		}	

	}
	
	
	/*
	//TODO sample depth
	private double sampleDepthPrior(int[] nSampledMoms, int[] nSampledDads){
		
		double toReturn = 0d;
		
		// P(max) = 1/totalUnits, P(max-1) = 2/totalUnits, ..., P(0)=(maxDepth+1)/totalUnits
		for(int i=0; i<nSampledMoms.length; i++){

			toReturn += (nSampledMoms[i] + nSampledDads[i]) * (Math.log(maxDepth+1-i));
			
		}
		
		toReturn = toReturn - numIndiv*Math.log(totalUnits);
		
		return toReturn;
	}
	*/
	

}

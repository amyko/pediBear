package likelihood;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import dataStructures.Node;
import dataStructures.Pedigree;
import utility.ArrayUtility;

//gasbarra prior
public class Prior {
	
	
	//gasbarra prior
	private int[] F;
	private int[] M;
	private int[] Fs;
	private int[] Ms;
	private Map<Integer, Integer> Cf; // number of children father f has so far
	private Map<Integer, HashMap<Integer, Integer>> Cfm; // number of children couple (f,m) has so far
	private int[] nFSampled;
	private int[] nMSampled;
	private int[] kSoFar;
	Set<Integer> oldF = new HashSet<Integer>();
	Set<Integer> oldM = new HashSet<Integer>();
	
	private int maxDepth;
	
	//ghost nodes
	Node ghostNode = new Node("ghost", "ghost", 0, false, 0, 0, -1);
	
	//sample prior
	double log2 = Math.log(2);


	//TODO when does this return NA? Implement alpha = inf, beta = inf
	public Prior(Random rGen, int maxDepth) {
				
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
		this.maxDepth = maxDepth;
	
		
	}
	
	

	
	
	
	
	public double computePrior(Pedigree currPedigree) {
		

		int N = currPedigree.getN();
		double alpha = currPedigree.getAlpha();
		double beta = currPedigree.getBeta();
		int nActiveNodes = currPedigree.getNActiveNodes();
		
		
		int Nf = N/2;
		int Nm = Nf;


		//clear
		ArrayUtility.clear(F);
		ArrayUtility.clear(M);
		ArrayUtility.clear(Fs);
		ArrayUtility.clear(Ms);
		Cf.clear();
		Cfm.clear();
		ArrayUtility.clear(nFSampled);
		ArrayUtility.clear(nMSampled);
		countSampledNodes(currPedigree, nFSampled, nMSampled);
		ArrayUtility.clear(kSoFar);
		oldF.clear();
		oldM.clear();
		

		double totalProb = 0;
		
		//for every node
		for(int i=0; i<nActiveNodes; i++) {
			
			Node child = currPedigree.getNode(i);
			
			double childProb = getProbForChild(child, N, Nf, Nm, alpha, beta);
			
			if(Double.isInfinite(childProb) || Double.isNaN(childProb)) {
				return Double.NEGATIVE_INFINITY;
			}
			
			
			totalProb += childProb;
			
		}
		

		//TODO sample depth?
		//return 0;
		return totalProb;
		
		
	}
	
	
	// returns log probability of this node choosing its parents
	private double getProbForChild(Node child, int N, int Nf, int Nm, double alpha, double beta) {	
		
		if(child.getDepth()>=maxDepth) return 0;
		
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
			else {
				
				//if (f,m) were never a couple
				if(!Cfm.containsKey(f) || !Cfm.get(f).containsKey(m)) {
					toReturn += Math.log(beta) - Math.log(Nm*beta + Cf.get(f));
				}
				
				else
					toReturn += Math.log(beta + Cfm.get(f).get(m)) - Math.log(Nm*beta + Cf.get(f));
			}
			
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
		if(f==-1) {
			ghostNode.setDepth(child.getDepth()+1);
			toReturn += getProbForChild(ghostNode, N, Nf, Nm, alpha, beta);
		}
		if(m==-1) {
			ghostNode.setDepth(child.getDepth()+1);
			toReturn += getProbForChild(ghostNode, N, Nf, Nm, alpha, beta);
		}
		
		return toReturn;
		
		
	}
	
	
	private static void countSampledNodes(Pedigree currPedigree, int[] nFSampled, int[] nMSampled) {
		

		for(int i=0; i<currPedigree.getNActiveNodes(); i++) {
			
			Node x = currPedigree.getNode(i);
			
			//count sampled nodes
			if(x.sampled) {
				
				if(x.getSex()==0)
					nMSampled[x.getDepth()]++;
				
				else
					nFSampled[x.getDepth()]++;
				
			}
				

			
			
		}	

	}
	
	
	
	//sampled depth = geometric (1/2)
	private double sampleDepthPrior(int[] nSampledMoms, int[] nSampledDads){
		
		double toReturn = 0d;
		
		// P(0) = 1/2, P(1) = 1/4, ...
		for(int i=0; i<nSampledMoms.length; i++){

			toReturn += -(nSampledMoms[i] + nSampledDads[i]) * (i+1) * log2;
			
		}
		
		return toReturn;
	}
	
	
	
	public static void main(String[] args) {
		
		
		// build pedigree
		Node c1 = new Node("fam","c1", 1, true, -1, 0, 0);
		Node c2 = new Node("fam", "c2", 0, true, -1, 0, 1);
		Node c3 = new Node("fam","c3", 1, true, -1, 0, 2);
		Node c4 = new Node("fam", "c4", 0, true, -1, 0, 3);
		Node c5 = new Node("fam","c5", 1, true, -1, 0, 4);
		Node c6 = new Node("fam", "c6", 0, true, -1, 0, 5);
		Node p1 = new Node("fam", "p1", 1, false, -1, 1, 6);
		Node p2 = new Node("fam", "p2", 0, false, -1, 1, 7);
		Node p3 = new Node("fam", "p3", 1, false, -1, 1, 8);

		
		p1.addChild(c1);
		p1.addChild(c2);
		c1.addParent(p1);
		c2.addParent(p1);
		
		p2.addChild(c3);
		p2.addChild(c4);
		c3.addParent(p2);
		c4.addParent(p2);
		
		p3.addChild(c5);
		p3.addChild(c6);
		c5.addParent(p3);
		c6.addParent(p3);
		
		 
		
		// mating parameters
		double alpha = .1;
		double beta = .01;
		int maxgen = 2;
		
		Random rgen = new Random();
		
		
		for(int N = 200; N < 201; N+=50) {
			
			int Nf = N/2;
			int Nm = Nf;
			
			Prior prior = new Prior(rgen, maxgen);
			double m1 = prior.getProbForChild(c1, N, Nf, Nm, alpha, beta);
			double m2 = prior.getProbForChild(c2, N, Nf, Nm, alpha, beta);
			double m3 = prior.getProbForChild(c3, N, Nf, Nm, alpha, beta);
			double m4 = prior.getProbForChild(c4, N, Nf, Nm, alpha, beta);
			double m5 = prior.getProbForChild(c5, N, Nf, Nm, alpha, beta);
			double m6 = prior.getProbForChild(c6, N, Nf, Nm, alpha, beta);
			double m7 = prior.getProbForChild(p1, N, Nf, Nm, alpha, beta);
			double m8 = prior.getProbForChild(p2, N, Nf, Nm, alpha, beta);
			double m9 = prior.getProbForChild(p3, N, Nf, Nm, alpha, beta);
			
			double prob = m1 + m2 + m3 + m4 + m5 + m6 + m7 + m8 + m9;
			
			System.out.print(String.format("%f, ", prob));
		}

		
		
		
	}
	

}

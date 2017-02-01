package mcmcMoves;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import dataStructures.Path;
import dataStructures.Pedigree;
import dataStructures.Node;

//generic MCMC move
public abstract class Move {

	public final static int maxDepth = 6;
	public final static int initSize = 10;

	
	//features for each move
	final int index; 
	public final String name;
	final double prob;
	
	//shared variables
	public static int indexCount;
	public final static int REJECT = 2;
	public final static double lnHalf = Math.log(.5);
	public final static double lnTwo = Math.log(2);
	private final static List<Double> powersOfHalf = new ArrayList<Double>();
	private final static List<Double> powersOfHalf2 = new ArrayList<Double>(); // 1/(2^i - 1)
	private final static List<Double> powersOfTwo = new ArrayList<Double>();
	private final static List<Double> logChooseTwo = new ArrayList<Double>();
	private final static List<Double> logChooseOne = new ArrayList<Double>();
	protected final static  Map<String, Double> moveProbs = new HashMap<String, Double>();
	
	public int nAccept = 0;
	public int nTried = 0;
	
	//misc
	List<Node> ancestors = new ArrayList<Node>();
	
	
	//store powers of twos and halves
	static{
		
		for(int i=0; i<initSize; i++){
			powersOfHalf.add(1d/Math.pow(2, i));
			powersOfHalf2.add(1d/(Math.pow(2, i) - 1)); // 1 / (2^d - 1)
			powersOfTwo.add(Math.pow(2, i));
			logChooseTwo.add(Math.log(2d/(i*(i-1))));
			logChooseOne.add(Math.log(1d/i));
		}
		
	}
	
	
	
	public Move(String name, double prob){
		
		this.index = indexCount++;
		this.name = name;
		this.prob = prob;
		moveProbs.put(name, prob);
		
	}

	//tries a move, tests it, and returns the old pedigree if it fails, and the new pedigree if it passes
	public final void mcmcMove(Pedigree currPedigree, double heat){

		
		
		double prevLkhd = currPedigree.getLogLikelihood();
		
		double acceptanceRatio = tryMove(currPedigree, heat);
		
		//banned moves
		if(acceptanceRatio==REJECT){
			//System.out.println("Skipped");
			return; //skip bad move
		}
		

		
		boolean reject = false;
		
		//check depth and age //TODO make this better
		for(int i=0; i<currPedigree.numIndiv; i++){
			
			Node myNode = currPedigree.getNode(i);
			
			//depth
			if(myNode.getDepth() > currPedigree.maxSampleDepth){
				reject = true;
				break;
			}
			
			
			//age
			if(myNode.getAge() != -1){
				ancestors.clear();
				myNode.getAncestors(ancestors);
				for(Node p : ancestors){
					
					if(p.getAge() != -1 && p.getAge() < myNode.getAge()){
						reject = true;
						break;
					}
					
				}
			}
			
			/*
			//ghost
			if(!myNode.sampled && myNode.getChildren().size()==0){
				
				System.out.println(this.name);
				reject = true;
				break;
				
			}
			*/
			
			
			
			
		}
		if(reject==true){
			reverseMove(currPedigree);
			return;
		}
		
		
		nTried++;
		double acceptanceAlpha = currPedigree.rGen.nextDouble();
		
		//reverse move
		if(acceptanceAlpha > acceptanceRatio){
			//System.out.println("REVERSED");
			reverseMove(currPedigree);
			
			if(Math.abs(prevLkhd - currPedigree.getLogLikelihood()) > 1){
				System.out.println("REVERSED NOT WORKING PROPERLY");
			}
			
		}
		
		//accept move and clean
		else{
			//System.out.println(this.name);	


			nAccept++;
			clean(currPedigree);
		}
			
	}

	
	//these should only ever be called by mcmcMove.  if they are called elsewhere, bad things might happen
	//(randomly) performs a proposed move, and returns the acceptance ratio for that move
	protected abstract double tryMove(Pedigree currPedigree, double heat);
	
	//reverts the last move that was looked at, changing both the graph structure and the adjacency matrix
	//in currPedigree
	protected abstract void reverseMove(Pedigree currPedigree);
	
	
	//deletes unnecessary ghost nodes
	protected abstract void clean(Pedigree currPedigree);


	
	//returns x ~ Geo(1/2); x=1,2, ...
	protected int geometricDist(Random rGen){

		double u = rGen.nextDouble();
		
		if(u==0) 
			return 1;
		else
			return (int) Math.ceil(Math.log(1-u) / lnHalf);		
		
	}
	
	protected double getPowersOfHalf(int idx){
		
		//grow list
		while(idx > powersOfHalf.size()-1){
			int oldCapacity = powersOfHalf.size();
			int newCapacity = (3*oldCapacity)/2 + 2;
			
			for(int i=oldCapacity; i<newCapacity; i++){
				powersOfHalf.add(1d/Math.pow(2, i));
			}	
		}
		
		return powersOfHalf.get(idx);
		
	}
	
	protected double getPowersOfHalf2(int idx){
		
		//grow list
		while(idx > powersOfHalf2.size()-1){
			int oldCapacity = powersOfHalf2.size();
			int newCapacity = (3*oldCapacity)/2 + 2;
			
			for(int i=oldCapacity; i<newCapacity; i++){
				powersOfHalf2.add(1d/(Math.pow(2, i)-1));
			}	
		}
		
		return powersOfHalf2.get(idx);
		
	}
	
	protected double getPowersOfTwo(int idx){
		
		//grow list
		while(idx > powersOfTwo.size()-1){
			int oldCapacity = powersOfTwo.size();
			int newCapacity = (3*oldCapacity)/2 + 2;
			
			for(int i=oldCapacity; i<newCapacity; i++){
				powersOfTwo.add(Math.pow(2, i));
			}	
		}
		
		return powersOfTwo.get(idx);
		
	}
	
	protected double getLogChooseOne(int idx){
		
		//grow list
		while(idx > logChooseOne.size()-1){
			int oldCapacity = logChooseOne.size();
			int newCapacity = (3*oldCapacity)/2 + 2;
			
			for(int i=oldCapacity; i<newCapacity; i++){
				logChooseOne.add(Math.log(1d/i));
			}	
		}
		
		return logChooseOne.get(idx);
		
	}
	
	protected double getLogChooseTwo(int idx){
		
		//grow list
		while(idx >= logChooseTwo.size()){
			int oldCapacity = logChooseTwo.size();
			int newCapacity = (3*oldCapacity)/2 + 2;
			
			//System.out.println(String.format("%d, %d", oldCapacity, newCapacity));
			
			for(int i=oldCapacity; i<newCapacity; i++){
				logChooseTwo.add(Math.log(2d/(i*(i-1))));
			}	
		}
		
		return logChooseTwo.get(idx);
		
	}
	
	
	public double getProb(){
		return this.prob;
	}
	
	
	
	
		
}
	


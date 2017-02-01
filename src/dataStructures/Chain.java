package dataStructures;

import java.util.ArrayList;
import java.util.List;

public class Chain {

	public final int index;
	private Pedigree ped;
	private double heat;
	private List<Integer> neighbors = new ArrayList<Integer>();

	
	
	public Chain(int index, Pedigree ped){
		
		this.index = index;
		this.ped = ped;
		
	}
	
	public void setPedigree(Pedigree ped){
		this.ped = ped;
	}
	
	public void setHeat(double deltaT){
		this.heat = 1d/(1+deltaT*index);
	}
	
	public void addNeighbor(int i){
		neighbors.add(i);
	}
	
	
	public Pedigree getPedigree(){
		return this.ped;
	}
	
	public double getLikelihood(){
		
		return this.ped.getLogLikelihood();
		
	}
	
	public double getHeat(){
		
		return this.heat;
		
	}
	
	public int getRandomNeighborIndex(){
		
		return neighbors.get(ped.rGen.nextInt(neighbors.size()));
		
	}
	
	
}


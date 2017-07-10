package dataStructures;


public class PedInfo{

	public double[] lkhd; //log lkhd P(G, N) for each N
	public int counts[]; //count for each N
	public int multiplier;

	
	public PedInfo(int minN, int maxN){
		
		//init lkhd
		lkhd = new double[maxN - minN + 1];
		for(int i=0; i<lkhd.length; i++){
			lkhd[i] = Double.NEGATIVE_INFINITY;
		}

		//init counts
		counts = new int[maxN - minN + 1]; 
		
	}
	

	
}

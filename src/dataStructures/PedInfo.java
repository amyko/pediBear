package dataStructures;


public class PedInfo{

	public double[] lkhd; //lkhd for each particular effective size
	public int count;
	public int multiplier;

	
	public PedInfo(int minN, int maxN, int stepSize){
		
		lkhd = new double[(maxN-minN)/stepSize];
		for(int i=0; i<lkhd.length; i++){
			lkhd[i] = 0;
		}
		
	}
	

	
}

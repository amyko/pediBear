package dataStructures;


public class PedInfo{

	public double lkhd;
	public int count;
	public int multiplier;
	public double[] logLkhdOfNe;

	
	
	public PedInfo(int minN, int maxN, int stepSize){
		
		logLkhdOfNe = new double[(maxN-minN)/stepSize];
		for(int i=0; i<logLkhdOfNe.length; i++){
			logLkhdOfNe[i] = 0;
		}
		
	}
	

	
}

package likelihood;
import utility.LogSum;

public class Prior {
	
	
	//returns stirling's numbers of the second kind
	public static double[][] logStirlingSecond(int N, int n, int g){
		
		int size = (int) Math.min(N, n*Math.pow(2, g-1)) + 1;
		
		double[][] toReturn = new double[size][size];
		
		//initialize
		for(int j=1; j<size; j++)
			toReturn[0][j] = Double.NEGATIVE_INFINITY;
		
		//recursion
		for(int i=1; i<size; i++){
			
			for(int j=0; j<size; j++){
				//first entry is always 0
				if(j==0)
					toReturn[i][j] = Double.NEGATIVE_INFINITY;
				
				else if(j<i)
					toReturn[i][j] = LogSum.addLogSummand(Math.log(j) + toReturn[i-1][j], toReturn[i-1][j-1]);
							
				else if(j==i)
					toReturn[i][j] = 0;
				else
					toReturn[i][j] = Double.NEGATIVE_INFINITY;
					
				
			}
		}

		
		return toReturn;
		
	}
	
	
	//returns falling log(n), log(n(n-1)), ..., log(n!)
	public static double[] logFallingFactorial(int n){
		
		double[] toReturn = new double[n];
		
		toReturn[0] = 0;
		toReturn[1] = Math.log(n);
		
		for(int i=2; i<n; i++){
			
			toReturn[i] = toReturn[i-1] + Math.log(n-i+1); 
		}
		
		return toReturn;
		
		
	}
	
	
	//returns 0!, 1!, ..., n!
	public static double[] logFactorial(int n){
		
		double[] toReturn = new double[n+1];
		
		toReturn[0] = 0;
		
		for(int i=1; i<toReturn.length; i++){
			
			toReturn[i] = toReturn[i-1] + Math.log(i);
			
		}
		
		return toReturn;
		
	}
	
	//returns matrix of ancestor probabilities; rows = generation, cols = number of ancestor
	public static double[][] ancestorLogProbs(int N, int n, int g){
		
		
		//constants
		double logN = Math.log(N/2);
		int upperBound = (int) Math.min(N, n*Math.pow(2,g));
		
		//initialize
		double[][] toReturn = new double[g+1][upperBound+1];
		for(int i=0; i<toReturn.length; i++){
			for(int j=0; j<toReturn[0].length; j++){
				toReturn[i][j] = Double.NEGATIVE_INFINITY;
			}
		}
		toReturn[0][n] = 0;
		
		//fill in each generation dynamically
		for(int i=1; i<=g; i++){
			
			//System.out.println(i);
			
			//max number of ancestors for this generation
			int jUpperBound = (int) Math.min(N, n*Math.pow(2,i));
			
			//pre-compute stirling's number, falling factorials
			double[][] logStirling = logStirlingSecond(N, n, g);
			double[] logFallingFact = logFallingFactorial(N/2);

						
			for(int j=2; j<=jUpperBound; j++){
				
				//boundaries for number of ancestors below this generation
				int kLower = (int) Math.ceil(j/2);
				int kUpper = (int) Math.min(N, n*Math.pow(2, g-1));
				
			
				//for every possible number of ancestors in g-1
				for(int k=kLower; k<=kUpper; k++){
					
					//skip term if the previous generation cannot have this number of ancestors
					if(toReturn[i-1][k]==Double.NEGATIVE_INFINITY) continue;
					
					double innerSum = Double.NEGATIVE_INFINITY;
					int maleUpper = Math.min(k, j-1);
					for(int m=1; m<=maleUpper; m++){
				
						int f = j-m;
						if(f>k) continue; //skip if number of boxes is greater than number of balls
						
						double nextTerm = logStirling[k][m] + logStirling[k][f] + logFallingFact[m] + logFallingFact[f];
						innerSum = LogSum.addLogSummand(innerSum, nextTerm);
						
					}
					
					double outerSum = -2*k*logN + toReturn[i-1][k] + innerSum;
					
					toReturn[i][j] = LogSum.addLogSummand(toReturn[i][j], outerSum);
					
				}
				
			}			
			
		}

		
		return toReturn;
	}
	
	
	public static double logPg(int N, int n, int g){
		
		if(g==0) return 0;
		
		//precompute ancestor probs
		double[][] ancLogProbN = ancestorLogProbs(N, n-1, g);
		double[][] ancLogProb1 = ancestorLogProbs(N, 1, g);
		double[] logFact = logFactorial(N);
		
		double prevLogP = 0;
		double currLogP = 0;
		double mySum;
		for(int i=1; i<=g; i++){
			
			//boundaries
			int b1Upper = (int) Math.min(N, (n-1)*Math.pow(2, i));
			int b2Upper = (int) Math.min(N, Math.pow(2,i));
			
			mySum = Double.NEGATIVE_INFINITY;
			
			for(int b1=2; b1<=b1Upper; b1++){
				
				double innerSum = Double.NEGATIVE_INFINITY;
				
				for(int b2=2; b2<=b2Upper; b2++){

					double logComb = logFact[N-b1] + logFact[N-b2] - logFact[N-b1-b2] - logFact[N];
					double nextTerm = ancLogProb1[i][b2] + logComb;
					innerSum = LogSum.addLogSummand(innerSum, nextTerm);

				}
				
				double outerSum = ancLogProbN[i][b1] + innerSum;
				
				mySum = LogSum.addLogSummand(mySum, outerSum);

				
			}
			
			
			
			currLogP = prevLogP + mySum;
			prevLogP = currLogP;
			
			
		}
		
		return currLogP;
		
	}
	
	
	public static void main(String[] args){

		int N = 10000;
		int n = 20;
		int g = 4;
	
		
		System.out.println(Math.exp(logPg(N,n,g)));
		
		
	}
	

}

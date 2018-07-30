package test;

import java.io.IOException;
import java.util.Random;



public class testMethods {
	
	
	
	public static void main(String[] args) throws IOException{

		int minN = 50;
		int maxN = 10000;
		int stepSize = 50;
		Random rgen = new Random(123);

		int t = 0;
		int T = 100000;
		int[] count = new int[maxN/stepSize];
		while(t < T) {
			int curr = rgen.nextInt((maxN - minN)/stepSize)*stepSize + minN;
			count[curr/stepSize]++;
			t++;
		}
		
		for(int x : count)
			System.out.println(x);

		

		

		
	}
	
	
}

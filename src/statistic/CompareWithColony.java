package statistic;

import java.io.BufferedReader;
import java.io.IOException;

import utility.DataParser;

public class CompareWithColony {
	
	// go through each posterior sample and compute P(FS | FS), P(HS|HS), P(UR | UR)
	public static int[][] accuracyMatrixTotalProbability(String truePath, String countPath) throws IOException{
		
		
		String truth = DataParser.openReader(truePath).readLine().split("\n")[0];
		
		// (FS, HS, UR, FC, HC) x (FS, HS, UR, FC, HC)
		int[][] acc = new int[5][5];
		
		
		BufferedReader reader = DataParser.openReader(countPath);
		String line;
		
		while((line = reader.readLine()) != null) {
			
			String[] fields = line.split("\\s");
			
			String pred = fields[0];
			int count = Integer.parseInt(fields[1]);
		
			//go through each pair
			for(int i=0; i<truth.length(); i+=3) {
				
				
				int true_idx = 0; //FS
				if(truth.charAt(i) == '1' && truth.charAt(i+1) == '1' && truth.charAt(i+2) == '1') //HS
					true_idx = 1;
				else if(truth.charAt(i) == '0' && truth.charAt(i+1) == '0' && truth.charAt(i+2) == '0') 
					true_idx = 2;
				else if(truth.charAt(i) == '2' && truth.charAt(i+1) == '2' && truth.charAt(i+2) == '2')
					true_idx = 3;
				else if(truth.charAt(i) == '2' && truth.charAt(i+1) == '2' && truth.charAt(i+2) == '1')
					true_idx = 4;
				
				
				if(pred.charAt(i) == '1' && pred.charAt(i+1) == '1' && pred.charAt(i+2) == '2') 
					acc[true_idx][0]+= count;
				else if(pred.charAt(i) == '1' && pred.charAt(i+1) == '1' && pred.charAt(i+2) == '1') 
					acc[true_idx][1]+= count;
				else if(pred.charAt(i) == '0' && pred.charAt(i+1) == '0' && pred.charAt(i+2) == '0') 
					acc[true_idx][2]+= count;
				else if(pred.charAt(i) == '2' && pred.charAt(i+1) == '2' && pred.charAt(i+2) == '2') 
					acc[true_idx][3]+= count;
				else
					acc[true_idx][4]+= count;
				
				
			}
			
			
		}
		
		
		
		return acc;
		
		
		
		
	}	
	

	// go through each posterior sample and compute P(FS | FS), P(HS|HS), P(UR | UR)
	public static double[][] accuracyMatrix(String truePath, String countPath, int n, double thresh) throws IOException{
		
		
		String truth = DataParser.openReader(truePath).readLine().split("\n")[0];

		
		// (FS, HS, UR, FC, HC) x (FS, HS, UR, FC, HC)
		int[][][] counts = new int[n][n][5];
		
		
		BufferedReader reader = DataParser.openReader(countPath);
		String line;
		
		while((line = reader.readLine()) != null) {
			
			String[] fields = line.split("\\s");
			
			String pred = fields[0];
			int count = Integer.parseInt(fields[1]);
		
			
			
			//for each pair (i,j), count FS, HS, UR
			int idx = 0;
			for(int i=0; i<n; i++) {
				for(int j=i+1; j<n; j++) {
					
					//FS
					if(pred.charAt(idx) == '1' && pred.charAt(idx+1) == '1' && pred.charAt(idx+2) == '2') 
						counts[i][j][0] += count;
					//HS
					else if(pred.charAt(idx) == '1' && pred.charAt(idx+1) == '1' && pred.charAt(idx+2) == '1') 
						counts[i][j][1] += count;
					//UR
					else if(pred.charAt(idx) == '0' && pred.charAt(idx+1) == '0' && pred.charAt(idx+2) == '0') 
						counts[i][j][2] += count;
					//FC
					else if(pred.charAt(idx) == '2' && pred.charAt(idx+1) == '2' && pred.charAt(idx+2) == '2') 
						counts[i][j][3] += count;
					//HC
					else
						counts[i][j][4] += count;
				
					
					idx += 3;
					
				}
			}
			
		}
		
		
		
		//assign
		String[][] assignment = new String[n][n];
		for(int i=0; i<n; i++) {
			for(int j=i+1; j<n; j++) {
				
				//find one with the max number of votes
				int max = -1;
				int maxIdx = -1;
				for(int k=0; k<5; k++) {
					//TODO skip unrelated
					//if(k == 2) continue;
					
					if(counts[i][j][k] > max) {
						max = counts[i][j][k];
						maxIdx = k;
					}
				}
				
				
				/*
				//if UR_count * thresh > sum of all other counts, assign as UR
				if(counts[i][j][2]*thresh > counts[i][j][0] + counts[i][j][1] + counts[i][j][3] + counts[i][j][4]) {
					maxIdx = 2;
				}
				*/
				
				
				
				if(maxIdx == 0)
					assignment[i][j] = "FS";
				else if(maxIdx == 1)
					assignment[i][j] = "HS";
				else if(maxIdx == 2)
					assignment[i][j] = "UR";
				else if(maxIdx == 3)
					assignment[i][j] = "FC";
				else
					assignment[i][j] = "HC";
				
				
			}
			
			
		}
		
		
		
		//compute accuracy
		double[][] acc = new double[5][5];
		
		int k = 0;
		for(int i=0; i<n; i++) {
			for(int j=i+1; j<n; j++) {
				
				//get true relationship
				String trueRel = "UR";
				if(truth.charAt(k)=='1' && truth.charAt(k+1)=='1' && truth.charAt(k+2)=='2')
					trueRel = "FS";
				else if(truth.charAt(k)=='1' && truth.charAt(k+1)=='1' && truth.charAt(k+2)=='1')
					trueRel = "HS";
				else if(truth.charAt(k)=='2' && truth.charAt(k+1)=='2' && truth.charAt(k+2)=='2')
					trueRel = "FC";
				else if(truth.charAt(k)=='2' && truth.charAt(k+1)=='2' && truth.charAt(k+2)=='1')
					trueRel = "HC";

				
				k+=3;

				
				String pred = assignment[i][j];
				if(trueRel.equals("FS")) {
					if(pred.equals("FS")) acc[0][0]++;
					else if(pred.equals("HS")) acc[0][1]++;
					else if(pred.equals("UR")) acc[0][2]++;
					else if(pred.equals("FC")) acc[0][3]++;
					else acc[0][4]++;
				}
				else if(trueRel.equals("HS")) {
					if(pred.equals("FS")) acc[1][0]++;
					else if(pred.equals("HS")) acc[1][1]++;
					else if(pred.equals("UR")) acc[1][2]++;
					else if(pred.equals("FC")) acc[1][3]++;
					else acc[1][4]++;
				}
				else if(trueRel.equals("UR")) {
					if(pred.equals("FS")) acc[2][0]++;
					else if(pred.equals("HS")) acc[2][1]++;
					else if(pred.equals("UR")) acc[2][2]++;
					else if(pred.equals("FC")) acc[2][3]++;
					else acc[2][4]++;
				}
				else if(trueRel.equals("FC")) {
					if(pred.equals("FS")) acc[3][0]++;
					else if(pred.equals("HS")) acc[3][1]++;
					else if(pred.equals("UR")) acc[3][2]++;
					else if(pred.equals("FC")) acc[3][3]++;
					else acc[3][4]++;
				}
				else{
					if(pred.equals("FS")) acc[4][0]++;
					else if(pred.equals("HS")) acc[4][1]++;
					else if(pred.equals("UR")) acc[4][2]++;
					else if(pred.equals("FC")) acc[4][3]++;
					else acc[4][4]++;
				}
					
				
				
			}
			
			
		}
		
		
		/*
		//normalize acc
		for(int i=0; i<acc.length; i++) {
			
			double rowSum = 0d;
			
			for(int j=0; j<acc[0].length; j++) {
				
				rowSum += acc[i][j];
				
			}
			
			for(int j=0; j<acc[0].length; j++) {
				
				acc[i][j] /= rowSum;
				
			}
			
		}
		*/
		
		
		return acc;
		
		
		
		
	}	
	
	
}
